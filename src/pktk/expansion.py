"""Dataset expansion and identity validation; no new model training."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from pktk.collectors.http import Collector,log
from pktk.collectors import pmc,pmc_tables,discovery,cebs
from pktk.qc.composite import audit as composite_audit
from pktk.qc.attrition import audit
from pktk.standardization.resolve import registry,resolve
from pktk.standardization.aliases import enrich
from pktk.standardization.core import identity


def write_json(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x),allow_nan=False)+'\n')


def build_series(endpoints):
    grouping=['molecular_group_id','study_id','protocol_id','species','strain','sex','route','matrix','regimen','vehicle','fasting','endpoint_origin']
    grouping += [c for c in ['infusion_duration_min','dose_volume_ml_kg'] if c in endpoints]
    rows=[]
    for vals,g in endpoints.groupby(grouping,dropna=False,sort=True):
        if g.dose_mg_kg.duplicated().any():raise ValueError('Duplicate dose in curated protocol; cannot silently average')
        if g.source_record_id.duplicated().any():raise ValueError('Duplicate source record')
        doses=sorted(float(v) for v in g.dose_mg_kg.unique())
        auc=g[g.high_confidence_auc & g.auc.gt(0)]
        cmax=g[g.high_confidence_cmax & g.cmax.gt(0)]
        auc_ok=len(auc)>=2 and auc.auc_definition.nunique()==1 and auc.auc_unit.nunique()==1
        cmax_ok=len(cmax)>=2 and cmax.cmax_unit.nunique()==1
        if g.publication_doi.nunique()!=1:raise ValueError('Cross-publication pooling')
        r=dict(zip(grouping,vals))
        r.update(series_id=identity(vals),compound_name=g.compound_name.iloc[0],n_doses=len(doses),
                 dose_min=min(doses),dose_max=max(doses),dose_ratio=max(doses)/min(doses),doses_json=json.dumps(doses),
                 n_groups=len(g),auc_available=bool(g.auc.notna().all()),cmax_available=bool(g.cmax.notna().all()),
                 auc_definition=';'.join(sorted(g.auc_definition.unique())),cmax_unit=';'.join(sorted(g.cmax_unit.unique())),
                 auc_unit=';'.join(sorted(g.auc_unit.unique())),eligible_auc=bool(auc_ok),eligible_cmax=bool(cmax_ok),
                 primary_eligible=bool(g.primary_eligible.all() and (auc_ok or cmax_ok) and len(doses)>=2),
                 protocol_verified=bool(g.protocol_verified.all()),has_numeric_curve=False,
                 resolved_structure=bool(g.inchikey.ne('').all()),publication_doi=g.publication_doi.iloc[0],
                 qc_flags=';'.join(sorted(set(';'.join(g.qc_flags).split(';'))- {''})),
                 source_record_ids=json.dumps(g.source_record_id.tolist()),
                 metadata_completeness=sum(bool(g[k].iloc[0]) for k in ['species','strain','sex','route','matrix','protocol_id'])/6)
        rows.append(r)
    return pd.DataFrame(rows)


def run(root,offline=False,resolve_identities=True):
    root=Path(root);collector=Collector(root,offline)
    log(root/'reports/progress.jsonl',phase='expansion',status='started',offline=offline)
    # Rebuild source-derived endpoints from pinned bytes, not from editable processed CSVs.
    config=json.loads((root/'configs/sources.json').read_text())
    original,_=pmc.collect(collector,config['pmc_sha256'])
    original=original.assign(protocol_verified=True,primary_eligible=True,casrn='',pubchem_cid='',infusion_duration_min=0)
    specs=json.loads((root/'configs/curation/pmc_expansion.json').read_text())
    round2_specs=json.loads((root/'configs/curation/pmc_round2.json').read_text())
    new,manifests,failures=pmc_tables.collect(collector,specs)
    if new.empty:raise ValueError('No expansion data collected')
    extra,extra_manifests,extra_failures=pmc_tables.collect(collector,round2_specs)
    cebs_data,cebs_audit=cebs.collect(collector)
    manifests.update(extra_manifests);failures.extend(extra_failures)
    x=pd.concat([original,new,cebs_data,extra],ignore_index=True)
    # Incomplete oral formulation evidence is retained, but withheld from the primary subset.
    unknown=x.qc_flags.str.contains('oral_vehicle_not_explicit')
    x.loc[unknown,['primary_eligible','high_confidence_auc','high_confidence_cmax']]=False
    obs=pd.read_csv(root/'data/processed/observations.csv',keep_default_na=False,low_memory=False)
    attrition=audit(root)
    composite_summary=composite_audit(root,cebs_data)
    entries=enrich(collector,registry(obs,x))
    candidate_table=pd.read_csv(root/'data/processed/candidate_series.csv')
    priority=set(candidate_table.loc[candidate_table.n_doses>=2,'compound_id'])
    if resolve_identities:
        resolved=resolve(collector,entries,priority_ids=priority,limit=config.get('identity_priority_limit',15))
    else:
        resolved=pd.read_csv(root/'data/processed/compound_identity_registry.csv',keep_default_na=False)
        if set(resolved.compound_id)!=set(r['compound_id'] for r in entries):raise ValueError('Identity registry scope mismatch')
    maps=resolved.set_index('compound_id')
    for col in ['canonical_smiles','inchikey','pubchem_cid','identity_status','chebi_id','resolved_cas_json','scaffold']:
        x[col]=x.compound_id.map(maps[col]).fillna('')
    x['molecular_group_id']=x.apply(lambda r:'InChIKey:'+r.inchikey if r.identity_status=='resolved_exact' else r.compound_id,axis=1)
    x['qc_flags']=x.apply(lambda r:';'.join(f for f in r.qc_flags.split(';') if not (f=='missing_structure' and r.inchikey)),axis=1)
    for col in ['auc','cmax','auc_sd','cmax_sd']:
        x[col]=pd.to_numeric(x[col],errors='coerce')
    if (x.dose_mg_kg<=0).any():raise ValueError('Nonpositive curated dose')
    if x.source_record_id.duplicated().any():raise ValueError('Repeated source table cell')
    for endpoint in ['auc','cmax']:
        if (x[endpoint].dropna()<=0).any():raise ValueError('Nonpositive curated exposure')
    series=build_series(x)
    primary=series[series.primary_eligible]
    accepted_ids={rid for ids in primary.source_record_ids for rid in json.loads(ids)}
    high=x[x.source_record_id.isin(accepted_ids)].copy()
    # DOI/identity overlaps are reported, never silently treated as independent replication.
    obs['molecular_group_id']=obs.compound_id.map(maps.inchikey).fillna('')
    obs['molecular_group_id']=obs.apply(lambda r:'InChIKey:'+r.molecular_group_id if r.molecular_group_id else r.compound_id,axis=1)
    overlap=[]
    for _,g in x.groupby(['publication_doi','molecular_group_id']):
        a=obs[((obs.publication_doi.str.strip()==g.publication_doi.iloc[0]) | obs.document_id.isin(g.document_id)) & (obs.molecular_group_id==g.molecular_group_id.iloc[0])]
        if len(a):overlap.append(dict(doi=g.publication_doi.iloc[0],molecular_group_id=g.molecular_group_id.iloc[0],cvtdb_rows=len(a),pmc_groups=len(g),decision='same_publication_not_independent'))
    for name,data in [('reported_endpoints_expanded',x),('crossdose_series_expanded',series),('high_confidence_endpoints',high)]:
        data.to_csv(root/'data/processed'/f'{name}.csv',index=False)
    write_json(root/'reports/expansion_overlap.json',overlap)
    write_json(root/'reports/expansion_source_manifest.json',manifests)
    write_json(root/'reports/expansion_collection_failures.json',failures)
    events=discovery.collect(collector)
    pre=obs[obs.species.ne('human')]
    # Source identifiers remain visible; molecular counts are conditional on incomplete resolution.
    identities=set(pre.compound_id)|set(x.compound_id)
    parent_ids=set(pre[pre.parent_analyte].compound_id)|set(x.compound_id)
    def document_key(r):
        return 'doi:'+r.publication_doi.strip().lower() if r.publication_doi.strip() else r.document_id
    document_alias={r.document_id:document_key(r) for _,r in x.iterrows()}
    documents={document_alias.get(r.document_id,document_key(r)) for _,r in pre.iterrows()}|set(x.apply(document_key,axis=1))
    source_documents=set(pre.document_id)|set(x.document_id)
    s=dict(cvtdb_original_observations=len(obs),preclinical_observations=len(pre),
           preclinical_analyte_identifiers_collected=len(identities),preclinical_parent_compound_identifiers=len(parent_ids),
           preclinical_document_studies=len(documents),preclinical_source_document_ids=len(source_documents),
           cvtdb_source_series_ids=obs.experiment_id.nunique(),cvtdb_usable_source_curves=attrition['strict_profiles'],
           cvtdb_candidate_crossdose_series=attrition['strict_crossdose_series'],
           collected_pmc_publications=x.loc[x.source.eq('pmc'),'publication_doi'].nunique(),curated_endpoint_groups=len(x),
           protocol_verified_series=len(series),primary_comparable_series=len(primary),
           primary_comparable_compound_identifiers=primary.molecular_group_id.nunique(),
           primary_independent_publications=primary.loc[~primary.study_id.str.startswith('CEBS:'),'publication_doi'].nunique(),
           primary_independent_studies=primary.study_id.nunique(),
           compounds_at_least_two_comparable_doses=primary.molecular_group_id.nunique(),
           compounds_at_least_three_comparable_doses=primary.loc[primary.n_doses>=3,'molecular_group_id'].nunique(),
           primary_series_at_least_three_doses=int((primary.n_doses>=3).sum()),
           primary_series_with_auc_and_cmax=int((primary.auc_available & primary.cmax_available).sum()),
           primary_auc_eligible_series=int(primary.eligible_auc.sum()),primary_cmax_eligible_series=int(primary.eligible_cmax.sum()),
           primary_both_endpoints_eligible_series=int((primary.eligible_auc & primary.eligible_cmax).sum()),
           primary_endpoint_groups=len(high),flagged_protocol_series=int((~series.primary_eligible).sum()),
           identity_registry_entries=len(resolved),identity_status_counts=resolved.identity_status.value_counts().to_dict(),
           resolved_inchikeys=resolved.loc[resolved.identity_status=='resolved_exact','inchikey'].nunique(),
           resolved_pubchem_cids=int(resolved.pubchem_cid.ne('').sum()),
           primary_compounds_with_resolved_structures=primary.loc[primary.resolved_structure,'molecular_group_id'].nunique(),
           primary_species=sorted(primary.species.unique()),primary_routes=sorted(primary.route.unique()),
           collected_species=sorted(set(pre.species)|set(x.species)),collected_routes=sorted(set(pre.route)|set(x.route)),
           fraction_curated_pmc_groups_primary=len(high[high.source.eq('pmc')])/len(x[x.source.eq('pmc')]),fraction_curated_groups_primary=len(high)/len(x),fraction_parent_identifiers_verified_comparable=primary.molecular_group_id.nunique()/len(parent_ids),
           candidate_counts_not_verified=True,identity_counts_not_fully_structure_deduplicated=True,
           extractions_failed=len(failures),source_overlap_groups=len(overlap),complex_models_trained=0)
    s['reviewed_composite_candidates']=composite_summary['reviewed_candidates']
    s['validated_composite_series']=composite_summary['status_counts'].get('validated_group_endpoints',0)
    s['primary_series_by_source']={source:len(build_series(high[high.source.eq(source)])) for source in sorted(high.source.unique())}
    s['resolved_unique_parent_structures']=resolved.loc[resolved.identity_status.eq('resolved_exact'),'inchikey'].nunique()
    s['unresolved_parent_identifiers']=int(resolved.identity_status.ne('resolved_exact').sum())
    s['primary_noncommercial_series']=int(high.loc[high.qc_flags.str.contains('license_noncommercial'),'protocol_id'].nunique())
    snapshots=[]
    for stage,subset in [('identity_resolution',x[~x.source.eq('cebs') & ~x.publication_doi.isin(extra.publication_doi)]),
                         ('composite_audit',x[~x.publication_doi.isin(extra.publication_doi)]),('additional_publications',x)]:
        st=build_series(subset);p=st[st.primary_eligible]
        scope_ids=set(pre.loc[pre.parent_analyte,'compound_id'])|set(subset.compound_id)
        scoped=resolved[resolved.compound_id.isin(scope_ids)]
        snapshots.append(dict(stage=stage,collected_parent_identifiers=len(scoped),resolved_unique_structures=scoped.loc[scoped.identity_status.eq('resolved_exact'),'inchikey'].nunique(),
            unresolved_parent_identifiers=int(scoped.identity_status.ne('resolved_exact').sum()),validated_compounds=p.molecular_group_id.nunique(),
            independent_studies=p.study_id.nunique(),validated_series=len(p),auc_series=int(p.eligible_auc.sum()),
            cmax_series=int(p.eligible_cmax.sum()),both_eligible_series=int((p.eligible_auc&p.eligible_cmax).sum())))
    write_json(root/'reports/stage_snapshots.json',snapshots)
    write_json(root/'reports/expansion_summary.json',s)
    write_json(root/'reports/summary.json',s)
    report(root,s,series,attrition,resolved,events)
    log(root/'reports/progress.jsonl',phase='expansion',status='succeeded',summary=s)
    print(json.dumps(s,indent=2,default=int))
    return s


def report(root,s,series,attrition,resolved,events):
    from pktk.reporting import expanded_report
    expanded_report(root,s,series,attrition,resolved)


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path('.'));p.add_argument('--offline',action='store_true')
    p.add_argument('--reuse-identities',action='store_true',help='Reuse audited local identity registry instead of resolving again')
    a=p.parse_args();run(a.root,a.offline,not a.reuse_identities)

if __name__=='__main__':main()
