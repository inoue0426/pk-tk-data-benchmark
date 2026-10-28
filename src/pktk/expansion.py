"""Dataset expansion and identity validation; no new model training."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from pktk.collectors.http import Collector,log
from pktk.collectors import pmc,pmc_tables,discovery
from pktk.qc.attrition import audit
from pktk.standardization.resolve import registry,resolve
from pktk.standardization.core import identity


def write_json(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x),allow_nan=False)+'\n')


def build_series(endpoints):
    grouping=['molecular_group_id','study_id','protocol_id','species','strain','sex','route','matrix','regimen','vehicle','fasting','endpoint_origin']
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
    new,manifests,failures=pmc_tables.collect(collector,specs)
    if new.empty:raise ValueError('No expansion data collected')
    x=pd.concat([original,new],ignore_index=True)
    # Incomplete oral formulation evidence is retained, but withheld from the primary subset.
    unknown=x.qc_flags.str.contains('oral_vehicle_not_explicit')
    x.loc[unknown,['primary_eligible','high_confidence_auc','high_confidence_cmax']]=False
    obs=pd.read_csv(root/'data/processed/observations.csv',keep_default_na=False,low_memory=False)
    attrition=audit(root)
    entries=registry(obs,x)
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
        a=obs[(obs.publication_doi.str.strip()==g.publication_doi.iloc[0]) & (obs.molecular_group_id==g.molecular_group_id.iloc[0])]
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
    documents=set(pre.apply(document_key,axis=1))|set(x.apply(document_key,axis=1))
    source_documents=set(pre.document_id)|set(x.document_id)
    s=dict(cvtdb_original_observations=len(obs),preclinical_observations=len(pre),
           preclinical_analyte_identifiers_collected=len(identities),preclinical_parent_compound_identifiers=len(parent_ids),
           preclinical_document_studies=len(documents),preclinical_source_document_ids=len(source_documents),
           cvtdb_source_series_ids=obs.experiment_id.nunique(),cvtdb_usable_source_curves=attrition['strict_profiles'],
           cvtdb_candidate_crossdose_series=attrition['strict_crossdose_series'],
           collected_pmc_publications=x.publication_doi.nunique(),curated_endpoint_groups=len(x),
           protocol_verified_series=len(series),primary_comparable_series=len(primary),
           primary_comparable_compound_identifiers=primary.molecular_group_id.nunique(),
           primary_independent_publications=primary.publication_doi.nunique(),
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
           fraction_curated_pmc_groups_primary=len(high)/len(x),fraction_parent_identifiers_verified_comparable=primary.molecular_group_id.nunique()/len(parent_ids),
           candidate_counts_not_verified=True,identity_counts_not_fully_structure_deduplicated=True,
           extractions_failed=len(failures),source_overlap_groups=len(overlap),complex_models_trained=0)
    write_json(root/'reports/expansion_summary.json',s)
    write_json(root/'reports/summary.json',s)
    report(root,s,series,attrition,resolved,events)
    log(root/'reports/progress.jsonl',phase='expansion',status='succeeded',summary=s)
    print(json.dumps(s,indent=2,default=int))
    return s


def report(root,s,series,attrition,resolved,events):
    primary=series[series.primary_eligible]
    lines=['| Compound | Route | Doses | Eligible endpoints | Study |','|---|---|---:|---|---|']
    for _,r in primary.iterrows():
        endpoints=','.join(n for n,b in [('AUC',r.eligible_auc),('Cmax',r.eligible_cmax)] if b)
        lines.append(f'| {r.compound_name} ({r.sex}) | {r.route} | {r.n_doses} | {endpoints} | {r.publication_doi} |')
    text=f'''# Verified expanded dataset — 2026-10-28

The curated primary subset grew from **2 to {s['primary_comparable_series']} comparable multi-dose series**, from **1 to {s['primary_comparable_compound_identifiers']} compound identifiers**, and from **1 to {s['primary_independent_publications']} independent publications**. Series within the same publication are correlated strata, not independent studies. This is a validated dataset expansion, not a claim of adequate ML diversity.

{chr(10).join(lines)}

All counts are computed from downloaded bytes. The pipeline retains {s['cvtdb_original_observations']:,} CvTdb observations ({s['preclinical_observations']:,} preclinical) and {s['curated_endpoint_groups']} curated PMC dose groups from {s['collected_pmc_publications']} publications. Across sources there are {s['preclinical_parent_compound_identifiers']} preclinical parent-compound identifiers, {s['preclinical_analyte_identifiers_collected']} analyte identifiers including metabolites, and {s['preclinical_document_studies']} document-level study proxies after exact DOI normalization. These are not fully resolved unique molecules or proven independent experiments.

- Compounds with ≥2 comparable doses: **{s['compounds_at_least_two_comparable_doses']}**; with ≥3: **{s['compounds_at_least_three_comparable_doses']}**.
- Primary series with reported AUC and Cmax: **{s['primary_series_with_auc_and_cmax']}**. AUC eligible: **{s['primary_auc_eligible_series']}**; Cmax eligible: **{s['primary_cmax_eligible_series']}**; both eligible: **{s['primary_both_endpoints_eligible_series']}**. AUC definitions are kept separate; no assumed AUC0-t windows.
- Primary dose groups: **{s['primary_endpoint_groups']}/{s['curated_endpoint_groups']} ({100*s['fraction_curated_pmc_groups_primary']:.1f}%)** of curated PMC groups. Verified comparable compound fraction: **{100*s['fraction_parent_identifiers_verified_comparable']:.2f}%** of collected parent identifiers. No mixed timepoint/group denominator is presented as an experiment fraction.
- Collected species: {', '.join(s['collected_species'])}; routes: {', '.join(s['collected_routes'])}. Primary subset: {', '.join(s['primary_species'])}, {', '.join(s['primary_routes'])}, plasma.

## Why CvTdb yielded only ten candidates

{attrition['singleton_series']:,} of {attrition['source_series_ids']:,} source series IDs have one observation, often individual animals in NTP datasets. A three-timepoint-per-ID requirement misses potentially reconstructable composite designs. ND/NQ omissions recover only {attrition['censor_tolerant_crossdose_series']} candidate groups versus ten under the original gate. There are {attrition['potential_crossdose_groups_if_composite_sampling_confirmed']} possible groups if composite sampling is established; none is promoted without source protocol review. Relaxed metadata grouping yields {attrition['strict_profiles_less_restrictive_metadata_groups']} groups rather than ten (some groups merge), so merely loosening grouping does not solve the problem. The parser did not discard observations. See [quantitative audit](cvtdb_attrition.md).

## Identity validation and collection

All {s['identity_registry_entries']} parent-identifier registry entries retain their source names/CAS values. {s['resolved_inchikeys']} distinct structures resolve through ChEBI exact CAS/name evidence and RDKit agreement with the full source InChIKey; {s['resolved_pubchem_cids']} entries have a single PubChem CID cross-reference through UniChem. Multiple structures, CAS/name conflicts, salts, stereoisomers and failed lookups remain explicit alternatives. {s['primary_compounds_with_resolved_structures']} primary compounds have resolved structures. No name-only fuzzy match or automatic parent/salt collapse is accepted.

The DA-6886 supplementary ZIP was downloaded and its PDF inspected: calibration/stability controls and chromatograms, not extra in vivo dose groups. The legacy PK-DB statistics and single-study metadata endpoints were accessible after reducing page size; the inspected record is human and closed-licence and was not imported. New PK-DB endpoints had TLS failures. PubChem's direct REST path and ChEBI's FTP listing were disallowed by robots; public ChEBI/UniChem APIs were used with their own documented terms. No restrictions were bypassed.

## Remaining limitations and benchmark feasibility

{ s['flagged_protocol_series'] } additional within-publication series remain outside the primary subset: oral vehicle composition changes, unreported oral vehicle, dose-volume changes/composite sampling, or suspect units. The original HY-071085 AUC terminal time and dog table remain unresolved. Some primary IV studies lack fasting details and AUC tail fractions; these are flagged, not filled in. Primary evidence includes group summaries and no new numerical individual curves.

No complex models were trained. The seven primary compounds and five publications are too few for reliable compound/scaffold-disjoint evaluation or tuning. Existing baseline results concern only the original HY-071085 pilot; they are not expanded-dataset results. Task A is structurally feasible for eligible endpoints, but no additional model development was prioritized. Tasks B–D remain unsupported for meaningful generalization/classification. Identity and protocol curation, especially NTP composite profiles and new independent publications, remain the next priorities.

See [source inventory](data_inventory.md), [QC](qc_report.md), [validation](validation.md), [machine-readable statistics](expansion_summary.json), [identity provenance](identity_source_manifest.json), and [unresolved issues](unresolved_issues.json).
'''
    (root/'reports/final_report.md').write_text(text)
    (root/'reports/collection_summary.md').write_text('# Expanded collection summary\n\nCurrent authoritative counts are in [the final report](final_report.md) and [summary.json](summary.json).\n\n```json\n'+json.dumps(s,indent=2)+'\n```\n')
    qc_path=root/'reports/qc_counts.json'
    qc=json.loads(qc_path.read_text())
    qc['expansion']={k:s[k] for k in ['curated_endpoint_groups','primary_endpoint_groups','primary_comparable_series','flagged_protocol_series','extractions_failed','source_overlap_groups','identity_status_counts']}
    write_json(qc_path,qc)
    (root/'reports/qc_report.md').write_text('# Current QC report\n\nAll original observations and all selected curated endpoints are retained. Exclusions apply to primary eligibility, not source deletion.\n\n```json\n'+json.dumps(qc,indent=2)+'\n```\n\nArticle/table hashes, exact cell addresses, dose columns, evidence phrases and units are checked against reviewed configurations. Censored CvTdb values remain censored; audit-only recovery does not imply validated AUC. Primary grouping requires publication, protocol, compound, species, strain, sex, route, matrix, regimen and formulation conditions. Source DOI overlaps are flagged; no records from independent publications are pooled. Chemical resolution accepts only unique full-InChIKey agreement; conflicts retain candidate lists. Exact CAS checksum is validated.\n\nFlags remain for unavailable fasting, AUC tail fractions, aggregate-only uncertainty and missing structures. Reported AUC0-inf is distinct from observed-window AUC. Rat HY-071085 AUC0-tn is ineligible. Oral formulation changes and ambiguous Cmax units in PMC6160927, incomplete tofacitinib oral formulation evidence, and memantine composite/dose-volume conditions are outside the primary subset.\n')
    (root/'reports/benchmark_feasibility.md').write_text('# Expanded benchmark feasibility\n\n'+
        f"{s['primary_comparable_series']} primary series span {s['primary_comparable_compound_identifiers']} compounds and {s['primary_independent_publications']} publications. {s['primary_auc_eligible_series']} series support comparable AUC definitions and {s['primary_cmax_eligible_series']} support Cmax.\n\n"+
        'No new models were fitted. Original results in results/predictions.csv and metrics.csv apply only to HY-071085; they do not characterize this expanded dataset.\n\nTask A: eligible endpoint-specific low-to-high dose episodes can be defined; keep complete study/route/sex/matrix groups and never mix AUC windows. Task B: too few independent compounds/publications for reliable evaluation. Task C: incomplete structures and too few scaffold families. Task D: no independent nonlinear labels/negative controls. No compound/scaffold split or complex ML model is claimed. Count shared-publication compounds and sex/route strata as correlated when evaluating future models.\n')


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path('.'));p.add_argument('--offline',action='store_true')
    p.add_argument('--reuse-identities',action='store_true',help='Reuse audited local identity registry instead of resolving again')
    a=p.parse_args();run(a.root,a.offline,not a.reuse_identities)

if __name__=='__main__':main()
