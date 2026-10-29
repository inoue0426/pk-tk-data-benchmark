"""Dataset A: auditable molecular PK endpoints. Never promotes records to Dataset B."""
import argparse,json,math,re,hashlib
from pathlib import Path
from collections import Counter
import pandas as pd
import numpy as np
from pktk.collectors.http import Collector,sha256,log
from pktk.collectors.molecular import collect
from pktk.standardization.molecular import structure,context,endpoint

ENDPOINTS={'CL':'clearance','Vd':'volume_distribution','T1/2':'half_life','AUC':'auc','Cmax':'cmax'}

def dump(path,obj):
    path.write_text(json.dumps(obj,indent=2,allow_nan=False,default=lambda x:x.item() if hasattr(x,'item') else str(x))+'\n')

def publication(doc):
    doi=(doc.get('doi') or '').strip().lower().removeprefix('https://doi.org/')
    if doi:return 'doi:'+doi
    if doc.get('pubmed_id'):return 'pmid:'+str(doc['pubmed_id'])
    return 'chembl_document:'+doc['document_chembl_id']

def chembl_records(rows,metadata,manifest):
    files={}
    for m in manifest['downloads']:
        if 'activity_' in m['file']:
            page=json.loads((manifest['_root']/'data/raw'/m['file']).read_text())
            for r in page['activities']:files[r['activity_id']]=m
    out=[]
    for r in rows:
        a=metadata['assay'][r['assay_chembl_id']];doc=metadata['document'][r['document_chembl_id']]
        description=a.get('description') or '';cx=context(description);st=structure(r.get('canonical_smiles') or '')
        normalized,unit_flag=endpoint(r['standard_type'],r['standard_units'],r['standard_value'],st.get('molecular_weight'))
        exclusions=[];flags=cx.pop('context_flags')
        if 'multiple_routes' in flags or 'multiple_matrices' in flags:exclusions.append('ambiguous_route_or_matrix')
        if r.get('target_organism')!='Rattus norvegicus' or a.get('assay_organism')!='Rattus norvegicus':exclusions.append('species_mismatch_or_missing')
        if r['assay_type']!='A' or a.get('assay_type')!='A':exclusions.append('not_adme_assay')
        if not cx['in_vivo_evidence']:exclusions.append('no_explicit_in_vivo_evidence_or_non_pk_context')
        if st['structure_status']!='resolved':exclusions.append(st['structure_status'])
        if normalized is None:exclusions.append(unit_flag)
        elif unit_flag:flags.append(unit_flag)
        if r.get('standard_relation')!='=' or r.get('standard_upper_value') is not None:exclusions.append('censored_or_interval_value')
        if r.get('data_validity_comment'):exclusions.append('source_data_validity_flag')
        if r.get('potential_duplicate'):flags.append('source_potential_duplicate')
        # A ChEMBL type alone cannot turn a timed concentration or PD response into Cmax/AUC.
        tokens={'Cmax':r'c\s*max|maximum|\bpeak\b','AUC':r'\bauc\b|area under','CL':r'clearance|\bcl\b',
                'Vd':r'volume|\bvd\b|\bvss\b','T1/2':r'half[- ]?life|t\s*1\s*/\s*2|t1/2'}
        if not re.search(tokens[r['standard_type']],description,re.I):exclusions.append('endpoint_description_not_confirmed')
        if cx['matrix'] in ['brain','csf','liver','urine']:flags.append('non_systemic_matrix')
        if r['standard_type']=='AUC' and cx['auc_definition']=='AUC_unspecified':flags.append('auc_definition_unspecified')
        if r['standard_type']=='CL' and cx['route'] not in ['intravenous','']:flags.append('apparent_clearance_possible')
        if r['standard_type']=='Vd':flags.append('volume_definition_unspecified')
        if r['standard_type']=='T1/2':flags.append('half_life_phase_unspecified')
        # Source units may carry normalization omitted by ChEMBL's standard unit field.
        raw_unit=str(r.get('units') or '')
        if re.search(r'/(?:mg|ug|μg)|(?:mg|ug)\^-?1',raw_unit,re.I) and r['standard_type'] in ['AUC','Cmax']:
            cx['dose_normalized']='yes';cx['dose_mg_kg']=None;flags.append('source_unit_dose_normalization')
        meta=files[r['activity_id']]
        record=dict(source='chembl',source_record_id='chembl:'+str(r['activity_id']),compound_id=r['molecule_chembl_id'],
          compound_name=r.get('molecule_pref_name') or '',species='rat',assay_id=r['assay_chembl_id'],
          publication_id=publication(doc),publication_doi=doc.get('doi') or '',pubmed_id=doc.get('pubmed_id') or '',
          original_study_id='',study_status='publication_proxy_experiment_not_verified',
          record_kind='literature_abstracted_endpoint',endpoint=ENDPOINTS[r['standard_type']],value=None,unit='',
          endpoint_definition=cx['auc_definition'] if r['standard_type']=='AUC' else r['standard_type']+'_source_defined',
          original_value=r.get('value'),original_unit=r.get('units'),source_standard_value=r.get('standard_value'),
          source_standard_unit=r.get('standard_units'),source_relation=r.get('relation'),
          source_url=meta['url'],raw_file=meta['file'],raw_sha256=meta['sha256'],retrieved_at=meta['retrieved_at'],
          dataset_version='ChEMBL_37',license='CC-BY-SA-3.0',description=description,
          assay_parameters_json=json.dumps(a.get('assay_parameters',[]),sort_keys=True),
          original_json=json.dumps(r,sort_keys=True),eligible_A=not exclusions,strict_B_eligible=False,
          exclusion_reasons=';'.join(sorted(set(exclusions))),qc_flags=';'.join(sorted(set(flags))),
          **cx,**st)
        if normalized:record.update(normalized)
        out.append(record)
    return out

def pksmart_records(root,offline):
    cfg=json.loads((root/'configs/molecular/sources.json').read_text())['pksmart'];c=Collector(root,offline)
    p,m=c.fetch(cfg['file'],cfg['url'],cfg['version'],cfg['sha256'])
    if hashlib.md5(p.read_bytes()).hexdigest()!=cfg['publisher_file_md5']:raise ValueError('Publisher supplementary checksum mismatch')
    x=pd.read_excel(p,sheet_name=cfg['sheet'],header=1,keep_default_na=False)
    if len(x)!=371 or 'rat_CL_mL_min_kg' not in x:raise ValueError('Unexpected PKSmart animal table schema')
    out=[]
    for i,r in x.iterrows():
        st=structure(r.smiles)
        for species in ['rat','dog','monkey']:
            for suffix,kind,unit in [('CL_mL_min_kg','CL','mL.min-1.kg-1'),('VDss_L_kg','Vd','L.kg-1')]:
                col=species+'_'+suffix;value=r[col]
                if value=='':continue
                norm,reason=endpoint(kind,unit,value)
                exclusions=[]
                if st['structure_status']!='resolved':exclusions.append(st['structure_status'])
                if norm is None:exclusions.append(reason)
                record=dict(source='pksmart',source_record_id=f'pksmart:S2:row{i+3}:{col}',compound_id=f'pksmart:S2:row{i+3}',
                  compound_name=r['name'],species=species,route='intravenous',matrix='',sex='',dose_mg_kg=None,
                  dose_normalized='not_reported',endpoint=ENDPOINTS[kind],value=None,unit='',
                  endpoint_definition='Vss' if kind=='Vd' else 'systemic_clearance',
                  original_value=value,original_unit=unit,source_standard_value=value,source_standard_unit=unit,source_relation='=',
                  publication_id='',publication_doi='',pubmed_id='',original_study_id='',
                  acquisition_publication='doi:10.1186/s13321-025-01066-5',study_status='original_publications_unmapped',
                  record_kind='source_aggregate_endpoint',description='Table S2 animal measured PK; raw-scale numeric columns. Not predicted animal values in Table S6.',
                  source_url=m['url'],raw_file=cfg['file'],raw_sha256=m['sha256'],retrieved_at=m['retrieved_at'],
                  dataset_version=cfg['version'],license=cfg['license'],original_json=json.dumps(r.to_dict(),sort_keys=True),
                  eligible_A=not exclusions,strict_B_eligible=False,exclusion_reasons=';'.join(exclusions),
                  qc_flags='aggregate_no_dose_or_original_study;missing_matrix;source_median_aggregation',**st)
                if norm:record.update(norm)
                out.append(record)
    return out,m

def benchmark_records(root):
    from pktk.benchmark.episodes import endpoint_table
    ep=pd.read_csv(root/'data/processed/high_confidence_endpoints.csv',keep_default_na=False)
    series=pd.read_csv(root/'data/processed/crossdose_series_expanded.csv',keep_default_na=False)
    table=endpoint_table(ep,series);original=ep.set_index('source_record_id');out=[]
    for _,r in table.iterrows():
        src=original.loc[r.source_record_id];st=structure(src.canonical_smiles)
        record=dict(source='strict_b',source_record_id='strict_b:'+r.source_record_id+':'+r.endpoint,
          compound_id=src.compound_id,compound_name=src.compound_name,species=src.species,sex=src.sex,
          route=src.route,matrix=src.matrix,dose_mg_kg=src.dose_mg_kg,dose_normalized='no',
          endpoint=r.endpoint,value=r.exposure,unit=r.exposure_unit,endpoint_definition=r.exposure_definition,
          publication_id=('doi:'+src.publication_doi.lower()) if src.publication_doi else src.study_id,
          publication_doi=src.publication_doi,pubmed_id=src.get('pmid',''),original_study_id=src.study_id,
          study_status='original_protocol_verified',record_kind='reviewed_group_endpoint',
          original_value=r.exposure,original_unit=r.exposure_unit,source_standard_value=r.exposure,source_standard_unit=r.exposure_unit,
          source_relation='=',source_url=src.source_url,raw_file='',raw_sha256=src.raw_sha256,retrieved_at=src.retrieved_at,
          dataset_version=src.dataset_version,license=src.get('license_url','') or 'see strict benchmark source inventory',
          original_json=src.original_json,eligible_A=st['structure_status']=='resolved',strict_B_eligible=True,
          exclusion_reasons='' if st['structure_status']=='resolved' else st['structure_status'],qc_flags=src.qc_flags,**st)
        out.append(record)
    return out

def run(root=Path('.'),offline=False):
    root=Path(root);out=root/'data/processed/molecular';out.mkdir(parents=True,exist_ok=True)
    log(root/'reports/progress.jsonl',phase='molecular_pk',status='started',offline=offline)
    rows,meta,manifest=collect(root,offline);manifest['_root']=root
    records=chembl_records(rows,meta,manifest);pksmart,pmeta=pksmart_records(root,offline);records+=pksmart+benchmark_records(root)
    x=pd.DataFrame(records).fillna('')
    if x.source_record_id.duplicated().any():raise ValueError('Duplicate source endpoint IDs')
    # Detect one source molecule ID associated with incompatible standardized structures.
    bad=x[x.structure_status.eq('resolved')].groupby(['source','compound_id']).inchikey.nunique()
    conflicts=set(bad[bad>1].index)
    for i,r in x.iterrows():
        if (r.source,r.compound_id) in conflicts:
            x.loc[i,'eligible_A']=False;x.loc[i,'exclusion_reasons']+=';conflicting_structures_for_identifier'
    good=x[x.eligible_A].copy()
    # Overlap is a review flag, never evidence of an independent new experiment.
    keys=['inchikey','publication_id','species','route','matrix','endpoint','endpoint_definition','dose_mg_kg','dose_normalized','value','unit']
    overlap=good[good.publication_id.ne('') & good.duplicated(keys,keep=False)].copy()
    overlap.to_csv(out/'potential_endpoint_overlap.csv',index=False)
    x.to_csv(out/'endpoint_records.csv',index=False,float_format='%.12g')
    good.to_csv(out/'eligible_endpoints.csv',index=False,float_format='%.12g')
    registry=[]
    for key,g in good.groupby('inchikey',sort=True):
        registry.append(dict(inchikey=key,canonical_smiles=g.canonical_smiles.iloc[0],connectivity_key=g.connectivity_key.iloc[0],
            scaffold=g.scaffold.iloc[0],sources_json=json.dumps(sorted(g.source.unique())),
            source_compound_ids_json=json.dumps(sorted(set(g.source+':'+g.compound_id))),
            names_json=json.dumps(sorted(set(g.compound_name)-{''})),endpoints_json=json.dumps(sorted(g.endpoint.unique())),
            species_json=json.dumps(sorted(g.species.unique())),records=len(g)))
    pd.DataFrame(registry).to_csv(out/'compounds.csv',index=False)
    # Dose candidates are screening leads only, never protocol-validated B additions.
    q=good[(good.source=='chembl') & good.endpoint.isin(['auc','cmax']) & good.dose_mg_kg.ne('') & good.route.ne('') & good.matrix.ne('') & good.dose_normalized.ne('yes')]
    candidates=[]
    for key,g in q.groupby(['inchikey','publication_id','species','sex','route','matrix','endpoint','endpoint_definition'],sort=True):
        doses=sorted(set(float(v) for v in g.dose_mg_kg))
        if len(doses)>=2:candidates.append(dict(zip(['inchikey','publication_id','species','sex','route','matrix','endpoint','endpoint_definition'],key))|dict(doses_json=json.dumps(doses),n_doses=len(doses),source_records_json=json.dumps(sorted(g.source_record_id)),status='unverified_original_protocol_required',strict_B_eligible=False))
    pd.DataFrame(candidates,columns=['inchikey','publication_id','species','sex','route','matrix','endpoint','endpoint_definition','doses_json','n_doses','source_records_json','status','strict_B_eligible']).to_csv(out/'dose_candidates_unverified.csv',index=False)
    from pktk.benchmark.molecular_feasibility import assess
    feasibility=assess(good,out)
    from pktk.molecular_reporting import reports
    summary=reports(root,x,good,meta,manifest,pmeta,feasibility,candidates,conflicts)
    log(root/'reports/progress.jsonl',phase='molecular_pk',status='succeeded',**summary['dataset_a'])
    print(json.dumps(summary['dataset_a'],indent=2))

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path('.'));p.add_argument('--offline',action='store_true');a=p.parse_args()
    try:run(a.root,a.offline)
    except Exception as exc:
        log(a.root/'reports/progress.jsonl',phase='molecular_pk',status='failed',error=str(exc))
        raise

if __name__=='__main__':main()
