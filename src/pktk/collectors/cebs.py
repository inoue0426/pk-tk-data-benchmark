"""Reviewed CEBS group endpoints; composite sampling is never individual-animal PK."""
import json,re
from pathlib import Path
import pandas as pd
from pktk.standardization.core import identity


def collect(collector):
    config=json.loads((collector.root/'configs/curation/cebs.json').read_text())
    records=[];audit=[]
    for spec in config:
        path,meta=collector.fetch(spec['file'],spec['url'],'CEBS-audit-2026-10-28',spec['sha256'])
        x=pd.read_excel(path,keep_default_na=False)
        for evidence in spec['evidence']:
            if not x.astype(str).apply(lambda c:c.str.contains(evidence,regex=False)).any().any():
                raise ValueError('Missing CEBS protocol evidence: '+evidence)
        for index,r in x.iterrows():
            if r['ROUTE'] not in ['Intravenous','Oral gavage','IV','Gavage']:continue
            if spec['study']=='S0592' and (r['ROUTE']!='Oral gavage' or r['Best Fit Model']!='noncompartmental model'):continue
            if spec['study']=='K10262' and r['ROUTE']!='Gavage':continue
            dose=float(r.DOSE);flags=['composite_group_endpoint','fasting_not_reported','auc_tail_fraction_unavailable']
            eligible=True
            if spec['study']=='S0592':
                m=re.fullmatch(r'Actual dose is ([0-9.]+) mg/kg\.?',str(r.COMMENT).strip())
                if not m:raise ValueError('Unreviewed actual dose comment')
                dose=float(m[1]);flags+=['actual_not_nominal_dose','source_summary_metadata_errors']
            if spec['study']=='K10262' and r.SPECIES_COMMON_NAME=='Rat':
                eligible=False;flags+=['rat_middle_dose_conflict_35_vs_70']
            if r['DOSE_UNIT']!='mg/kg' or r['TISSUE']!='Plasma' or r['Is_Parent_Compound']!='T':
                raise ValueError('Unreviewed CEBS dose, matrix or analyte')
            auc_column='AUC_0-T' if spec['study']=='K10262' else 'AUCinf_pred'
            raw_auc=r[auc_column];raw_cmax=r['Cmax_obs']
            factor=1/60 if spec['study'] in ['S0592','K10262'] else 1
            expected={'S0592':'ug*min/mL','S0613':'ug*hr/mL','K10262':'ug mL^-1 min'}[spec['study']]
            if raw_auc!='' and r[auc_column+'[DEPOSITOR_UNIT]']!=expected:raise ValueError('Unreviewed AUC unit')
            if raw_cmax!='' and r['Cmax_obs[DEPOSITOR_UNIT]']!='ug/mL':raise ValueError('Unreviewed Cmax unit')
            if spec['study']=='K10262':flags+=['auc_terminal_time_unresolved']
            if r['Cmax_pred']!='':flags+=['predicted_cmax_not_observed']
            if r['EXCEPTION_COMMENTS']:flags+=['source_exception_preserved']
            rid=f'CEBS:{spec["study"]}:Sheet1:row{index+2}'
            species=r.SPECIES_COMMON_NAME.lower();route='oral' if r.ROUTE in ['Oral gavage','Gavage'] else 'intravenous'
            protocol=f'CEBS:{spec["study"]}:{species}:{route}:single'
            records.append(dict(source='cebs',source_record_id=rid,experiment_id=rid,compound_id=spec['compound_id'],
                compound_name=r.CHEMICAL,casrn=r.CAS_NO,pubchem_cid='',canonical_smiles='',inchikey='',
                species=species,strain=r.STRAIN,sex=r.SEX.lower(),route=route,matrix='plasma',dose_mg_kg=dose,
                dose_unit='mg/kg',dose_basis='analytically measured actual' if spec['study']=='S0592' else 'reported nominal',
                dose_volume_ml_kg=float(r.VOLUME_PER_ADMIN),auc=float(raw_auc)*factor if raw_auc!='' else None,
                cmax=float(raw_cmax) if raw_cmax!='' else None,auc_sd=None,cmax_sd=None,
                auc_definition='AUC0-t_unspecified' if spec['study']=='K10262' else 'AUC0-inf',auc_unit='ug*h/mL',cmax_unit='ug/mL',
                n_animals=r.SUBJECTS_PER_GROUP,n_animals_per_timepoint=3,regimen='once',vehicle=r.VEHICLE_DESC,
                fasting='not reported',infusion_duration_min=None,protocol_id=protocol,study_id='CEBS:'+spec['study'],
                document_id=spec['cvtdb_document_id'],publication_doi=spec['doi'],pmid='',publication_url=spec['landing_url'],
                source_url=meta['url'],retrieved_at=meta['retrieved_at'],dataset_version=meta['version'],raw_sha256=meta['sha256'],
                original_json=json.dumps(r.to_dict(),default=str),protocol_evidence_json=json.dumps(spec['evidence']),
                endpoint_origin='reported_composite_profile',protocol_verified=True,primary_eligible=eligible,
                high_confidence_auc=eligible and spec['study']!='K10262',high_confidence_cmax=eligible and raw_cmax!='',
                has_numeric_curve=False,qc_flags=';'.join(flags),curation_spec_id=identity([spec['study'],protocol])))
        audit.append(dict(study_id=spec['study'],source_sha256=meta['sha256'],cvtdb_document_id=spec['cvtdb_document_id']))
    out=pd.DataFrame(records)
    out.to_csv(collector.root/'data/interim/cebs_endpoints.csv',index=False)
    return out,audit
