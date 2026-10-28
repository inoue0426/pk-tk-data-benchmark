"""Pinned, explicitly reviewed CEBS rows. No automatic promotion by chemical match."""
import json
import re
import pandas as pd
from pktk.standardization.core import identity

AUC_FACTORS={'mg*min/L':1/60,'ug/mL*min':1/60,'ug*min/mL':1/60,
             'ug/mL*hr':1,'ug mL^-1 min':1/60,'ng*min*mL^-1':1/60000,'ng/mL*hr':.001}
CMAX_FACTORS={'mg/L':1,'ug/mL':1,'ng*mL^-1':.001,'ng/mL':.001}

def numeric(value):
    if value=='':return None
    # Only conventional thousands separators are accepted, not censored measurements.
    text=str(value)
    if ',' in text and not re.fullmatch(r'\d{1,3}(,\d{3})+(\.\d+)?',text):
        raise ValueError('Ambiguous numeric cell: '+text)
    return float(text.replace(',',''))

def collect(collector):
    specs=json.loads((collector.root/'configs/curation/cebs_expansion.json').read_text())
    specs+=json.loads((collector.root/'configs/curation/cebs_targeted.json').read_text())
    records=[];ledger=[];manifests=[]
    checks={'group_name':'GROUP_NAME','species':'SPECIES_COMMON_NAME','sex':'SEX','route':'ROUTE',
            'dose':'DOSE','matrix':'TISSUE','analyte':'ANALYTE'}
    for spec in specs:
        path,meta=collector.fetch(spec['file'],spec['url'],spec['version'],spec['sha256'])
        x=pd.read_excel(path,keep_default_na=False)
        selected={v['row'] for v in spec['rows']}
        for i,r in x.iterrows():
            ledger.append(dict(file=spec['file'],source_row=i+2,group_name=r.get('GROUP_NAME',''),
                selected=i+2 in selected,review_note=spec['review_note'],raw_sha256=meta['sha256']))
        for expected in spec['rows']:
            r=x.iloc[expected['row']-2]
            for key,column in checks.items():
                if r[column]!=expected[key]:raise ValueError(f'CEBS reviewed cell mismatch: {spec["study"]} {column}')
            if str(r.NTP_CHEMTRACK_NUMBER).strip()!=spec['study']:
                raise ValueError('Experiment identifier mismatch')
            if r.VOLUME_PER_ADMIN_UNIT!='mL/kg' or r.DOSE_FREQUENCY!=1 or r.DOSE_FREQUENCY_UNIT not in ['per study','study']:
                raise ValueError('Unreviewed administration volume or dosing regimen')
            if r.DOSE_UNIT!='mg/kg' or str(r.Is_Parent_Compound).strip()!='T':raise ValueError('Unsupported dose/analyte')
            if spec['model_evidence'] and spec['model_evidence'] not in r['Best Fit Model']:
                raise ValueError('Reviewed endpoint method/window missing')
            protocol_evidence={k:r.get(k,'') for k in ['PROTOCOL_DESC','PROTOCOL_DETAILS','Best Fit Model','EXCEPTION_COMMENTS']}
            if not r.PROTOCOL_DESC:raise ValueError('Missing experimental protocol')
            dose=float(r.DOSE)
            flags=['source_exception_preserved','aggregate_endpoint','fasting_not_reported','sampling_schedule_dose_dependent']
            if spec['dose_basis']=='actual':
                m=re.search(r'Actual administered dose is ([0-9.]+) mg/kg',r.EXCEPTION_COMMENTS)
                if not m:raise ValueError('Actual dose not verified')
                dose=float(m[1]);flags.append('actual_not_nominal_dose')
            raw_auc=numeric(r[spec['auc_column']]);raw_cmax=numeric(r.get('Cmax_obs',''))
            auc=raw_auc*AUC_FACTORS[r[spec['auc_column']+'[DEPOSITOR_UNIT]']] if raw_auc is not None else None
            cmax=raw_cmax*CMAX_FACTORS[r['Cmax_obs[DEPOSITOR_UNIT]']] if raw_cmax is not None else None
            if not spec['eligible_auc']:flags.append('auc_definition_or_method_not_eligible')
            elif spec['auc_definition']=='AUC0-inf':flags.append('auc_tail_fraction_unavailable')
            if r.get('Cmax_pred','')!='':flags.append('predicted_cmax_not_observed')
            species=r.SPECIES_COMMON_NAME.lower();route={'Gavage':'oral','Oral gavage':'oral','IV':'intravenous','Intravenous':'intravenous'}[r.ROUTE]
            matrix={'Plasma':'plasma','Whole Blood':'blood'}[r.TISSUE]
            origin='reported_individual_parameter_mean' if spec['study']=='S0571' and species=='rat' else 'reported_composite_profile'
            # Formulation is part of the protocol key, not a post-hoc covariate to ignore.
            protocol='CEBS:'+spec['study']+':'+identity([species,route,matrix,r.VEHICLE_DESC,r.VOLUME_PER_ADMIN,spec['auc_definition']])
            rid=f'CEBS:{spec["study"]}:Sheet1:row{expected["row"]}'
            records.append(dict(source='cebs',source_record_id=rid,experiment_id=rid,compound_id=spec['compound_id'],
                compound_name=str(r.CHEMICAL).strip(),casrn=str(r.CAS_NO).strip(),pubchem_cid='',canonical_smiles='',inchikey='',species=species,
                strain=r.STRAIN,sex=r.SEX.lower(),route=route,matrix=matrix,dose_mg_kg=dose,dose_unit='mg/kg',
                dose_basis=spec['dose_basis'],dose_volume_ml_kg=numeric(r.VOLUME_PER_ADMIN),auc=auc,cmax=cmax,
                auc_sd=None,cmax_sd=None,auc_definition=spec['auc_definition'],auc_unit='ug*h/mL',cmax_unit='ug/mL',
                n_animals=None,n_animals_per_timepoint=None,regimen='once',vehicle=r.VEHICLE_DESC,fasting='not reported',
                infusion_duration_min=None,protocol_id=protocol,study_id='CEBS:'+spec['study'],document_id=spec['cvtdb_document_id'],
                publication_doi=spec['doi'],pmid='',publication_url=spec['landing_url'],source_url=meta['url'],
                retrieved_at=meta['retrieved_at'],dataset_version=meta['version'],raw_sha256=meta['sha256'],
                original_json=json.dumps(r.to_dict(),default=str),protocol_evidence_json=json.dumps(protocol_evidence),
                endpoint_origin=origin,protocol_verified=True,primary_eligible=True,
                high_confidence_auc=spec['eligible_auc'] and auc is not None,high_confidence_cmax=cmax is not None,
                has_numeric_curve=False,qc_flags=';'.join(flags),curation_spec_id=identity([spec['study'],protocol]),
                curation_round=spec.get('curation_round','composite_expansion')))
        manifests.append(dict(file=spec['file'],study_id=spec['study'],selected_rows=len(selected),total_rows=len(x),**meta))
    pd.DataFrame(ledger).to_csv(collector.root/'reports/cebs_row_selection.csv',index=False)
    (collector.root/'reports/cebs_expansion_manifest.json').write_text(json.dumps(manifests,indent=2)+'\n')
    return pd.DataFrame(records)
