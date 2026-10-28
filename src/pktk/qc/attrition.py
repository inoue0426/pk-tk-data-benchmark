"""Audit design granularity and filtering without promoting exploratory recoveries."""
import json
from pathlib import Path
import pandas as pd
from pktk.benchmark.series import GROUP_FIELDS
from pktk.qc.checks import BLOCKING


def count_series(frame, keys):
    counts=frame.groupby(keys,dropna=False).dose_mg_kg.nunique()
    return int((counts>=2).sum())


def audit(root):
    root=Path(root)
    x=pd.read_csv(root/'data/processed/observations.csv',low_memory=False,keep_default_na=False)
    raw=pd.read_csv(root/'data/interim/cvtdb_original.csv',low_memory=False,keep_default_na=False)
    for c in ['dose_mg_kg','concentration_ug_ml','time_h']:
        x[c]=pd.to_numeric(x[c],errors='coerce')
    x['censor_code']=raw.conc.astype(str).str.upper().to_numpy()
    x['curation_set_tag']=raw.curation_set_tag.to_numpy()
    size=x.groupby('experiment_id').size()
    stages={
        'input_observations':len(x),'source_series_ids':len(size),
        'singleton_series':int((size==1).sum()),
        'series_with_at_least_three_rows':int((size>=3).sum()),
        'ND_NQ_observations':int(x.censor_code.isin(['ND','NQ']).sum()),
        'missing_normalized_concentrations':int(x.censor_code.eq('').sum()),
        'source_study_ids_with_more_than_one_dose':int((x.groupby('study_id').dose_mg_kg.nunique()>1).sum()),
    }
    profiles=[]
    for exp,g in x.groupby('experiment_id'):
        flags=[set(filter(None,s.split(';')))&BLOCKING for s in g.qc_flags]
        usable=g.loc[[not f for f in flags]]
        # ND/NQ removal is an exploratory sensitivity only. Unknown missingness stays blocking.
        censored=g.censor_code.isin(['ND','NQ'])
        safe=all(not f or (f=={'invalid_concentration'} and is_censored) for f,is_censored in zip(flags,censored))
        profile=g.iloc[0].to_dict()
        profile.update(experiment_id=exp,n_input=len(g),n_quantified=len(usable),n_censored=int(censored.sum()),
                       recoverable_censored_profile=bool(safe and usable.time_h.nunique()>=3 and not usable.time_h.duplicated().any()),
                       strict_profile=bool(g.eligible_for_curve.all() and g.time_h.nunique()>=3 and not g.time_h.duplicated().any() and g.concentration_ug_ml.max()>0))
        profiles.append({k:profile[k] for k in GROUP_FIELDS+['experiment_id','study_id','dose_mg_kg','n_input','n_quantified','n_censored','recoverable_censored_profile','strict_profile']})
    p=pd.DataFrame(profiles)
    strict=p[p.strict_profile]
    recovered=p[p.recoverable_censored_profile]
    stages.update(strict_profiles=len(strict),strict_crossdose_series=count_series(strict,GROUP_FIELDS),
                  censor_tolerant_profiles=len(recovered),censor_tolerant_crossdose_series=count_series(recovered,GROUP_FIELDS),
                  censor_tolerant_compounds=recovered.compound_id.nunique())
    essential=['compound_id','species','sex','route','matrix','document_id','regimen']
    stages['strict_profiles_less_restrictive_metadata_groups']=count_series(strict,essential)
    stages['censor_tolerant_less_restrictive_metadata_groups']=count_series(recovered,essential)
    # Count potential destructive-sampling cells, never fabricate longitudinal curves from animal IDs.
    valid=x[x.eligible_for_curve.astype(bool)]
    cells=valid.groupby(['study_id']+essential+['dose_mg_kg'],dropna=False).agg(n_time=('time_h','nunique'),n_rows=('source_record_id','size'),n_subjects=('subject_id','nunique')).reset_index()
    cells=cells[cells.n_time>=3]
    stages['dose_study_cells_with_three_quantified_times']=len(cells)
    stages['potential_crossdose_groups_if_composite_sampling_confirmed']=count_series(cells,essential)
    by_tag=x.assign(singleton=x.experiment_id.map(size).eq(1)).groupby('curation_set_tag').agg(observations=('source_record_id','size'),singleton_observations=('singleton','sum')).reset_index()
    alias=raw.groupby('analyzed_chem_dtxsid',dropna=False).agg(source_ids=('fk_analyzed_chemical_id','nunique'),cas_values=('analyzed_chem_casrn',lambda s:len(set(s)-{''})))
    stages['dtxsid_with_multiple_source_chemical_ids']=int((alias.loc[alias.index!='','source_ids']>1).sum())
    stages['dtxsid_with_multiple_cas_values']=int((alias.loc[alias.index!='','cas_values']>1).sum())
    stages['interpretation']='Sensitivity counts are not validated new series. Never pool across publication IDs; subject-level and composite profiles remain distinct.'
    p.to_csv(root/'data/processed/cvtdb_profile_audit.csv',index=False)
    cells.to_csv(root/'data/processed/cvtdb_composite_candidates.csv',index=False)
    by_tag.to_csv(root/'reports/cvtdb_attrition_by_tag.csv',index=False)
    (root/'reports/cvtdb_attrition.json').write_text(json.dumps(stages,indent=2)+'\n')
    (root/'reports/cvtdb_attrition.md').write_text('# Why 35,669 observations produced ten candidates\n\n'+
        '```json\n'+json.dumps(stages,indent=2)+'\n```\n\n'+
        'The original parser retains every row. Source series IDs often identify individual animals with one sample; requiring three times per ID cannot reconstruct a destructive-sampling population curve. Many multi-time series contain ND/NQ concentrations, causing the original all-or-nothing curve gate to fail. ND/NQ must not be silently converted to zero; numeric-only maxima may miss the true peak and truncated AUC changes its definition. The audit measures recovery potential but does not approve these endpoints.\n\n'+
        'The metadata sensitivity removes vehicle/fasting/age/duration from grouping while retaining publication, compound, species, sex, route, matrix and regimen. Its increase is an upper bound needing protocol review, not a new verified count. Dose-specific source study IDs contain no multi-dose series; publication-level matching is necessary but insufficient. Multiple source chemical IDs can share a DTXSID, and some DTXSIDs carry multiple CAS values; identity conflicts need explicit resolution.\n\n'+
        'Proposed composite-sampling cells preserve source study and animal IDs and are saved separately. No independent publications or individual profiles are merged into verified experiments. Prioritize original NTP protocol/strain/formulation retrieval before approving these cells.\n')
    return stages
