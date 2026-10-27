import json
import numpy as np
import pandas as pd
from pktk.standardization.core import identity

GROUP_FIELDS=['compound_id','species','strain','sex','route','matrix','document_id','regimen','vehicle','fasting','age','age_unit','dose_duration','dose_duration_unit']

def curves_and_candidates(observations):
    curves=[]
    for exp,g in observations.groupby('experiment_id',sort=True):
        r=g.iloc[0]
        flags=[]
        if not g.eligible_for_curve.all(): flags.append('blocked_observation')
        if g.time_h.nunique()<3: flags.append('insufficient_sampling_times')
        if g.time_h.duplicated().any(): flags.append('duplicate_sampling_times')
        if g.concentration_ug_ml.max()<=0: flags.append('nonpositive_exposure')
        base={c:r[c] for c in GROUP_FIELDS}
        base.update(experiment_id=exp,study_id=r.study_id,dose_mg_kg=r.dose_mg_kg,n_observations=len(g),
                    candidate_group_id=identity(r[c] for c in GROUP_FIELDS),n_animals=r.n_animals,
                    qc_flags=';'.join(flags),usable_curve=not flags,high_confidence=False,
                    time_start_h=float(g.time_h.min()),time_end_h=float(g.time_h.max()),
                    auc=None,cmax=None,auc_definition='',auc_unit='ug*h/mL',cmax_unit='ug/mL',
                    endpoint_origin='computed_from_source_curve',source_record_ids=json.dumps(g.source_record_id.tolist()))
        if not flags:
            s=g.sort_values('time_h')
            base['auc']=float(np.trapezoid(s.concentration_ug_ml.to_numpy(float),s.time_h.to_numpy(float)))
            base['cmax']=float(s.concentration_ug_ml.max())
            base['auc_definition']=f'AUC_observed_{base["time_start_h"]:.12g}_{base["time_end_h"]:.12g}_h_linear_trapezoid'
        curves.append(base)
    curves=pd.DataFrame(curves)
    candidates=[]
    for key,g in curves[curves.usable_curve].groupby('candidate_group_id'):
        doses=sorted(float(x) for x in g.dose_mg_kg.unique())
        r=g.iloc[0]
        candidate={c:r[c] for c in GROUP_FIELDS}
        candidate.update(candidate_series_id=key,n_doses=len(doses),dose_min=min(doses),dose_max=max(doses),
                         dose_ratio=max(doses)/min(doses),doses_json=json.dumps(doses),n_curves=len(g),
                         has_auc=True,has_cmax=True,has_numeric_curve=True,
                         comparable_auc_window=g.auc_definition.nunique()==1,
                         metadata_completeness=sum(bool(r[c]) for c in GROUP_FIELDS)/len(GROUP_FIELDS),
                         high_confidence=False,suitability='candidate_only_missing_strain_and_protocol',
                         experiment_ids=json.dumps(g.experiment_id.tolist()),
                         study_ids=json.dumps(sorted(g.study_id.unique())),
                         qc_flags='unverified_protocol;missing_strain'+(';insufficient_dose_levels' if len(doses)<2 else ''))
        candidates.append(candidate)
    return curves,pd.DataFrame(candidates)

def reported_series(endpoints):
    keys=['compound_id','species','strain','sex','route','matrix','study_id','protocol_id','regimen','vehicle','fasting']
    records=[]
    for vals,g in endpoints.groupby(keys):
        doses=sorted(float(x) for x in g.dose_mg_kg.unique())
        row=dict(zip(keys,vals))
        row.update(candidate_series_id=identity(vals),n_doses=len(doses),dose_min=min(doses),dose_max=max(doses),
                   dose_ratio=max(doses)/min(doses),doses_json=json.dumps(doses),has_auc=g.auc.notna().all(),
                   has_cmax=g.cmax.notna().all(),has_numeric_curve=False,metadata_completeness=1.,
                   high_confidence=bool(g.high_confidence_cmax.all() and len(doses)>=2),
                   suitability='Cmax_pilot_only;AUC_terminal_time_unverified',
                   experiment_ids=json.dumps(g.experiment_id.tolist()),qc_flags='auc_terminal_time_unverified;missing_structure')
        records.append(row)
    return pd.DataFrame(records)
