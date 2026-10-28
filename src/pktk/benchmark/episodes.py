"""Endpoint-specific episodes and compound-study connected components."""
import json
import numpy as np
import pandas as pd
from pktk.standardization.core import identity

STRATA=['endpoint','species','route','matrix','exposure_definition','exposure_unit','regimen']

def components(x):
    parent={}
    def find(a):
        parent.setdefault(a,a)
        if parent[a]!=a:parent[a]=find(parent[a])
        return parent[a]
    def join(a,b):
        a,b=find(a),find(b)
        if a!=b:parent[max(a,b)]=min(a,b)
    for _,r in x.iterrows():
        join('m:'+r.molecular_group_id,'s:'+r.study_id)
        join('m:'+r.molecular_group_id,'source:'+r.compound_id)
    return {m:identity(['dependency_component',find('m:'+m)]) for m in sorted(x.molecular_group_id.unique())}

def endpoint_table(endpoints,series):
    mapping=components(endpoints);rows=[]
    accepted=series[series.primary_eligible]
    for _,s in accepted.sort_values('series_id').iterrows():
        ids=json.loads(s.source_record_ids);g=endpoints[endpoints.source_record_id.isin(ids)]
        if len(g)!=len(ids):raise ValueError('Missing selected source record')
        for endpoint in ['auc','cmax']:
            if not s['eligible_'+endpoint]:continue
            use=g[g['high_confidence_'+endpoint]].copy()
            use['exposure']=pd.to_numeric(use[endpoint],errors='coerce')
            if not (np.isfinite(use.exposure).all() and use.exposure.gt(0).all()):
                raise ValueError('Nonpositive or nonfinite eligible exposure; resolve source QC before benchmarking')
            if not (np.isfinite(use.dose_mg_kg).all() and use.dose_mg_kg.gt(0).all()):
                raise ValueError('Nonpositive or nonfinite eligible dose; resolve source QC before benchmarking')
            definition_col='auc_definition' if endpoint=='auc' else None
            definitions=use[definition_col] if definition_col else pd.Series('observed_Cmax',index=use.index)
            # Eligibility is endpoint-specific; never compare unlike windows/units.
            if definitions.nunique()!=1 or use[endpoint+'_unit'].nunique()!=1:
                raise ValueError('Incompatible eligible endpoint definitions')
            if use.dose_mg_kg.duplicated().any():raise ValueError('Repeated dose group')
            if len(use)<2:continue
            for _,r in use.sort_values('dose_mg_kg').iterrows():
                d={k:r[k] for k in ['source_record_id','compound_id','molecular_group_id','compound_name','study_id','protocol_id',
                    'source','species','strain','sex','route','matrix','regimen','vehicle','endpoint_origin','qc_flags','scaffold','inchikey']}
                d.update(series_id=s.series_id,endpoint_series_id=identity([s.series_id,endpoint]),endpoint=endpoint,
                    component_id=mapping[r.molecular_group_id],dose=float(r.dose_mg_kg),exposure=float(r.exposure),
                    exposure_definition=definitions.iloc[0],exposure_unit=r[endpoint+'_unit'],
                    source_url=r.source_url,raw_sha256=r.raw_sha256)
                rows.append(d)
    return pd.DataFrame(rows).sort_values(['endpoint_series_id','dose']).reset_index(drop=True)

def generate(table):
    episodes=[];membership=[]
    for sid,g in table.groupby('endpoint_series_id',sort=True):
        g=g.sort_values('dose').reset_index(drop=True);n=len(g)
        plans=[]
        for target in range(1,n):plans.append(('one_extreme','low_to_high',[0],target))
        for target in range(n-1):plans.append(('one_extreme','high_to_low',[n-1],target))
        # For two-dose series, one-extreme already is leave-one-dose-out; avoid duplicates.
        if n>=3:
            for target in range(n):
                support=[j for j in range(n) if j!=target]
                task='high_to_low' if target==0 else 'low_to_high' if target==n-1 else 'interpolation'
                plans.append(('leave_one_dose_out',task,support,target))
        for policy,task,support,target in plans:
            obs=g.iloc[target];train=g.iloc[support]
            ep=identity([sid,policy,task,obs.source_record_id])
            d=obs.to_dict();d.pop('exposure');d.pop('dose')
            ratio=obs.dose/train.dose.max() if task=='low_to_high' else train.dose.min()/obs.dose if task=='high_to_low' else 1.
            d.update(episode_id=ep,policy=policy,task=task,n_support=len(train),target_dose=float(obs.dose),
                     observed=float(obs.exposure),target_record_id=obs.source_record_id,
                     dose_ratio=float(ratio),log_dose_leverage=float(1/len(train)+(np.log(obs.dose)-np.log(train.dose).mean())**2/np.sum((np.log(train.dose)-np.log(train.dose).mean())**2)) if len(train)>=2 else None,support_span=float(train.dose.max()/train.dose.min()),
                     support_ids_json=json.dumps(train.source_record_id.tolist()))
            episodes.append(d)
            for _,r in g.iloc[support+[target]].iterrows():
                membership.append(dict(episode_id=ep,endpoint_series_id=sid,source_record_id=r.source_record_id,
                    role='target' if r.source_record_id==obs.source_record_id else 'support',dose=r.dose))
    out=pd.DataFrame(episodes);members=pd.DataFrame(membership)
    if out.episode_id.duplicated().any() or members.duplicated(['episode_id','source_record_id']).any():
        raise ValueError('Duplicate evaluation episode or membership')
    return out,members

def audit_episode(ep,support,reference):
    if ep.target_record_id in set(support.source_record_id):raise ValueError('Target leakage into support')
    if set(support.endpoint_series_id)!={ep.endpoint_series_id}:raise ValueError('Cross-series support')
    if ep.component_id in set(reference.component_id):raise ValueError('Compound/study leakage into pooled training')
    if ep.molecular_group_id in set(reference.molecular_group_id) or ep.study_id in set(reference.study_id):
        raise ValueError('Identity/study reference leakage')
    if ep.task=='low_to_high' and not ep.target_dose>support.dose.max():raise ValueError('Not low-to-high extrapolation')
    if ep.task=='high_to_low' and not ep.target_dose<support.dose.min():raise ValueError('Not high-to-low extrapolation')
    if ep.task=='interpolation' and not support.dose.min()<ep.target_dose<support.dose.max():raise ValueError('Not interpolation')
    for k in STRATA:
        if len(reference) and not reference[k].eq(ep[k]).all():raise ValueError('Incompatible pooled reference stratum')
