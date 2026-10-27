import json
import numpy as np
import pandas as pd
from pktk.models.scaling import linear_scaling,power_law
from pktk.standardization.core import identity

def check_split(support,target):
    assert set(support.experiment_id).isdisjoint(set(target.experiment_id)), 'Record leakage'
    assert support.dose_mg_kg.max()<target.dose_mg_kg.min(), 'Not extrapolation'
    for c in ['compound_id','species','strain','sex','route','matrix','study_id','protocol_id']:
        assert support[c].nunique()==target[c].nunique()==1 and support[c].iloc[0]==target[c].iloc[0], 'Incompatible protocol'
    assert not support.dose_mg_kg.duplicated().any(), 'Replicate grouping requires adjudication'

def evaluate(endpoints):
    """Task A episodes only; no global training, tuning or test-informed calibration."""
    splits=[]; predictions=[]; nonlinear=[]
    keys=['compound_id','species','strain','sex','route','matrix','study_id','protocol_id']
    eligible=endpoints[endpoints.high_confidence_cmax & (endpoints.cmax>0)].copy()
    for vals,g in eligible.groupby(keys):
        g=g.sort_values('dose_mg_kg')
        if g.dose_mg_kg.nunique()<2: continue
        group_id=identity(vals)
        for n_support in range(1,len(g)):
            support=g.iloc[:n_support]; target=g.iloc[n_support:]
            check_split(support,target)
            episode=group_id+':k'+str(n_support)
            for role,part in [('support',support),('target',target)]:
                for _,r in part.iterrows():
                    splits.append(dict(task='A',episode_id=episode,series_id=group_id,role=role,
                        experiment_id=r.experiment_id,compound_id=r.compound_id,study_id=r.study_id,
                        dose_mg_kg=r.dose_mg_kg,endpoint='Cmax',regime='extrapolation'))
            for name,model in [('linear_scaling',linear_scaling),('power_law',power_law)]:
                if name=='power_law' and n_support<2: continue
                for _,r in target.iterrows():
                    prediction=model(support.dose_mg_kg.tolist(),support.cmax.tolist(),r.dose_mg_kg)
                    predictions.append(dict(episode_id=episode,series_id=group_id,model=name,n_support=n_support,
                        experiment_id=r.experiment_id,compound_id=r.compound_id,species=r.species,sex=r.sex,
                        endpoint='Cmax',dose_ratio=r.dose_mg_kg/support.dose_mg_kg.max(),
                        observed=r.cmax,predicted=prediction,log_error=float(np.log(prediction/r.cmax)),
                        absolute_relative_error=abs(prediction/r.cmax-1)))
        if len(g)>=3:
            slope=float(np.polyfit(np.log(g.dose_mg_kg.to_numpy(float)),np.log(g.cmax.to_numpy(float)),1)[0])
            nonlinear.append(dict(series_id=group_id,endpoint='Cmax',log_log_slope=slope,
                dose_normalized_fold_range=float((g.cmax/g.dose_mg_kg).max()/(g.cmax/g.dose_mg_kg).min()),
                interpretation='descriptive_only_no_independent_labels_or_uncertainty_test'))
    predictions=pd.DataFrame(predictions)
    metrics=[]
    for dimensions in [['model','n_support'],['model','n_support','dose_ratio'],['model','n_support','compound_id'],['model','n_support','species']]:
        if predictions.empty: break
        for vals,g in predictions.groupby(dimensions):
            row=dict(zip(dimensions,vals))
            row.update(aggregation='/'.join(dimensions),n_predictions=len(g),mae_log=float(g.log_error.abs().mean()),
                rmse_log=float(np.sqrt((g.log_error**2).mean())),mean_absolute_relative_error=float(g.absolute_relative_error.mean()),
                n_independent_compounds=g.compound_id.nunique(),interval_coverage=None)
            metrics.append(row)
    status={'A':'executed_Cmax_pilot' if splits else 'unsupported',
            'B':'not_generated: only one high-confidence compound; molecular identity unresolved',
            'C':'not_generated: no verified molecular structures; no scaffold coverage',
            'D':'descriptive slopes only; no independent nonlinear labels or negative controls',
            'leakage_checks':'support/target IDs disjoint; target doses strictly higher; protocol and identity match',
            'dependence':'episodes reuse groups and targets; never treat predictions as independent replicates',
            'intervals':'not estimated: no independent calibration set',
            'complex_models':'RF/XGBoost/ChemProp and saturable models withheld: one compound, at most two support doses'}
    return pd.DataFrame(splits),predictions,pd.DataFrame(metrics),pd.DataFrame(nonlinear),status
