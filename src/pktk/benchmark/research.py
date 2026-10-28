"""Executed empirical benchmark with explicit episodes, disjoint pooling and clustered summaries."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from pktk.collectors.http import log,sha256
from pktk.benchmark.episodes import endpoint_table,generate,audit_episode,STRATA
from pktk.models.empirical import predict,exponent,NotIdentifiable

MODELS=['linear','power','pooled_power','shrinkage_power','saturating','superlinear']

def dump(path,value):
    path.write_text(json.dumps(value,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x),allow_nan=False)+'\n')

def reference_exponent(reference,minimum):
    slopes=[]
    for _,g in reference.groupby('endpoint_series_id'):
        if len(g)>=2:slopes.append(dict(component_id=g.component_id.iloc[0],slope=exponent(g.dose,g.exposure)))
    if not slopes:return None,0
    means=pd.DataFrame(slopes).groupby('component_id').slope.mean()
    return (float(means.mean()) if len(means)>=minimum else None),len(means)

def predict_episodes(table,episodes,config):
    outputs=[];training=[]
    for _,ep in episodes.iterrows():
        support=table[(table.endpoint_series_id==ep.endpoint_series_id)&table.source_record_id.isin(json.loads(ep.support_ids_json))]
        reference=table[table.component_id!=ep.component_id]
        for key in STRATA:reference=reference[reference[key]==ep[key]]
        audit_episode(ep,support,reference)
        pooled,nref=reference_exponent(reference,config['minimum_pooled_components'])
        training.append(dict(episode_id=ep.episode_id,target_component=ep.component_id,
            reference_components_json=json.dumps(sorted(reference.component_id.unique())),
            reference_record_ids_json=json.dumps(sorted(reference.source_record_id.unique())),
            reference_endpoint_series_json=json.dumps(sorted(reference.endpoint_series_id.unique())),
            n_reference_components=nref,pooled_exponent=pooled,leakage_check='passed'))
        for model in MODELS:
            r=ep.to_dict();r.update(model=model,n_reference_components=nref,predicted=None,status='not_identifiable',reason='',parameters_json='{}')
            try:
                value,params=predict(model,support.dose,support.exposure,ep.target_dose,pooled,
                    config['shrinkage_local_weight'],config['nonlinear_max_jacobian_condition'])
                r.update(predicted=value,status='ok',parameters_json=json.dumps(params,sort_keys=True))
            except NotIdentifiable as e:r['reason']=str(e)
            outputs.append(r)
    out=pd.DataFrame(outputs)
    out['log_error']=np.log(pd.to_numeric(out.predicted)/out.observed)
    out['absolute_log_error']=out.log_error.abs();out['squared_log_error']=out.log_error**2
    out['absolute_relative_error']=(pd.to_numeric(out.predicted)/out.observed-1).abs()
    out['fold_error']=np.exp(out.absolute_log_error)
    out['within_twofold']=out.fold_error.le(2).where(out.status.eq('ok'))
    return out,pd.DataFrame(training)

LOSSES=['absolute_log_error','squared_log_error','absolute_relative_error','within_twofold']

def summarize(g):
    # Average episodes within series, then series within dependency component.
    means=g.groupby(['component_id','endpoint_series_id'])[LOSSES].mean().groupby('component_id').mean()
    return dict(n_predictions=len(g),n_series=g.endpoint_series_id.nunique(),n_compounds=g.molecular_group_id.nunique(),
        n_studies=g.study_id.nunique(),n_components=len(means),mae_ln=float(means.absolute_log_error.mean()),
        rmse_ln=float(np.sqrt(means.squared_log_error.mean())),mean_absolute_relative_error=float(means.absolute_relative_error.mean()),
        within_twofold=float(means.within_twofold.mean()))

def metrics(predictions):
    good=predictions[predictions.status=='ok'].copy();out=[]
    good['ratio_bin']=pd.cut(good.dose_ratio,[0,2,4,10,np.inf],labels=['<=2','>2-4','>4-10','>10']).astype(str)
    keys=['endpoint','policy','task','model','n_support']
    for dimension in [None,'compound_name','study_id','species','route','ratio_bin']:
        groupkeys=keys+([dimension] if dimension else [])
        for values,g in good.groupby(groupkeys,sort=True):
            r=dict(zip(groupkeys,values));r.update(aggregation=dimension or 'overall',group=r.pop(dimension) if dimension else 'all')
            r.update(summarize(g));out.append(r)
    return pd.DataFrame(out)

def paired_comparisons(predictions,config):
    good=predictions[predictions.status=='ok'];out=[];rng=np.random.default_rng(config['seed'])
    keys=['endpoint','policy','task','n_support']
    for vals,g in good.groupby(keys,sort=True):
        base=g[g.model=='linear']
        for model in MODELS[1:]:
            other=g[g.model==model]
            pair=base.merge(other[['episode_id','absolute_log_error']],on='episode_id',suffixes=('_linear','_other'),validate='one_to_one')
            if pair.empty:continue
            pair['difference']=pair.absolute_log_error_other-pair.absolute_log_error_linear
            # Paired, equal-series then equal-component comparison on common episodes only.
            by=pair.groupby(['component_id','endpoint_series_id']).difference.mean().groupby('component_id').mean()
            r=dict(zip(keys,vals));r.update(model=model,n_common_predictions=len(pair),n_components=len(by),
                delta_mae_ln=float(by.mean()),ci_low=None,ci_high=None,interval_status='insufficient_components',
                leave_one_component_out_min=float(min((by.sum()-v)/(len(by)-1) for v in by)) if len(by)>1 else None,
                leave_one_component_out_max=float(max((by.sum()-v)/(len(by)-1) for v in by)) if len(by)>1 else None)
            if len(by)>=config['minimum_bootstrap_components']:
                sample=rng.choice(by.to_numpy(),size=(config['bootstrap_replicates'],len(by)),replace=True).mean(axis=1)
                r.update(ci_low=float(np.quantile(sample,.025)),ci_high=float(np.quantile(sample,.975)),interval_status='exploratory_component_bootstrap')
            out.append(r)
    return pd.DataFrame(out)

def describe_nonlinearity(table):
    rows=[]
    for sid,g in table.groupby('endpoint_series_id',sort=True):
        g=g.sort_values('dose');norm=g.exposure/g.dose
        r={k:g[k].iloc[0] for k in ['endpoint_series_id','series_id','endpoint','compound_name','molecular_group_id','component_id','study_id','species','route','sex','matrix','regimen','exposure_definition','exposure_unit']}
        ratio=float(norm.iloc[-1]/norm.iloc[0])
        r.update(n_doses=len(g),dose_min=float(g.dose.min()),dose_max=float(g.dose.max()),dose_span=float(g.dose.max()/g.dose.min()),
            power_exponent=exponent(g.dose,g.exposure),normalized_high_low_ratio=ratio,
            normalized_max_min_ratio=float(norm.max()/norm.min()),
            adjacent_slopes_json=json.dumps((np.diff(np.log(g.exposure))/np.diff(np.log(g.dose))).tolist()),
            descriptive_trend='superproportional' if ratio>1.25 else 'subproportional' if ratio<.8 else 'near_proportional',
            inferential_status='descriptive_only_no_individual_covariance',target_use='all_doses_descriptive_not_predictor')
        rows.append(r)
    return pd.DataFrame(rows)

def feasibility(table,episodes,predictions):
    rows=[]
    for vals,g in episodes.groupby(['endpoint','policy','task','n_support'],sort=True):
        rows.append(dict(zip(['endpoint','policy','task','n_support'],vals))|dict(episodes=len(g),series=g.endpoint_series_id.nunique(),
            compounds=g.molecular_group_id.nunique(),studies=g.study_id.nunique(),independence_components=g.component_id.nunique()))
    chemistry=table.drop_duplicates('molecular_group_id');resolved=chemistry[chemistry.inchikey.ne('')]
    scaffold_counts=resolved.groupby('scaffold').molecular_group_id.nunique().to_dict()
    return dict(task_coverage=rows,total_compounds=table.molecular_group_id.nunique(),total_studies=table.study_id.nunique(),
        independence_components=table.component_id.nunique(),endpoint_series=table.endpoint_series_id.nunique(),
        resolved_compounds=len(resolved),scaffold_group_counts=scaffold_counts,
        acyclic_compounds=int(resolved.scaffold.eq('').sum()),
        pooled_generalization='Leave-compound-study-component-out slopes; target compound support remains permitted. Not zero-shot molecular prediction.',
        scaffold_generalization='Descriptive scaffold counts only; unknown structures and small endpoint-specific strata preclude robust scaffold-disjoint tuning.',
        uncertainty='Paired component bootstrap describes uncertainty of fixed cross-validated errors. Small, purposively selected source cohort; no population or individual prediction-interval guarantee.',
        nonlinear_models=predictions[predictions.model.isin(['saturating','superlinear'])].groupby(['model','status']).size().reset_index(name='attempts').to_dict('records'))

def run(root=Path('.')):
    root=Path(root);out=root/'results/research';out.mkdir(parents=True,exist_ok=True)
    cfg=json.loads((root/'configs/analysis/benchmark.json').read_text())
    log(root/'reports/progress.jsonl',phase='empirical_research',status='started')
    ep=pd.read_csv(root/'data/processed/high_confidence_endpoints.csv',keep_default_na=False)
    series=pd.read_csv(root/'data/processed/crossdose_series_expanded.csv',keep_default_na=False)
    table=endpoint_table(ep,series);episodes,membership=generate(table)
    predictions,training=predict_episodes(table,episodes,cfg)
    summary=metrics(predictions);paired=paired_comparisons(predictions,cfg);nonlinear=describe_nonlinearity(table)
    sensitivity=[]
    for label,mask in [('all',predictions.status.eq('ok')),('resolved_structures',predictions.inchikey.ne('')),
                       ('exclude_noncommercial',~predictions.qc_flags.str.contains('license_noncommercial'))]:
        v=metrics(predictions[mask]);v=v[v.aggregation.eq('overall')].copy();v['sensitivity']=label;sensitivity.append(v)
    tables={'sensitivity_metrics':pd.concat(sensitivity,ignore_index=True),'endpoint_observations':table,'episodes':episodes,'episode_membership':membership,'predictions':predictions,
        'reference_training':training,'metrics':summary,'paired_comparisons':paired,'nonlinearity':nonlinear}
    for name,data in tables.items():data.to_csv(out/(name+'.csv'),index=False,float_format='%.12g')
    f=feasibility(table,episodes,predictions);dump(out/'feasibility.json',f)
    manifest=dict(config=cfg,inputs={str(p.relative_to(root)):sha256(p) for p in [root/'data/processed/high_confidence_endpoints.csv',root/'data/processed/crossdose_series_expanded.csv']},
        outputs={p.name:sha256(p) for p in sorted(out.glob('*.csv'))},episodes=len(episodes),successful_predictions=int(predictions.status.eq('ok').sum()),
        attempts=len(predictions),leakage_checks_passed=len(training))
    dump(root/'reports/research_manifest.json',manifest)
    from pktk.benchmark.visualize import figures
    figures(root,table,episodes,predictions,summary,nonlinear)
    from pktk.benchmark.scientific_report import write_reports
    write_reports(root,table,episodes,predictions,summary,paired,nonlinear,f)
    log(root/'reports/progress.jsonl',phase='empirical_research',status='succeeded',episodes=len(episodes),components=f['independence_components'],successful_predictions=manifest['successful_predictions'])
    print(json.dumps({k:manifest[k] for k in ['episodes','successful_predictions','attempts','leakage_checks_passed']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path('.'));a=p.parse_args();run(a.root)
