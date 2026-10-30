"""Resumable experiments. Validation selects models; calibration/test never tune them."""
import argparse,json,time,platform,os,sys,traceback,shutil
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import joblib,torch
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
from .data import ROOT,OUT,sha,dump,connected_groups
from .features import load_features,nearest_similarity
from .neural import fit_neural

NEURAL=['mlp','gcn','gin','gine']
RESULTS=ROOT/'results/ml'


def now():return datetime.now(timezone.utc).isoformat()
def event(**kwargs):
    RESULTS.mkdir(exist_ok=True,parents=True)
    with (RESULTS/'events.jsonl').open('a') as f:f.write(json.dumps(dict(time=now(),**kwargs),default=str)+'\n')


def weighted_metrics(d,pred):
    z=d[['family','publication_id','y','value']].copy();z['pred']=pred
    z['ae']=abs(pred-z.y);z['se']=(pred-z.y)**2
    raw=10.**np.asarray(pred);z['raw_ae']=abs(raw-z.value);z['raw_se']=(raw-z.value)**2
    z['are']=abs(raw/z.value-1)
    m=z.groupby('family')[['ae','se','raw_ae','raw_se','are']].mean().mean()
    weight=z.family.map(1/z.family.value_counts()).to_numpy()
    return dict(mae_log10=float(m.ae),rmse_log10=float(np.sqrt(m.se)),r2_log10=float(r2_score(z.y,pred,sample_weight=weight)),
                mae_original=float(m.raw_ae),rmse_original=float(np.sqrt(m.raw_se)),mean_absolute_relative_error=float(m.are),
                r2_original=float(r2_score(z.value,raw,sample_weight=weight)),
                n_records=len(d),n_compounds=d.inchikey.nunique(),n_families=d.family.nunique(),n_publications=d.publication_id.nunique(),
                study_macro_mae_log10=float(z.groupby(['publication_id','family']).ae.mean().groupby('publication_id').mean().mean()))


def grids(name):
    if name in ['mean','median']:return [{}]
    if name in ['ridge','dose_ridge']:return [{'alpha':a} for a in [.1,1.,10.,100.]]
    if name=='rf':return [dict(n_estimators=300,min_samples_leaf=2,max_features=.3),dict(n_estimators=300,min_samples_leaf=5,max_features=1.)]
    if name=='xgb':return [dict(max_depth=2,n_estimators=600,learning_rate=.03,min_child_weight=3),dict(max_depth=4,n_estimators=600,learning_rate=.03,min_child_weight=5)]
    return [dict(hidden=32,lr=.001),dict(hidden=64,lr=.0005)]


def experiment(table,desc,fps,graphs,split,task,kind,seed,name,cfg,train_fraction=1.,ablate=False,member=0,fixed_params=None):
    suffix=(f'__fraction{train_fraction:g}' if train_fraction!=1 else '')+('__ablated' if ablate else '')+(f'__member{member}' if member else '')
    runid=f'{task}__{kind}__s{seed}__{name}'+suffix
    folder=RESULTS/'runs'/runid;folder.mkdir(parents=True,exist_ok=True)
    manifest=dict(task=task,split=kind,seed=seed,model=name,train_fraction=train_fraction,ablate=ablate,member=member,
                  task_sha256=sha(OUT/'tasks.csv'),split_sha256=sha(OUT/'splits.csv'),config_sha256=sha(ROOT/'configs/ml/benchmark.json'),
                  code_sha256={p.name:sha(p) for p in (ROOT/'src/pktk/ml').glob('*.py') if p.name in ['run.py','features.py','neural.py','data.py']})
    if (folder/'metrics.json').exists():
        assert json.loads((folder/'manifest.json').read_text())==manifest,'Stale checkpoint: '+runid
        return json.loads((folder/'metrics.json').read_text())
    dump(folder/'manifest.json',manifest);event(run=runid,status='started')
    start=time.monotonic();random_seed=seed+member*10000
    indices={r:table.index[table.row_id.isin(split[split.partition==r].row_id)].to_numpy() for r in ['train','validation','calibration','test']}
    if train_fraction<1:
        subset=table.iloc[indices['train']];groups=connected_groups(subset,'scaffold_study');unique=sorted(groups.unique());rng=np.random.default_rng(seed);rng.shuffle(unique)
        keep=set(unique[:max(1,int(len(unique)*train_fraction))]);indices['train']=np.array([i for i in subset.index if groups[i] in keep])
    tr=indices['train'];va=indices['validation'];cal=indices['calibration'];test=indices['test']
    assert set(tr).isdisjoint(set(va)|set(cal)|set(test))
    assert all(set(table.iloc[tr].family).isdisjoint(table.iloc[ids].family) for ids in [va,cal,test])
    dose=table.log_dose.to_numpy(float).reshape(-1,1)
    # Learned preprocessing only sees train rows. Checkpoints include fitted objects.
    imputer=SimpleImputer(strategy='median',keep_empty_features=True);imputer.fit(desc[tr]);clean=imputer.transform(desc)
    scaler=StandardScaler().fit(np.column_stack([clean[tr],dose[tr]]));dx=scaler.transform(np.column_stack([clean,dose]))
    dose_scaler=StandardScaler().fit(dose[tr]);dos=dose_scaler.transform(dose)
    x=dx if name=='ridge' else dos if name=='dose_ridge' else np.column_stack([dx,fps]) if name=='mlp' else np.column_stack([fps,dos])
    joblib.dump(dict(imputer=imputer,scaler=scaler,dose_scaler=dose_scaler,train_row_ids=table.iloc[tr].row_id.tolist()),folder/'preprocessing.joblib')
    y=table.y.to_numpy(float);counts=table.iloc[tr].family.value_counts();weights=table.iloc[tr].family.map(1/counts).to_numpy();weights=weights/weights.mean()
    candidates=[];best=float('inf');best_predictions=None;best_info=None;best_params=None
    params_list=[fixed_params] if fixed_params is not None else grids(name)
    for ci,params in enumerate(params_list):
        cstart=time.monotonic();info={};predictions={}
        if name in NEURAL:
            predictions,info=fit_neural(name,graphs,x,table,indices,params,random_seed,cfg,folder/f'candidate_{ci}.pt',ablate)
        else:
            if name in ['mean','median']:
                means=table.iloc[tr].groupby('family').y.mean()
                constant=float(means.mean() if name=='mean' else means.median())
                model=None;predictions={r:np.full(len(ids),constant) for r,ids in indices.items() if r!='train'}
            else:
                if name in ['ridge','dose_ridge']:model=Ridge(**params)
                elif name=='rf':model=RandomForestRegressor(**params,n_jobs=cfg['threads'],random_state=random_seed)
                elif name=='xgb':
                    from xgboost import XGBRegressor
                    model=XGBRegressor(**params,objective='reg:squarederror',tree_method='hist',n_jobs=cfg['threads'],random_state=random_seed,
                                       subsample=.85,colsample_bytree=.8,reg_lambda=10,early_stopping_rounds=25,eval_metric='mae')
                if name=='xgb':
                    vw=table.iloc[va].family.map(1/table.iloc[va].family.value_counts()).to_numpy();vw=vw/vw.mean()
                    model.fit(x[tr],y[tr],sample_weight=weights,eval_set=[(x[va],y[va])],sample_weight_eval_set=[vw],verbose=False)
                    info['best_iteration']=int(model.best_iteration)
                else:model.fit(x[tr],y[tr],sample_weight=weights)
                predictions={r:model.predict(x[ids]) for r,ids in indices.items() if r!='train'}
            joblib.dump(model if model is not None else {'constant':constant},folder/f'candidate_{ci}.joblib')
        val=weighted_metrics(table.iloc[va],predictions['validation'])['mae_log10']
        info.update(params=params,validation_mae_log10=val,seconds=time.monotonic()-cstart,candidate=ci)
        candidates.append(info)
        if val<best:
            best=val;best_predictions=predictions;best_info=info;best_params=params
    dump(folder/'tuning.json',candidates)
    records=[]
    for role in ['validation','calibration','test']:
        idx=indices[role];d=table.iloc[idx].copy();d['partition']=role;d['prediction_log10']=best_predictions[role]
        d['abs_error_log10']=abs(d.prediction_log10-d.y);d['nearest_train_tanimoto']=nearest_similarity(fps,tr,idx)
        d['split_group']=d.row_id.map(split.set_index('row_id').split_group)
        records.append(d[['row_id','inchikey','family','publication_id','split_group','species','source','unit','log_dose','y','value','partition','prediction_log10','abs_error_log10','nearest_train_tanimoto']])
    pred=pd.concat(records);pred.to_csv(folder/'predictions.csv',index=False,float_format='%.12g')
    metrics=weighted_metrics(table.iloc[test],best_predictions['test'])
    metrics.update(run_id=runid,task=task,split=kind,seed=seed,model=name,train_fraction=train_fraction,ablate=ablate,member=member,
                   selected_params=best_params,selected_candidate=best_info['candidate'],validation_mae_log10=best,
                   parameters=best_info.get('parameters'),epochs=best_info.get('epochs'),best_epoch=best_info.get('best_epoch'),
                   seconds=time.monotonic()-start,n_train_records=len(tr),n_train_families=table.iloc[tr].family.nunique(),
                   train_row_ids=table.iloc[tr].row_id.tolist(),test_dependency_groups=pred[pred.partition=='test'].split_group.nunique(),status='completed')
    dump(folder/'metrics.json',metrics);event(run=runid,status='completed',seconds=metrics['seconds'])
    print(runid,'completed',round(metrics['seconds'],1),'seconds',flush=True)
    return metrics


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['baselines','graphs','ablations','ensembles','all'],default='all');parser.add_argument('--seed',type=int);parser.add_argument('--task');parser.add_argument('--models',nargs='+');parser.add_argument('--limit',type=int)
    args=parser.parse_args();cfg=json.loads((ROOT/'configs/ml/benchmark.json').read_text());torch.set_num_threads(cfg['threads'])
    table=pd.read_csv(OUT/'tasks.csv',keep_default_na=False);splits=pd.read_csv(OUT/'splits.csv',keep_default_na=False);audit=json.loads((OUT/'split_audit.json').read_text())
    desc,fps,graphs=load_features(table)
    dump(RESULTS/'hardware.json',dict(platform=platform.platform(),python=sys.version,cpu_count=os.cpu_count(),torch=torch.__version__,mps_available=torch.backends.mps.is_available(),cuda_available=torch.cuda.is_available(),device_used='cpu',threads=cfg['threads']))
    models=args.models or [m for m in cfg['models'] if args.stage=='all' or (m in NEURAL if args.stage=='graphs' else m not in NEURAL)]
    jobs=[]
    for meta in audit:
        if not meta['eligible']:continue
        task,kind,seed=meta['task'],meta['split'],meta['seed']
        if args.seed is not None and seed!=args.seed:continue
        if args.task and task!=args.task:continue
        split=splits[(splits.task==task)&(splits.split==kind)&(splits.seed==seed)]
        # Study and scaffold+study coincide for clearance: exact partition reuse, not a new experiment.
        if kind=='study' and task=='rat_iv_plasma_clearance':continue
        if args.stage in ['ablations','ensembles']:
            # Predeclared focused follow-ups on fixed seed 42, study-disjoint clearance.
            if seed!=42 or kind!='scaffold_study' or task!='rat_iv_plasma_clearance':continue
            for model in ['xgb','gine']:
                base=RESULTS/'runs'/f'{task}__{kind}__s{seed}__{model}'/'metrics.json'
                if not base.exists():raise ValueError('Run primary models before follow-ups')
                params=json.loads(base.read_text())['selected_params']
                if args.stage=='ensembles':
                    for member in [1,2]:jobs.append((split,task,kind,seed,model,1.,False,member,params))
                else:
                    for fraction in [.35,.65]:jobs.append((split,task,kind,seed,model,fraction,False,0,None))
                    if model=='gine':jobs.append((split,task,kind,seed,model,1.,True,0,None))
        else:
            for model in models:jobs.append((split,task,kind,seed,model,1.,False,0,None))
    queue=[dict(task=j[1],split=j[2],seed=j[3],model=j[4],fraction=j[5],ablate=j[6],member=j[7]) for j in jobs]
    dump(RESULTS/f'queue_{args.stage}.json',queue)
    n=0;failures=[]
    for split,task,kind,seed,model,fraction,ablate,member,params in jobs:
        if args.limit is not None and n>=args.limit:break
        n+=1
        try:
            experiment(table,desc,fps,graphs,split,task,kind,seed,model,cfg,fraction,ablate,member,params)
        except Exception as exc:
            failures.append(dict(task=task,split=kind,seed=seed,model=model,error=repr(exc)));event(status='failed',**failures[-1]);traceback.print_exc()
    dump(RESULTS/f'failures_{args.stage}.json',failures)
    if failures:raise SystemExit('Some runs failed; inspect failures and resume completed checkpoints')

if __name__=='__main__':main()
