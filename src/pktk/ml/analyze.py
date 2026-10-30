"""Read-only empirical analysis of locked experiments and independent calibration."""
import json,math,os,shutil
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from .data import ROOT,OUT,dump,sha,connected_groups

RESULTS=ROOT/'results/ml'


def quantile_radius(scores,alpha):
    """Finite-sample split-conformal rank; infinity if calibration cannot support it."""
    a=np.asarray(scores,float)
    rank=math.ceil((len(a)+1)*(1-alpha))
    return float(np.sort(a)[rank-1]) if len(a) and rank<=len(a) else None


def intervals(pred,alpha=.1):
    cal=pred[pred.partition=='calibration'].copy();test=pred[pred.partition=='test'].copy()
    cal['error']=abs(cal.prediction_log10-cal.y);test['error']=abs(test.prediction_log10-test.y)
    scores=cal.groupby('split_group').error.max()
    radius=quantile_radius(scores,alpha)
    result=dict(nominal=1-alpha,n_calibration_groups=len(scores),n_test_groups=test.split_group.nunique(),radius_log10=radius)
    if radius is None:return result|dict(status='insufficient_calibration_groups',group_coverage=None,family_coverage=None,width_log10=None)
    group=test.groupby('split_group').error.max();family=test.groupby('family').error.max()
    return result|dict(status='empirical_shift_coverage_no_exchangeability_guarantee',group_coverage=float((group<=radius).mean()),family_coverage=float((family<=radius).mean()),width_log10=2*radius)


def primary_metrics():
    rows=[]
    for path in sorted((RESULTS/'runs').glob('*/metrics.json')):
        r=json.loads(path.read_text());r['path']=str(path.parent.relative_to(ROOT));rows.append(r)
    all_metrics=pd.DataFrame(rows)
    main=all_metrics[(all_metrics.train_fraction==1)&(~all_metrics.ablate)&(all_metrics.member==0)].copy()
    main.to_csv(RESULTS/'metrics.csv',index=False)
    columns=['mae_log10','rmse_log10','r2_log10','mae_original','rmse_original','r2_original','mean_absolute_relative_error','study_macro_mae_log10','seconds']
    summary=main.groupby(['task','split','model'])[columns].agg(['mean','std','count'])
    summary.columns=['_'.join(c) for c in summary.columns];summary=summary.reset_index();summary.to_csv(RESULTS/'model_summary.csv',index=False)
    return all_metrics,main,summary


def paired(main,table):
    rows=[];rng=np.random.default_rng(3030)
    for (task,split,seed),g in main.groupby(['task','split','seed']):
        base=g[g.model=='xgb']
        if len(base)!=1:continue
        bp=pd.read_csv(ROOT/base.iloc[0].path/'predictions.csv').query("partition == 'test'")
        for model in ['ridge','rf','mlp','gcn','gin','gine']:
            other=g[g.model==model]
            if len(other)!=1:continue
            op=pd.read_csv(ROOT/other.iloc[0].path/'predictions.csv').query("partition == 'test'")
            p=bp.merge(op[['row_id','prediction_log10']],on='row_id',validate='one_to_one',suffixes=('_base','_other'))
            p['delta']=abs(p.prediction_log10_other-p.y)-abs(p.prediction_log10_base-p.y)
            p['evaluation_group']=p.row_id.map(table.set_index('row_id').dependency_group)
            by=p.groupby(['evaluation_group','family']).delta.mean().groupby('evaluation_group').mean()
            r=dict(task=task,split=split,seed=int(seed),model=model,n_test_families=p.family.nunique(),n_dependency_groups=len(by),
                   delta_family_mae_log10=float(p.groupby('family').delta.mean().mean()),delta_component_mae_log10=float(by.mean()),ci_low=None,ci_high=None,
                   interval_status='withheld_shared_studies_or_small_n')
            if 'study' in split and len(by)>=8:
                sample=rng.choice(by.to_numpy(),size=(2000,len(by)),replace=True).mean(axis=1)
                r.update(ci_low=float(np.quantile(sample,.025)),ci_high=float(np.quantile(sample,.975)),interval_status='exploratory_component_bootstrap_fixed_predictions')
            rows.append(r)
    out=pd.DataFrame(rows);out.to_csv(RESULTS/'paired_comparisons.csv',index=False);return out


def reliability(main):
    intervals_rows=[];risk=[];associations=[];ensemble_rows=[]
    for _,r in main.iterrows():
        p=pd.read_csv(ROOT/r.path/'predictions.csv')
        if 'study' in r.split:
            for alpha in [.5,.2,.1]:intervals_rows.append(dict(task=r.task,split=r.split,seed=int(r.seed),model=r.model,**intervals(p,alpha)))
        test=p[p.partition=='test'].copy()
        f=test.groupby('family')[['abs_error_log10','nearest_train_tanimoto']].mean()
        rho=float(spearmanr(f.abs_error_log10,f.nearest_train_tanimoto).statistic) if f.nearest_train_tanimoto.nunique()>1 else None
        associations.append(dict(task=r.task,split=r.split,seed=int(r.seed),model=r.model,n_families=len(f),similarity_error_spearman=rho))
        for uncertainty,label in [(1-f.nearest_train_tanimoto,'one_minus_nearest_similarity')]:
            order=np.argsort(uncertainty.to_numpy(),kind='stable')
            for coverage in [.2,.4,.6,.8,1.]:
                keep=order[:max(1,math.ceil(len(f)*coverage))]
                risk.append(dict(task=r.task,split=r.split,seed=int(r.seed),model=r.model,uncertainty=label,coverage=len(keep)/len(f),mae_log10=float(f.iloc[keep].abs_error_log10.mean())))
        member_paths=[ROOT/(str(r.path)+f'__member{i}')/'predictions.csv' for i in [1,2]]
        if r.model not in ['xgb','gine'] or not all(path.exists() for path in member_paths):continue
        joined=p.set_index('row_id').copy();preds=[joined.prediction_log10.to_numpy()]
        for path in member_paths:
            member=pd.read_csv(path).set_index('row_id').loc[joined.index]
            assert np.allclose(member.y,joined.y) and (member.partition==joined.partition).all()
            preds.append(member.prediction_log10.to_numpy())
        joined['prediction_log10']=np.mean(preds,axis=0);joined['ensemble_sd']=np.std(preds,axis=0,ddof=1)
        joined['abs_error_log10']=abs(joined.prediction_log10-joined.y);joined=joined.reset_index()
        joined.to_csv(RESULTS/f'ensemble_{r.model}.csv',index=False)
        for alpha in [.5,.2,.1]:intervals_rows.append(dict(task=r.task,split=r.split,seed=int(r.seed),model=r.model+'_ensemble',**intervals(joined,alpha)))
        f=joined[joined.partition=='test'].groupby('family')[['abs_error_log10','ensemble_sd']].mean()
        corr=float(spearmanr(f.abs_error_log10,f.ensemble_sd).statistic) if f.ensemble_sd.nunique()>1 else None
        ensemble_rows.append(dict(task=r.task,split=r.split,seed=int(r.seed),model=r.model,n_members=3,n_test_families=len(f),mae_log10=float(f.abs_error_log10.mean()),uncertainty_error_spearman=corr))
        order=np.argsort(f.ensemble_sd.to_numpy(),kind='stable')
        for coverage in [.2,.4,.6,.8,1.]:
            keep=order[:max(1,math.ceil(len(f)*coverage))]
            risk.append(dict(task=r.task,split=r.split,seed=int(r.seed),model=r.model+'_ensemble',uncertainty='ensemble_sd',coverage=len(keep)/len(f),mae_log10=float(f.iloc[keep].abs_error_log10.mean())))
    for name,rows in [('intervals',intervals_rows),('risk_coverage',risk),('similarity_associations',associations),('ensemble_summary',ensemble_rows)]:pd.DataFrame(rows).to_csv(RESULTS/(name+'.csv'),index=False)
    return pd.DataFrame(intervals_rows),pd.DataFrame(risk),pd.DataFrame(ensemble_rows)


def figures(main,summary,all_metrics,interval,risk):
    os.environ.setdefault('MPLCONFIGDIR','/tmp/pktk-matplotlib');import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.hashsalt':'pktk-ml-v1','axes.spines.top':False,'axes.spines.right':False})
    folder=RESULTS/'figures';folder.mkdir(exist_ok=True)
    def save(fig,name):
        fig.tight_layout();fig.savefig(folder/(name+'.pdf'),bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None});fig.savefig(folder/(name+'.png'),dpi=300,bbox_inches='tight');plt.close(fig)
        presentation=ROOT/'reports/ml/figures';presentation.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(folder/(name+'.pdf'),presentation/(name+'.pdf'))
    audit=json.loads((ROOT/'reports/ml/audit.json').read_text())
    fig,ax=plt.subplots(figsize=(7,3.2));tasks=audit['tasks'];x=np.arange(len(tasks));ax.bar(x-.17,[r['compounds'] for r in tasks],.34,label='Structure-resolved compounds');ax.bar(x+.17,[r['records'] for r in tasks],.34,label='Collapsed endpoint rows');ax.set_xticks(x,['IV plasma clearance','Oral plasma Cmax']);ax.set_ylabel('Count');ax.set_title('ML-eligible conditional rat PK tasks (not cross-dose validation)');ax.legend();save(fig,'01_ml_coverage')
    for task in sorted(main.task.unique()):
        label='Clearance' if 'clearance' in task else 'Cmax'
        g=summary[summary.task==task];fig,axes=plt.subplots(1,2,figsize=(12,4.5))
        for ax,kind in zip(axes,['compound','scaffold_study']):
            v=g[g.split==kind].set_index('model').sort_values('mae_log10_mean')
            ax.barh(v.index,v.mae_log10_mean,xerr=v.mae_log10_std.fillna(0),color=['#c76d2e' if m in ['gcn','gin','gine'] else '#286184' for m in v.index]);ax.set_xlabel('Family-macro MAE (log10 endpoint)');ax.set_title(f'{label}: {kind}; {int(v.mae_log10_count.iloc[0])} split seed(s)')
        fig.suptitle('Error bars: descriptive SD across overlapping splits, not confidence intervals');save(fig,'02_models_'+task)
        fig,ax=plt.subplots(figsize=(7,4));v=g[g.model.isin(['xgb','rf','gine'])]
        piv=v.pivot(index='split',columns='model',values='mae_log10_mean');piv.plot.bar(ax=ax,color=[{'xgb':'#1f77b4','rf':'#2ca02c','gine':'#ff7f0e'}[m] for m in piv.columns]);ax.legend(loc='upper left',bbox_to_anchor=(1,1));ax.set_ylabel('Family-macro MAE (log10)');ax.set_xlabel('Partition design');ax.set_title(label+': compound/scaffold diagnostics versus study exclusion');ax.tick_params(axis='x',rotation=20);save(fig,'03_generalization_'+task)
        # Fixed seed chosen prospectively; do not choose a favorable test split.
        seed=int(main[(main.task==task)&(main.split=='scaffold_study')].seed.min())
        fig,axes=plt.subplots(1,2,figsize=(8,4))
        for ax,model in zip(axes,['xgb','gine']):
            r=main[(main.task==task)&(main.split=='scaffold_study')&(main.seed==seed)&(main.model==model)].iloc[0]
            p=pd.read_csv(ROOT/r.path/'predictions.csv').query("partition == 'test'");ax.scatter(p.y,p.prediction_log10,s=20,alpha=.7)
            lo=min(p.y.min(),p.prediction_log10.min());hi=max(p.y.max(),p.prediction_log10.max());ax.plot([lo,hi],[lo,hi],'k--',lw=1)
            unit='mL/min/kg' if 'clearance' in task else 'µg/mL';ax.set_xlabel(f'Observed log10({unit})');ax.set_ylabel(f'Predicted log10({unit})');ax.set_title(f'{model}; seed {seed}; n={p.family.nunique()} families')
        fig.suptitle(label+': scaffold + study held out; each point is an endpoint row');save(fig,'04_predictions_'+task)
        fig,axes=plt.subplots(1,2,figsize=(8,3.5))
        for ax,model in zip(axes,['xgb','gine']):
            selected=main[(main.task==task)&(main.split=='scaffold_study')&(main.model==model)]
            for _,r in selected.iterrows():
                p=pd.read_csv(ROOT/r.path/'predictions.csv').query("partition == 'test'").groupby('family')[['nearest_train_tanimoto','abs_error_log10']].mean();ax.scatter(p.nearest_train_tanimoto,p.abs_error_log10,s=18,alpha=.6,label=f'seed {r.seed}')
            ax.set_xlabel('Maximum training Morgan Tanimoto');ax.set_ylabel('Absolute log10 error');ax.set_title(model);ax.legend(fontsize=7)
        fig.suptitle(label+': repeated split observations are correlated');save(fig,'05_similarity_'+task)
    # Full error distributions and auditable outliers remain local; no post-hoc removal.
    fig,axes=plt.subplots(1,2,figsize=(9,3.6));outliers=[]
    for ax,task in zip(axes,sorted(main.task.unique())):
        for model in ['xgb','gine']:
            errors=[]
            for _,r in main[(main.task==task)&(main.split=='scaffold_study')&(main.model==model)].iterrows():
                p=pd.read_csv(ROOT/r.path/'predictions.csv').query("partition == 'test'").copy()
                errors.extend(p.groupby('family').abs_error_log10.mean().tolist())
                worst=p.nlargest(5,'abs_error_log10').copy();worst['task']=task;worst['model']=model;worst['seed']=r.seed;outliers.append(worst)
            ax.hist(errors,bins=12,alpha=.45,label=model)
        ax.set_xlabel('Absolute log10 error');ax.set_ylabel('Family / split observations');ax.set_title('Clearance' if 'clearance' in task else 'Cmax');ax.legend()
    fig.suptitle('Strict-split error distributions; repeated seeds are correlated');save(fig,'08_error_distributions')
    pd.concat(outliers).to_csv(RESULTS/'outliers.csv',index=False)
    follow=all_metrics[(all_metrics.task=='rat_iv_plasma_clearance')&(all_metrics.split=='scaffold_study')&(all_metrics.seed==42)&all_metrics.model.isin(['xgb','gine'])&(~all_metrics.ablate)&(all_metrics.member==0)]
    fig,ax=plt.subplots(figsize=(6,3.5))
    for model,g in follow.groupby('model'):
        g=g.sort_values('n_train_families');ax.plot(g.n_train_families,g.mae_log10,'o-',label=model,color={'xgb':'#1f77b4','gine':'#ff7f0e'}[model])
    ax.set_xlabel('Training chemical families');ax.set_ylabel('Test family-macro MAE (log10)');ax.set_title('Clearance learning curve; one fixed scaffold/study split');ax.legend();save(fig,'06_learning_curve')
    fig,axes=plt.subplots(1,2,figsize=(9,3.6))
    for model in ['xgb_ensemble','gine_ensemble']:
        v=interval[interval.model==model].sort_values('nominal');axes[0].plot(v.nominal,v.group_coverage,'o-',label=model)
        q=risk[(risk.model==model)&(risk.uncertainty=='ensemble_sd')];axes[1].plot(q.coverage,q.mae_log10,'o-',label=model)
    axes[0].plot([.5,1],[.5,1],'k--',lw=1);axes[0].set_xlabel('Nominal group coverage');axes[0].set_ylabel('Observed held-out group coverage');axes[0].set_ylim(0,1.05)
    axes[1].set_xlabel('Fraction of test families retained');axes[1].set_ylabel('MAE (log10)');axes[0].legend(fontsize=7);axes[1].legend(fontsize=7)
    fig.suptitle('Clearance: small-sample empirical calibration / selective risk; no shift guarantee');save(fig,'07_uncertainty')


def md(frame,columns=None,digits=3):
    if frame.empty:return '*No supported results.*'
    d=frame[columns] if columns else frame
    def fmt(v):
        if isinstance(v,(float,np.floating)):return '' if not np.isfinite(v) else f'{v:.{digits}f}'
        return str(v)
    return '| '+' | '.join(d.columns)+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(fmt(v) for v in r)+' |' for r in d.itertuples(index=False,name=None))


def reports(all_metrics,main,summary,pairs,interval,ensembles):
    audit=json.loads((ROOT/'reports/ml/audit.json').read_text());splits=json.loads((OUT/'split_audit.json').read_text())
    table=pd.read_csv(OUT/'tasks.csv',keep_default_na=False)
    events=[json.loads(line) for line in (RESULTS/'events.jsonl').read_text().splitlines()]
    scope='These are molecular-plus-known-dose predictions of reported rat plasma PK endpoints, not dose extrapolation, not a validation of the full Dataset A, and not Dataset B transfer. All trained records come from ChEMBL. Original A/B tables and original empirical cross-dose results are unchanged.'
    common='Scores use log10 endpoints (unlike the earlier Dataset B natural-log errors). Equal chemical-family weighting prevents repeated records from dominating. Standard deviations across split seeds are descriptive; test sets overlap and are not independent replications. Values in original units and family/study counts are in results/ml/metrics.csv. No target outliers were removed.'
    def score(task,split,model,column='mae_log10_mean'):
        return float(summary[(summary.task==task)&(summary.split==split)&(summary.model==model)].iloc[0][column])
    cl='rat_iv_plasma_clearance';cm='rat_oral_plasma_cmax'
    findings=(f"**Measured findings.** Strict scaffold+study mean log10 MAE: clearance GINE {score(cl,'scaffold_study','gine'):.3f}, XGBoost {score(cl,'scaffold_study','xgb'):.3f}, descriptor Ridge {score(cl,'scaffold_study','ridge'):.3f}; "
              f"Cmax GINE {score(cm,'scaffold_study','gine'):.3f}, XGBoost {score(cm,'scaffold_study','xgb'):.3f}, descriptor Ridge {score(cm,'scaffold_study','ridge'):.3f}. "
              "These averages describe overlapping splits, not independent replications. Model rankings are exploratory; paired component comparisons below determine the limits of inference.\n\n"
              f"Random-compound Cmax: Random Forest {score(cm,'compound','rf'):.3f}, GINE {score(cm,'compound','gine'):.3f}; "
              f"scaffold-only Cmax: XGBoost {score(cm,'scaffold','xgb'):.3f}, GINE {score(cm,'scaffold','gine'):.3f}. "
              "These diagnostic splits share publications across partitions and do not establish unseen-study generalization. "
              "Negative R² compares squared error with the held-out test mean, not with the separately trained mean predictor. "
              "Ensemble uncertainty and empirical interval coverage are reported separately; small group counts and distribution shift prevent coverage guarantees.")
    strict_pairs=pairs[(pairs.model=='gine')&(pairs.split=='scaffold_study')]
    excluded_zero=int(((strict_pairs.ci_high<0)|(strict_pairs.ci_low>0)).sum())
    findings+=(f"\n\n**Interpretation.** The training-mean control has clearance MAE {score(cl,'scaffold_study','mean'):.3f}, below all molecular models on average. "
               "All models have negative mean strict-split log-space R² for both endpoints. Molecular inputs therefore do not establish useful unseen-study clearance prediction here. "
               f"Of {len(strict_pairs)} prespecified GINE–XGBoost strict comparisons, {excluded_zero} exploratory fixed-prediction component-bootstrap intervals exclude zero in favor of GINE; "
               "the remaining intervals include zero. These unadjusted, overlapping-split results do not establish general graph superiority, and outperforming XGBoost alone does not establish improvement over simpler controls. "
               "Cmax becomes substantially harder when both studies and scaffolds are excluded; partition populations differ, so this is not a causal estimate of a scaffold effect.")
    for _,r in interval[interval.model.str.endswith('_ensemble')&np.isclose(interval.nominal,.9)].iterrows():
        findings+=(f"\n\n{r.model}: 90% nominal intervals use {r.n_calibration_groups} calibration groups and {r.n_test_groups} test groups; "
                   f"observed group coverage {r.group_coverage:.1%}, radius {r.radius_log10:.3f} log10 (×/÷{10**r.radius_log10:.0f}). "
                   "Such wide intervals limit usefulness despite high observed coverage; distribution-shift exchangeability is unverified.")
    source='Source: [ChEMBL](https://www.ebi.ac.uk/chembl/) 37, CC BY-SA 3.0; source publication IDs and hashes remain in the local task tables. Source-derived aggregate tables/figures in these reports are supplied under [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/) with attribution to ChEMBL and the linked source publications; they are transformed analyses, not original database releases; numerical endpoint tables, predictions and model checkpoints are not committed. Code remains under its existing license.'
    columns=['task','split','model','mae_log10_mean','mae_log10_std','rmse_log10_mean','r2_log10_mean','mae_log10_count']
    primary=summary[summary.split=='scaffold_study']
    write=lambda name,s:(ROOT/'reports'/name).write_text(s+'\n')
    write('molecular_ml_benchmark.md','# Molecular PK benchmark\n\n'+scope+'\n\n## Executed task coverage\n\n'+md(pd.DataFrame(audit['tasks']))+'\n\n'+common+'\n\n'+findings+'\n\n## Models and selection\n\nTraining-log mean/median, dose-only Ridge, 17 RDKit descriptors + Ridge, 2048-bit radius-2 chirality-aware Morgan + Random Forest / XGBoost, descriptor+fingerprint MLP, GCN, GIN and GINE were evaluated. Dose is an input for every fitted predictive model; constants are deliberately unconditioned controls. XGBoost is the prespecified strong comparator, GINE the primary graph comparator. Two RF configurations, two XGBoost configurations, four Ridge penalties and two small neural configurations are selected using validation family-macro MAE only. XGBoost uses validation early stopping, maximum 600 trees; neural models use at most 100 epochs and patience 15. No refitting on calibration or test data. Exact grids and hashes are in experiment code/manifests.\n\n## Scaffold + study excluded\n\n'+md(primary,columns)+'\n\n## All partition designs\n\n'+md(summary,columns)+'\n\nStudy-only clearance partitions are identical to scaffold+study and are executed once, not double-counted. Compound/scaffold-only diagnostics retain cross-partition publication overlap and cannot establish unseen-study performance. The sparse Cmax strict result has two admissible split seeds; do not interpret it as a robust ranking.\n\n'+source)
    neural=main[main.model.isin(['mlp','gcn','gin','gine'])]
    arch=neural.groupby('model').agg(runs=('run_id','size'),min_parameters=('parameters','min'),max_parameters=('parameters','max'),median_best_epoch=('best_epoch','median'),training_seconds=('seconds','sum')).reset_index()
    write('gnn_experiments.md','# Graph neural-network experiments\n\n'+scope+'\n\n'+findings+'\n\n## Executed architectures\n\n'+md(arch)+'\n\nPyTorch Geometric implements two-layer GCN, GIN and edge-aware GINE with hidden width 32 or 64, ReLU, concatenated mean/sum pooling, and a dose-conditioned prediction head. Sum pooling is divided by a fixed 20; this is not a fitted normalization. Atom features encode element, degree, formal charge, hybridization, CIP chirality, aromaticity, ring membership, hydrogen count and mass. Bond features encode bond type, stereo, conjugation and ring membership. GCN/GIN use connectivity but not bond attributes; GINE uses attributes. The explicit GINE ablation retains only atom type and connectivity with constant edge features.\n\nAll runs use CPU, two threads, AdamW, weight decay 0.001, dropout 0.1, gradient clipping and validation-only early stopping. Target centering/scaling and descriptor transforms use training observations only. Training losses weight repeated chemical families equally. Parameter counts, runtimes and learning histories are recorded. Tests verify atom-order invariance and that changing held-out labels cannot change fitted graph predictions.\n\n## Prespecified GINE versus XGBoost\n\n'+md(pairs[(pairs.model=='gine')&(pairs.split=='scaffold_study')])+ '\n\nDifferences are GINE minus XGBoost; negative favors GINE. Component-bootstrap intervals resample fixed paired errors over held-out linked study/identity components, only with at least eight groups. They do not refit models, are exploratory and unadjusted for multiple comparisons. Too-small splits have no interval. No graph-superiority claim follows from a single favorable score.\n\n## Atom/bond feature ablation\n\n'+md(all_metrics[(all_metrics.ablate)|(all_metrics.task.eq('rat_iv_plasma_clearance')&all_metrics.split.eq('scaffold_study')&all_metrics.seed.eq(42)&all_metrics.model.eq('gine')&all_metrics.train_fraction.eq(1)&all_metrics.member.eq(0))],['task','split','seed','model','ablate','mae_log10','rmse_log10','n_train_families'])+'\n\nThis one-seed comparison cannot establish a general feature effect. No pretrained encoder or multitask model was used. Small, heterogeneous endpoint cohorts do not justify increased network complexity merely to pursue positive results.\n\nImplementation references: [PyG GINEConv](https://pytorch-geometric.readthedocs.io/en/latest/generated/torch_geometric.nn.conv.GINEConv.html), [PyG installation](https://pytorch-geometric.readthedocs.io/en/2.7.0/install/installation.html).\n\n'+source)
    raw=pd.read_csv(ROOT/'data/processed/molecular/eligible_endpoints.csv',keep_default_na=False)
    meta=[]
    for endpoint,g in raw.groupby('endpoint'):
        selected=table[table.endpoint==endpoint]
        meta.append(dict(endpoint=endpoint,dataset_A_records=len(g),numeric_dose_field_records=int(pd.to_numeric(g.dose_mg_kg,errors='coerce').gt(0).sum()),
                         known_matrix_records=int(g.matrix.ne('').sum()),known_route_records=int(g.route.ne('').sum()),publication_field_records=int(g.publication_id.ne('').sum()),
                         ML_original_records=int(selected.original_records.sum()),ML_rows=len(selected),ML_compounds=selected.inchikey.nunique()))
    pd.DataFrame(meta).to_csv(RESULTS/'metadata_coverage.csv',index=False)
    audit_path=ROOT/'reports/ml_dataset_audit.md'
    base=audit_path.read_text().split('\n## Executed final coverage')[0]
    audit_path.write_text(base+'\n## Executed final coverage\n\n'+md(pd.DataFrame(meta))+'\n\n'+md(pd.DataFrame(audit['tasks']))+
                         f"\n\nUnion: {audit['compounds']} full-InChIKey compounds; {audit['records_eligible']} eligible original records; {len(table)} task rows. "
                         f"{audit['missing_stereochemistry_rows']} rows contain unspecified stereocenters, preserved as reported and linked conservatively for splitting. "
                         "Explicit method subtype (e.g. NCA versus compartment fitting) is not available for every source endpoint; these are reported PK labels, not independently reconstructed animal-level measurements. No assumptions about such missing method metadata are silently encoded.\n")
    split_frame=pd.DataFrame(splits)
    invalid=split_frame[~split_frame.eligible]
    learning=all_metrics[(all_metrics.task=='rat_iv_plasma_clearance')&(all_metrics.split=='scaffold_study')&(all_metrics.seed==42)&all_metrics.model.isin(['xgb','gine'])&(~all_metrics.ablate)&(all_metrics.member==0)]
    write('generalization_analysis.md','# Generalization analysis\n\n'+scope+'\n\n'+findings+'\n\n## Split constraints\n\nCompound families link full connectivity and canonical neutral tautomer identities irrespective of stereochemistry. Study components link original publications transitively, and scaffold components use stereo-free Bemis–Murcko scaffolds (all acyclic compounds grouped). No target values determine partitions. All scales, imputers and calibration functions are fit within permitted partitions.\n\n'+md(split_frame[['task','split','seed','eligible','partition_compounds','n_groups','largest_group_records','studies_crossing_partitions','scaffolds_crossing_partitions']])+'\n\n'+str(len(invalid))+' requested task/split seeds fail the minimum partition-size rule and were not trained. They were not replaced with more favorable random seeds. Random-compound and scaffold-only scores are diagnostics with explicit publication overlap, not leakage-free study generalization claims.\n\n## Generalization contrasts\n\n'+md(summary[summary.model.isin(['xgb','rf','gine'])],columns)+'\n\nThere is one eligible source (ChEMBL) and one species (rat). Source- and species-stratified scores therefore equal each task score; genuine cross-source or cross-species tests are unsupported. Differences across partition schemes mix changes in chemical distance, study composition and test difficulty; they do not isolate a causal scaffold effect.\n\n## Learning curves\n\n'+md(learning,['model','train_fraction','n_train_families','mae_log10','rmse_log10','validation_mae_log10'])+'\n\nLearning curves use nested subsets of whole training study/scaffold groups, with the same held-out validation/calibration/test sets. Each fraction is retuned on validation with the same declared grid. One split and two smaller fractions are exploratory; curves need not be monotone.\n\n## Applicability domain and transfer\n\nError histograms and the five largest absolute errors per strict model/split are saved in results/ml/figures/08_error_distributions and results/ml/outliers.csv. These post-hoc diagnostic selections never change task eligibility or refit models. Nearest-training Tanimoto similarity is calculated with cached chirality-aware Morgan fingerprints. Family-averaged error correlations and risk curves are descriptive, not proof of a universal similarity threshold. Source metadata features are not tested because source is constant and publication IDs are unavailable as meaningful features for a new study.\n\nDataset A→B transfer was not executed. The current molecular tasks are conditional CL/Cmax, not AUC dose-response representations; Dataset B has only 19 dependency components overall and 11/12 units in the established two-support comparison. A valid transfer test would need identity-disjoint representation training, disjoint model selection and calibration, plus enough target components to compare against linear scaling. This study does not establish that feasibility, and no transfer benefit is claimed.\n\n'+source)
    write('uncertainty_analysis.md','# Uncertainty and reliability\n\n'+scope+'\n\n'+findings+'\n\n## Calibration protocol\n\nFinite-sample construction follows [Angelopoulos and Bates](https://arxiv.org/abs/2107.07511), applied here to group-maximum residuals; exchangeability assumptions remain essential.\n\nA separate calibration partition was reserved before any fitting; it is never used for hyperparameters or early stopping. For study-excluded splits, nonconformity scores are the maximum absolute log10 residual within each calibration dependency group. At nominal coverage 1−α, radius q is the ceil((n+1)(1−α))-th ordered score. If that rank exceeds n, no finite interval is reported. Intervals are symmetric in log10 space and multiplicative in original units. Coverage is evaluated for **all observations in each held-out group**, alongside family coverage. This targets group coverage, not animal-level variability.\n\nUnder exchangeable calibration/test groups this is the usual split-conformal construction. Scaffold/source-like distribution shifts, unequal study designs and purposive public-data selection do not guarantee exchangeability. Accordingly, these results are empirical coverage assessments, not calibrated guarantees on unseen scaffolds. Compound/scaffold-only splits with shared publications are excluded from this group-calibration claim.\n\n## Three-member ensembles (fixed clearance scaffold/study split, seed 42)\n\n'+md(ensembles)+'\n\nGINE and XGBoost ensembles use the selected base-model hyperparameters with three random initializations/seeds on the same training partition. The calibration and test molecules remain excluded from every member. Variation captures optimization/subsampling variability, not all epistemic or measurement uncertainty. Ensemble size is small. Spearman associations and selective-risk curves are descriptive and may be unstable. No threshold was chosen using test errors.\n\n## Empirical 90% interval coverage\n\n'+md(interval[np.isclose(interval.nominal,.9)])+'\n\nVery small calibration/test group counts produce coarse estimates and wide intervals. The same test data appear at several nominal levels; these are not independent trials. No animal-level prediction intervals, conformal guarantees under covariate shift, or prospective dose-selection decisions are supported.\n\nRisk–coverage curves rank families by ensemble SD or by one minus nearest-training similarity. They report error among retained test families retrospectively; they are not validated deployment policies.\n\n'+source)
    primary_gine=pairs[(pairs.model=='gine')&(pairs.split=='scaffold_study')]
    summary_obj=dict(dataset=audit,completed_primary_runs=len(main),completed_runs=len(all_metrics),training_seconds=float(all_metrics.seconds.sum()),
                     first_event=events[0]['time'],last_event=events[-1]['time'],primary_summary=primary[columns].replace({np.nan:None}).to_dict('records'),
                     primary_gine_comparisons=primary_gine.replace({np.nan:None}).to_dict('records'),ensembles=ensembles.replace({np.nan:None}).to_dict('records'),
                     original_artifacts_preserved=all(sha(ROOT/p)==h for p,h in json.loads((ROOT/'reports/ml/original_artifacts.json').read_text()).items()),
                     figures=sorted(str(p.relative_to(ROOT)) for p in (RESULTS/'figures').glob('*.pdf')))
    dump(ROOT/'reports/ml/experiment_summary.json',summary_obj)
    write('final_ml_research_summary.md','# Final molecular PK research summary\n\n'+scope+'\n\n'+findings+'\n\n## Actual scope and research questions\n\n'+f"The conservative ML layer uses **{audit['compounds']} unique structure-resolved compounds**, {audit['records_eligible']} original eligible endpoint records and {len(table)} assay-specific publication/identity/dose cells. It does not train on all 2,713 Dataset A compounds. {len(main)} primary model runs and {len(all_metrics)-len(main)} focused follow-up runs completed. Actual first/last experiment events: {events[0]['time']} / {events[-1]['time']}; summed fit/tuning wall times {all_metrics.seconds.sum():.1f} seconds (overlapping processes make this different from elapsed session time). This was not 72 hours of execution.\n\n"+md(pd.DataFrame(audit['tasks']))+'\n\nRQ1 is assessed only for conditional rat plasma clearance/Cmax. AUC, half-life and volume are withheld for incompatible/insufficient definitions. RQ2 compares small graph networks with prespecified XGBoost and other baselines. RQ3 includes compound, scaffold and study-aware designs but cannot assess genuine cross-source generalization. RQ4 has a limited independent-calibration assessment. RQ5 transfer is unsupported and was not forced.\n\n## Strict scaffold + study results\n\n'+md(primary,columns)+'\n\n'+common+'\n\n## Evidence for or against graph improvement\n\n'+md(primary_gine)+'\n\nThe prespecified GINE comparison, rather than selecting a winning architecture after seeing test outcomes, is the basis for graph-versus-fingerprint interpretation. A negative difference favors GINE. Two-seed Cmax strict evaluation is particularly underpowered. No positive result from a single split is treated as established improvement. Full model rankings are descriptive, not selection for deployment.\n\n## Limitations and recommendations\n\nThe original sources contain far more endpoints than the scientifically compatible ML cohorts. Unknown original matrix/study metadata, unspecified AUC windows and half-life phases prevent training several requested targets. Retained source-reported assays still vary in sampling schedule, formulation, sex and strain; study exclusion reduces leakage but cannot remove all label heterogeneity. Unspecified stereocenters are preserved and grouped conservatively. Single-source rat cohorts do not establish broader species/source generalization. Small held-out study groups and overlapping split seeds limit precision. No calibrated-interval guarantee under scaffold shift, mechanistic saturation inference or Dataset B transfer claim is made.\n\nPrefer additional independent, metadata-complete studies and prospective validation over larger networks. Keep the known dose input and a dose-only control when measuring the incremental value of molecular features. Investigate failures with profile-level and experimental variance information.\n\n## Deliverables and reproducibility\n\nSee [dataset audit](ml_dataset_audit.md), [full benchmark](molecular_ml_benchmark.md), [GNN experiments](gnn_experiments.md), [generalization](generalization_analysis.md), [uncertainty](uncertainty_analysis.md) and [execution log](experiment_log.md). Exact split/task tables, local model checkpoints, prediction rows and metrics are under data/processed/ml and results/ml. Figures are vector PDF plus 300-dpi PNG under results/ml/figures. The reproducible launcher uses no paid hardware; caches permit offline experiment replay.\n\n'+source)


def main():
    for p,h in json.loads((ROOT/'reports/ml/original_artifacts.json').read_text()).items():assert sha(ROOT/p)==h,'Original artifact changed: '+p
    table=pd.read_csv(OUT/'tasks.csv',keep_default_na=False);table['dependency_group']=''
    for _,g in table.groupby('task'):table.loc[g.index,'dependency_group']=connected_groups(g,'scaffold_study')
    all_metrics,main_metrics,summary=primary_metrics();pairs=paired(main_metrics,table);interval,risk,ensembles=reliability(main_metrics)
    figures(main_metrics,summary,all_metrics,interval,risk);reports(all_metrics,main_metrics,summary,pairs,interval,ensembles)
    aggregates=ROOT/'reports/ml/tables';aggregates.mkdir(exist_ok=True)
    for name in ['model_summary','paired_comparisons','intervals','ensemble_summary','similarity_associations','metadata_coverage','risk_coverage']:
        shutil.copyfile(RESULTS/(name+'.csv'),aggregates/(name+'.csv'))
    print('Generated reports and figures from',len(all_metrics),'completed runs; original A/B preserved')

if __name__=='__main__':main()
