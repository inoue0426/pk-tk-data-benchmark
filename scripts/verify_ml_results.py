"""Reconcile every executed run with saved splits, predictions and validation selection."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]


def verify():
    cfg=json.loads((ROOT/'configs/ml/benchmark.json').read_text())
    audit=json.loads((ROOT/'data/processed/ml/split_audit.json').read_text())
    table=pd.read_csv(ROOT/'data/processed/ml/tasks.csv').set_index('row_id')
    splits=pd.read_csv(ROOT/'data/processed/ml/splits.csv')
    expected={(m['task'],m['split'],m['seed'],model) for m in audit if m['eligible'] and not (m['task']=='rat_iv_plasma_clearance' and m['split']=='study') for model in cfg['models']}
    found=set();checked=[];maximum_difference=0.
    for path in sorted((ROOT/'results/ml/runs').glob('*/metrics.json')):
        r=json.loads(path.read_text());p=pd.read_csv(path.parent/'predictions.csv')
        assert p.row_id.is_unique
        split=splits[(splits.task==r['task'])&(splits.split==r['split'])&(splits.seed==r['seed'])].set_index('row_id')
        assert p.set_index('row_id').partition.to_dict()==split[split.partition!='train'].partition.to_dict()
        train=r['train_row_ids'];assert set(train)<=set(split[split.partition=='train'].index)
        assert set(table.loc[train].family).isdisjoint(table.loc[p.row_id].family)
        if 'study' in r['split']:assert set(table.loc[train].publication_id).isdisjoint(table.loc[p.row_id].publication_id)
        if 'scaffold' in r['split']:assert set(table.loc[train].ml_scaffold).isdisjoint(table.loc[p.row_id].ml_scaffold)
        np.testing.assert_allclose(p.y,table.loc[p.row_id].y,atol=1e-10,rtol=1e-10)
        test=p[p.partition=='test'].copy();test['ae']=abs(test.prediction_log10-test.y);test['se']=(test.prediction_log10-test.y)**2
        mae=float(test.groupby('family').ae.mean().mean());rmse=float(np.sqrt(test.groupby('family').se.mean().mean()))
        np.testing.assert_allclose([mae,rmse],[r['mae_log10'],r['rmse_log10']],atol=1e-9,rtol=1e-9)
        maximum_difference=max(maximum_difference,abs(mae-r['mae_log10']),abs(rmse-r['rmse_log10']))
        tuning=json.loads((path.parent/'tuning.json').read_text())
        assert r['selected_candidate']==min(tuning,key=lambda c:c['validation_mae_log10'])['candidate']
        val=p[p.partition=='validation'].copy();val['loss']=abs(val.prediction_log10-val.y)
        np.testing.assert_allclose(val.groupby('family').loss.mean().mean(),r['validation_mae_log10'],atol=1e-8,rtol=1e-8)
        if r['train_fraction']==1 and not r['ablate'] and r['member']==0:found.add((r['task'],r['split'],r['seed'],r['model']))
        checked.append(r['run_id'])
    assert expected==found,(expected-found,found-expected)
    assert len(checked)==len(expected)+9,'Expected five ablation and four ensemble-member runs'
    result=dict(primary_runs=len(expected),focused_followups=9,total_runs_checked=len(checked),maximum_metric_difference=maximum_difference,
                all_selected_candidates_match_validation_minimum=True,all_split_memberships_verified=True,
                requested_split_cases=len(audit),withheld_split_cases=sum(not m['eligible'] for m in audit),
                source_count=table.source.nunique(),species_count=table.species.nunique(),
                scope='All saved run scores and memberships checked; this is not a fresh refit of every model.')
    (ROOT/'reports/ml/results_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':verify()
