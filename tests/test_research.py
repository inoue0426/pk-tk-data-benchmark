import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from pktk.models.empirical import predict,NotIdentifiable
from pktk.benchmark.episodes import components,endpoint_table,generate,audit_episode
from pktk.benchmark.research import predict_episodes,summarize,paired_comparisons

ROOT=Path(__file__).resolve().parents[1]

@pytest.mark.parametrize('target',[.5,3,16])
def test_empirical_known_power_and_units(target):
    d=np.array([1.,2.,8.]);y=3*d**1.5
    pred,p=predict('power',d,y,target)
    assert pred==pytest.approx(3*target**1.5)
    assert p['b']==pytest.approx(1.5)
    scaled,_=predict('power',d,y*1000,target)
    assert scaled==pytest.approx(pred*1000)
    linear,_=predict('linear',d,4*d,target)
    assert linear==pytest.approx(4*target)

def test_curvature_requires_identifiability():
    with pytest.raises(NotIdentifiable):predict('saturating',[1,2],[1,2],4)
    with pytest.raises(NotIdentifiable):predict('saturating',[1,2,4],[2,4,8],8)
    d=np.array([1.,4.,9.]);y=5*d/(1+d/3)
    value,p=predict('saturating',d,y,16)
    assert value==pytest.approx(5*16/(1+16/3),rel=1e-6)
    assert p['k']==pytest.approx(3,rel=1e-6)
    with pytest.raises(ValueError):predict('power',[1,2],[1,2],2)
    with pytest.raises(ValueError):predict('linear',[1,2],[0,2],3)

def real_table():
    e=pd.read_csv(ROOT/'data/processed/high_confidence_endpoints.csv',keep_default_na=False)
    s=pd.read_csv(ROOT/'data/processed/crossdose_series_expanded.csv',keep_default_na=False)
    return endpoint_table(e,s)

def test_real_episode_membership_and_target_invariance():
    t=real_table();eps,members=generate(t)
    assert len(t.endpoint_series_id.unique())==73
    assert not members.duplicated(['episode_id','source_record_id']).any()
    assert members[members.role=='target'].groupby('episode_id').size().eq(1).all()
    cfg=json.loads((ROOT/'configs/analysis/benchmark.json').read_text())
    ep=eps[(eps.task=='low_to_high') & (eps.n_support==2)].head(1)
    before,_=predict_episodes(t,ep,cfg)
    changed=t.copy();mask=(changed.endpoint_series_id==ep.iloc[0].endpoint_series_id)&(changed.source_record_id==ep.iloc[0].target_record_id)
    changed.loc[mask,'exposure']*=10000
    after,_=predict_episodes(changed,ep,cfg)
    np.testing.assert_allclose(before.predicted,after.predicted,equal_nan=True)
    # Direct leakage must fail rather than merely be logged.
    row=ep.iloc[0];g=t[t.endpoint_series_id==row.endpoint_series_id]
    with pytest.raises(ValueError,match='Target leakage'):audit_episode(row,g,t.iloc[:0])

def test_compound_study_transitive_closure():
    x=pd.DataFrame([dict(molecular_group_id='a',study_id='s1',compound_id='a'),
        dict(molecular_group_id='b',study_id='s1',compound_id='b'),
        dict(molecular_group_id='b',study_id='s2',compound_id='b'),
        dict(molecular_group_id='c',study_id='s2',compound_id='c'),
        dict(molecular_group_id='d',study_id='s3',compound_id='d')])
    m=components(x)
    assert m['a']==m['b']==m['c'] and m['a']!=m['d']

def test_macro_weighting_does_not_count_pairs_as_replicates():
    records=[]
    for c,n,loss in [('a',100,1.),('b',1,3.)]:
        for i in range(n):records.append(dict(component_id=c,endpoint_series_id=c,molecular_group_id=c,study_id=c,
            absolute_log_error=loss,squared_log_error=loss**2,absolute_relative_error=loss,within_twofold=0.))
    s=summarize(pd.DataFrame(records))
    assert s['mae_ln']==2 and s['rmse_ln']==pytest.approx(np.sqrt(5)) and s['n_components']==2

def test_no_incompatible_auc_pooling():
    t=real_table();eps,_=generate(t);row=eps[eps.endpoint=='auc'].iloc[0]
    support=t[(t.endpoint_series_id==row.endpoint_series_id)&t.source_record_id.isin(json.loads(row.support_ids_json))]
    ref=t[t.component_id!=row.component_id].head(1).copy();ref['exposure_definition']='AUC0-999h'
    with pytest.raises(ValueError,match='Incompatible'):audit_episode(row,support,ref)
