import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
pytest.importorskip('torch')
pytest.importorskip('torch_geometric')
from pktk.ml.data import ROOT,OUT,structure,connected_groups,allocate,sha,exclusion
from pktk.ml.features import graph,nearest_similarity
from pktk.ml.neural import GraphModel,fit_neural,seed_all
import torch
from torch_geometric.data import Batch


def test_original_data_preservation():
    snapshot=json.loads((ROOT/'reports/ml/original_artifacts.json').read_text())
    for path,expected in snapshot.items():assert sha(ROOT/path)==expected,path


def test_split_families_preserve_stereo_and_tautomers():
    assert structure('C[C@H](O)C(=O)O')[0]==structure('C[C@@H](O)C(=O)O')[0]
    assert structure('CC(=O)C')[0]==structure('CC(O)=C')[0]
    assert structure('CCC')[1]=='__ACYCLIC__'
    with pytest.raises(ValueError):structure('CC.[Na+]')


def test_actual_task_gates_and_split_leakage():
    d=pd.read_csv(OUT/'tasks.csv',keep_default_na=False)
    assert d.species.eq('rat').all() and d.matrix.eq('plasma').all()
    assert d.dose_mg_kg.gt(0).all() and set(d.endpoint)=={'clearance','cmax'}
    assert d.source.eq('chembl').all() and d.publication_id.ne('').all()
    s=pd.read_csv(OUT/'splits.csv',keep_default_na=False)
    for (_,kind,_),g in s.groupby(['task','split','seed']):
        for col in ['family','connectivity_key']+(['publication_id'] if 'study' in kind else [])+(['ml_scaffold'] if 'scaffold' in kind else []):
            assert g.groupby(col).partition.nunique().max()==1


def test_splits_ignore_targets_and_transitive_links():
    d=pd.DataFrame(dict(family=['a','b','b','c'],connectivity_key=['a','b','b','c'],publication_id=['x','x','y','y'],ml_scaffold=['s1','s2','s2','s3']))
    assert connected_groups(d,'study').nunique()==1
    assert connected_groups(d,'compound').nunique()==3
    actual=pd.read_csv(OUT/'tasks.csv',keep_default_na=False).query("endpoint == 'clearance'")
    cfg=json.loads((ROOT/'configs/ml/benchmark.json').read_text())
    s,_=allocate(actual,'scaffold_study',42,cfg)
    changed=actual.copy();changed['y']=np.random.default_rng(1).normal(size=len(actual))*100
    t,_=allocate(changed,'scaffold_study',42,cfg)
    pd.testing.assert_frame_equal(s,t)


def test_graph_features_and_permutation_invariance():
    g=graph('C[C@H](O)C(=O)O');other=graph('C[C@@H](O)C(=O)O')
    assert not torch.equal(g.x,other.x)
    seed_all(1);model=GraphModel('gine',g.num_node_features,g.num_edge_features).eval()
    g.dose=torch.tensor([0.]);perm=torch.arange(g.num_nodes-1,-1,-1);inv=torch.argsort(perm)
    h=g.clone();h.x=g.x[perm];h.edge_index=inv[g.edge_index]
    with torch.no_grad():a=model(Batch.from_data_list([g]));b=model(Batch.from_data_list([h]))
    torch.testing.assert_close(a,b,atol=1e-6,rtol=1e-5)


def test_heldout_exposure_cannot_change_graph_fit(tmp_path):
    n=16;gs=[graph('C'*(i%5+2)) for i in range(n)]
    t=pd.DataFrame(dict(family=[str(i) for i in range(n)],y=np.linspace(-1,1,n),log_dose=np.zeros(n)))
    idx=dict(train=np.arange(8),validation=np.arange(8,11),calibration=np.arange(11,13),test=np.arange(13,16))
    cfg=dict(threads=1,batch_size=4,max_epochs=3,patience=2)
    args=('gine',gs,np.zeros((n,1)),t,idx,dict(hidden=8,lr=.001),12,cfg)
    a,_=fit_neural(*args,tmp_path/'a.pt')
    t.loc[11:,'y']=1000
    b,_=fit_neural(*args,tmp_path/'b.pt')
    np.testing.assert_array_equal(a['test'],b['test'])


def test_similarity_uses_training_only():
    fps=np.array([[1,0,1],[0,1,0],[1,0,1]],dtype=float)
    assert nearest_similarity(fps,np.array([0]),np.array([1,2])).tolist()==[0.,1.]
