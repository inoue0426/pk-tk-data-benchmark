"""Rebuild eligibility/splits and independently refit fixed baseline/GINE experiments."""
import os,sys,importlib.util
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
if sys.platform=='darwin' and os.environ.get('PKTK_ML_VERIFY_STARTED')!='1':
    lib=Path(importlib.util.find_spec('torch').origin).parent/'lib'
    env=dict(os.environ,PKTK_ML_VERIFY_STARTED='1',DYLD_LIBRARY_PATH=str(lib),PYTHONPATH=str(ROOT/'src'),OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2')
    os.execve(sys.executable,[sys.executable,__file__],env)
import hashlib,json,subprocess,tempfile
import numpy as np
import pandas as pd
from pktk.ml import data,run
from pktk.ml.features import load_features


def main():
    artifacts=sorted(data.OUT.glob('*'))
    before={str(p.relative_to(ROOT)):data.sha(p) for p in artifacts if p.is_file()}
    subprocess.run([sys.executable,'-m','pktk.ml.data'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
    after={p:data.sha(ROOT/p) for p in before}
    assert before==after,'Eligibility/split replay differs'
    original=json.loads((ROOT/'reports/ml/original_artifacts.json').read_text())
    assert all(data.sha(ROOT/p)==h for p,h in original.items())
    cfg=json.loads((ROOT/'configs/ml/benchmark.json').read_text())
    t=pd.read_csv(data.OUT/'tasks.csv',keep_default_na=False);s=pd.read_csv(data.OUT/'splits.csv',keep_default_na=False)
    desc,fps,graphs=load_features(t)
    task='rat_iv_plasma_clearance';kind='scaffold_study';seed=42
    split=s[(s.task==task)&(s.split==kind)&(s.seed==seed)]
    # Fresh output directory guarantees a real refit, never checkpoint reuse.
    run.RESULTS=Path(tempfile.mkdtemp(prefix='replay_',dir=ROOT/'results/ml'))
    fits=[]
    for model in ['xgb','gine']:
        info=run.experiment(t,desc,fps,graphs,split,task,kind,seed,model,cfg)
        old=ROOT/'results/ml/runs'/info['run_id'];new=run.RESULTS/'runs'/info['run_id']
        a=pd.read_csv(old/'predictions.csv');b=pd.read_csv(new/'predictions.csv')
        assert a.row_id.tolist()==b.row_id.tolist()
        np.testing.assert_allclose(a.prediction_log10,b.prediction_log10,rtol=1e-9,atol=1e-10)
        fits.append(dict(model=model,refit=True,predictions_byte_identical=data.sha(old/'predictions.csv')==data.sha(new/'predictions.csv'),max_absolute_prediction_difference=float(abs(a.prediction_log10-b.prediction_log10).max())))
    # Reports/figures must be deterministic, while timestamps and training seconds are preserved logs.
    subprocess.run([sys.executable,'-m','pktk.ml.analyze'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
    files=list((ROOT/'results/ml/figures').glob('*'))+[ROOT/'reports'/name for name in ['molecular_ml_benchmark.md','gnn_experiments.md','generalization_analysis.md','uncertainty_analysis.md','final_ml_research_summary.md']]
    files+=list((ROOT/'reports/ml/figures').glob('*.pdf'))+list((ROOT/'reports/ml/tables').glob('*.csv'))
    hashes={str(p.relative_to(ROOT)):data.sha(p) for p in files}
    subprocess.run([sys.executable,'-m','pktk.ml.analyze'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
    assert all(data.sha(ROOT/p)==h for p,h in hashes.items()),'Analysis rendering replay differs'
    result=dict(task_artifacts_identical=len(before),original_artifacts_preserved=len(original),fresh_refits=fits,analysis_artifacts_identical=len(files),
                input_sha256=after,analysis_sha256=hashes,scope='Fresh refit verification covers one fixed XGBoost and one fixed GINE study/scaffold split. All other runs have resumable manifests; no claim of fresh refitting every run.')
    data.dump(ROOT/'reports/ml/reproducibility.json',result)
    print(json.dumps({k:v for k,v in result.items() if 'sha256' not in k},indent=2))

if __name__=='__main__':main()
