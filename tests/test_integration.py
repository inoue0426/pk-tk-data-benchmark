import json
from pathlib import Path
import pandas as pd
import pytest

@pytest.mark.skipif(not Path('reports/summary.json').exists(),reason='Run pipeline to create real-data outputs')
def test_actual_pipeline_outputs():
    summary=json.loads(Path('reports/summary.json').read_text())
    assert summary['cvtdb_original_observations']==35669
    obs=pd.read_csv('data/processed/observations.csv',low_memory=False)
    assert len(obs)==35669
    assert obs.original_json.notna().all()
    high=pd.read_csv('data/processed/high_confidence_cmax.csv')
    assert len(high)==6 and high.compound_id.nunique()==1
    splits=pd.read_csv('results/task_a_splits.csv')
    for _,g in splits.groupby('episode_id'):
        a=g[g.role=='support'];b=g[g.role=='target']
        assert set(a.experiment_id).isdisjoint(b.experiment_id)
        assert a.dose_mg_kg.max()<b.dose_mg_kg.min()
    predictions=pd.read_csv('results/predictions.csv')
    assert (predictions.predicted>0).all()
    status=json.loads(Path('results/benchmark_status.json').read_text())
    assert status['B'].startswith('not_generated')
    assert status['C'].startswith('not_generated')
