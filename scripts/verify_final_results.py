"""Independently reconcile final claims with local endpoint rows and predictions."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def macro_errors(predictions):
    """Recalculate from predictions, never reuse saved loss columns."""
    d = predictions.copy()
    ratio = d.predicted / d.observed
    d['mae_ln'] = np.abs(np.log(ratio))
    d['mse_ln'] = np.log(ratio)**2
    d['mean_absolute_relative_error'] = np.abs(ratio - 1)
    cols = ['mae_ln', 'mse_ln', 'mean_absolute_relative_error']
    m = d.groupby(['component_id', 'endpoint_series_id'])[cols].mean().groupby('component_id').mean().mean()
    return dict(mae_ln=float(m.mae_ln), rmse_ln=float(np.sqrt(m.mse_ln)),
                mean_absolute_relative_error=float(m.mean_absolute_relative_error))


def verify(root=ROOT):
    inputs = {}
    def read(name):
        p = root/name
        inputs[name] = hashlib.sha256(p.read_bytes()).hexdigest()
        return pd.read_csv(p)
    a = read('data/processed/molecular/eligible_endpoints.csv')
    all_a = read('data/processed/molecular/endpoint_records.csv')
    compounds = read('data/processed/molecular/compounds.csv')
    s = read('data/processed/crossdose_series_expanded.csv')
    s = s[s.primary_eligible]
    t = read('results/research/endpoint_observations.csv')
    ep = read('results/research/episodes.csv')
    p = read('results/research/predictions.csv')
    metrics = read('results/research/metrics.csv')
    pairs = read('results/research/paired_comparisons.csv')
    expected = json.loads((root/'reports/molecular_summary.json').read_text())
    sets = {k:set(g.inchikey) for k,g in a.groupby('endpoint')}
    ac = dict(unique_structured_compounds=a.inchikey.nunique(), eligible_endpoint_records=len(a),
              excluded_endpoint_records=int((~all_a.eligible_A).sum()), collected_endpoint_records=len(all_a),
              compounds_with_auc=len(sets['auc']), compounds_with_cmax=len(sets['cmax']),
              compounds_with_both_auc_cmax=len(sets['auc'] & sets['cmax']),
              compounds_with_clearance=len(sets['clearance']), compounds_with_volume=len(sets['volume_distribution']),
              compounds_with_half_life=len(sets['half_life']), known_source_publications=a.publication_id.nunique(),
              dose_normalized_records=int(a.dose_normalized.eq('yes').sum()))
    assert set(compounds.inchikey)==set(a.inchikey) and compounds.inchikey.is_unique
    assert a.canonical_smiles.notna().all() and a.inchikey.notna().all()
    assert len(a) == int(all_a.eligible_A.sum())
    bc = dict(validated_series=s.series_id.nunique(), compounds=s.molecular_group_id.nunique(),
              independent_studies=s.study_id.nunique(), auc_series=int(s.eligible_auc.sum()),
              cmax_series=int(s.eligible_cmax.sum()), both_endpoint_series=int((s.eligible_auc & s.eligible_cmax).sum()),
              compounds_two_comparable_doses=s[s.n_doses>=2].molecular_group_id.nunique(),
              compounds_three_comparable_doses=s[s.n_doses>=3].molecular_group_id.nunique())
    for label, actual in [('dataset_a',ac),('dataset_b',bc)]:
        for key, value in actual.items():
            assert value == expected[label][key], (label,key,value,expected[label][key])
    assert ep.episode_id.is_unique
    assert set(ep.episode_id) == set(p.episode_id)
    assert len(ep)==513 and t.component_id.nunique()==19
    main = p[(p.status=='ok') & (p.policy=='leave_one_dose_out') & (p.task=='low_to_high') &
             (p.n_support==2) & p.model.isin(['linear','power'])]
    results=[]
    for (endpoint,model),g in main.groupby(['endpoint','model']):
        r=dict(endpoint=endpoint,model=model,n_components=g.component_id.nunique(),n_series=g.endpoint_series_id.nunique(),
               n_compounds=g.molecular_group_id.nunique(), **macro_errors(g))
        saved=metrics[(metrics.endpoint==endpoint)&(metrics.model==model)&(metrics.aggregation=='overall')&
                      (metrics.task=='low_to_high')&(metrics.n_support==2)&(metrics.policy=='leave_one_dose_out')]
        assert len(saved)==1
        for k in ['mae_ln','rmse_ln','mean_absolute_relative_error','n_components']:
            np.testing.assert_allclose(r[k],saved.iloc[0][k],rtol=1e-9,atol=1e-10)
        results.append(r)
    paired=pairs[(pairs.model=='power')&(pairs.task=='low_to_high')&(pairs.n_support==2)&(pairs.policy=='leave_one_dose_out')]
    for _,r in paired.iterrows():
        by={v['model']:v for v in results if v['endpoint']==r.endpoint}
        np.testing.assert_allclose(by['power']['mae_ln']-by['linear']['mae_ln'],r.delta_mae_ln,atol=1e-10)
        assert r.ci_low < 0 < r.ci_high
    return dict(dataset_a=ac,dataset_b=bc,episodes=len(ep),independent_components=t.component_id.nunique(),
                successful_predictions=int(p.status.eq('ok').sum()),baseline_results=results,
                paired_comparisons=paired.to_dict('records'),input_sha256=inputs,
                verification='Counts recalculated from endpoint rows; losses recalculated from predictions with equal series/component weighting.')

if __name__=='__main__':
    result=verify()
    (ROOT/'reports/final_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='input_sha256'},indent=2))
