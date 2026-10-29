"""Final presentation claims must remain tied to actual rows and predictions."""
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('verify_final_results',ROOT/'scripts/verify_final_results.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_final_counts_and_metrics_reconcile():
    v=module.verify()
    assert v['dataset_a']['unique_structured_compounds'] > v['dataset_b']['compounds']
    assert [r['n_components'] for r in v['baseline_results']]==[11,11,12,12]


def test_final_verification_recalculates_losses_and_weights_units():
    # A stale saved loss column and many repeats must not inflate precision/weight.
    d=pd.DataFrame([dict(component_id='a',endpoint_series_id='a',predicted=2.,observed=1.,absolute_log_error=999.)]*20+
                   [dict(component_id='b',endpoint_series_id='b',predicted=1.,observed=1.,absolute_log_error=999.)])
    r=module.macro_errors(d)
    assert np.isclose(r['mae_ln'],np.log(2)/2)
    assert np.isclose(r['mean_absolute_relative_error'],.5)


def test_research_brief_has_two_pages_and_verified_primary_table():
    from pypdf import PdfReader
    pdf=PdfReader(ROOT/'reports/pk_tk_research_brief.pdf')
    assert len(pdf.pages)==2
    first=pdf.pages[0].extract_text();second=pdf.pages[1].extract_text()
    assert '2,713' in first and '6,459' in first and '59 comparable series' in first
    for row in module.verify()['baseline_results']:
        assert f"{row['mae_ln']:.3f}" in second
        assert f"{row['rmse_ln']:.3f}" in second
        assert f"{100*row['mean_absolute_relative_error']:.1f}%" in second
    assert 'Both include zero' in second
