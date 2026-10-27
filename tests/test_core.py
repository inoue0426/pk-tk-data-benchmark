import math
import numpy as np
import pandas as pd
import pytest
from pktk.standardization.core import convert,structure
from pktk.models.scaling import linear_scaling,power_law
from pktk.benchmark.evaluate import check_split
from pktk.benchmark.series import curves_and_candidates,GROUP_FIELDS

def test_units_are_dimension_safe():
    assert convert(1000,'ng/mL','concentration')==1
    assert convert(60,'min','time')==1
    assert convert(1,'g/kg','dose')==1000
    assert math.isnan(convert(20,'mg','dose'))
    assert math.isnan(convert(3,'uM','concentration'))
    assert math.isnan(convert('<LOQ','ng/mL','concentration'))

def test_rdkit_identity_and_invalid_structure():
    assert structure('OCC')['inchikey']==structure('CCO')['inchikey']
    assert structure('bad smiles')['structure_flag']=='invalid_structure'
    assert structure('')['structure_flag']=='missing_structure'
    assert '.' in structure('CC(=O)[O-].[Na+]')['canonical_smiles']

def test_model_identifiability_and_extrapolation():
    assert linear_scaling([1],[2],4)==pytest.approx(8)
    assert power_law([1,2],[3,12],4)==pytest.approx(48)
    with pytest.raises(ValueError): power_law([1],[2],4)
    with pytest.raises(ValueError): linear_scaling([2],[4],1)
    with pytest.raises(ValueError): linear_scaling([1],[0],2)

def test_split_rejects_matrix_and_dose_leakage():
    base={c:'x' for c in ['compound_id','species','strain','sex','route','matrix','study_id','protocol_id']}
    a=pd.DataFrame([{**base,'experiment_id':'a','dose_mg_kg':1.}])
    b=pd.DataFrame([{**base,'experiment_id':'b','dose_mg_kg':2.}])
    check_split(a,b)
    for key,val in [('experiment_id','a'),('matrix','serum'),('dose_mg_kg',.5)]:
        bad=b.copy();bad.loc[0,key]=val
        with pytest.raises(AssertionError): check_split(a,bad)

def test_partial_auc_no_zero_time_extrapolation_and_duplicates():
    base={c:'x' for c in GROUP_FIELDS}
    base.update(experiment_id='curve',study_id='study',dose_mg_kg=2,n_animals=3,eligible_for_curve=True)
    x=pd.DataFrame([{**base,'time_h':t,'concentration_ug_ml':c,'source_record_id':str(t)} for t,c in [(1,2),(2,4),(3,2)]])
    curves,_=curves_and_candidates(x)
    assert curves.iloc[0].auc==6
    assert curves.iloc[0].time_start_h==1
    bad=pd.concat([x,x.iloc[[0]]],ignore_index=True)
    curves,_=curves_and_candidates(bad)
    assert not curves.iloc[0].usable_curve
    assert 'duplicate_sampling_times' in curves.iloc[0].qc_flags

def test_qc_preserves_bad_and_duplicate_records():
    from pktk.qc.checks import check_observations
    base=dict(source='test',source_record_id='1',experiment_id='a',compound_id='x',species='human',
              sex='',matrix='plasma',route='unclear',dose_mg_kg=-1.,study_id='s',document_id='d',subject_id='u',
              qc_flags='',parent_analyte=False,concentration_ug_ml=float('nan'),time_h=float('nan'),
              regimen='unknown',radiolabeled='',publication_doi='')
    x=check_observations(pd.DataFrame([base,base]))
    assert len(x)==2
    assert not x.eligible_for_curve.any()
    for f in ['duplicate_record','invalid_time','invalid_dose','nonparent_analyte','not_preclinical','ambiguous_route']:
        assert x.qc_flags.str.contains(f).all()
