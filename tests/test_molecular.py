import json
from pathlib import Path
import pandas as pd
import pytest
from pktk.standardization.molecular import context,structure,endpoint
from pktk.benchmark.molecular_feasibility import linked_components
from pktk.molecular import pksmart_records
ROOT=Path(__file__).resolve().parents[1]

def test_dose_normalized_is_not_dose_specific():
    x=context('Cmax in rat blood normalized to a 1 mg/kg dose after oral administration')
    assert x['dose_normalized']=='yes' and x['dose_mg_kg'] is None
    assert context('AUC in rat after oral administration of 1 and 10 mg/kg')['dose_mg_kg'] is None
    assert context('AUC in rat plasma after a single oral dose of 10 mg/kg')['dose_mg_kg']==10
    assert context('Area under curve (0-6 hr) after 10 mg/kg oral administration in rat')['dose_mg_kg']==10
    assert context('Cmax in rat at dose 8.98 (mg/kg iv)')['dose_mg_kg']==8.98

def test_in_vitro_and_pd_are_not_in_vivo_pk():
    for d in ['Intrinsic clearance in rat hepatocytes in vitro',
              'Microsomal clearance after intravenous administration in rats',
              'Diuresis half-time in rat after oral dosing',
              'Area under the effect time curve after 10 mg/kg oral dosing in rat',
              'Predicted Cmax after oral administration in rat']:
        assert not context(d)['in_vivo_evidence']

def test_ambiguous_context_not_forced():
    x=context('AUC after oral and intravenous dosing in rat plasma and brain')
    assert x['route']=='' and x['matrix']==''
    assert 'multiple_routes' in x['context_flags']

def test_structures_preserve_stereo_and_reject_mixtures():
    a=structure('C[C@H](O)C(=O)O');b=structure('C[C@@H](O)C(=O)O')
    assert a['inchikey']!=b['inchikey'] and a['connectivity_key']==b['connectivity_key']
    assert structure('CC.[Na+]')['structure_status']=='multicomponent_structure_unresolved'
    assert structure('not a smiles')['structure_status']=='invalid_or_missing_structure'
    assert structure('CC(=O)[O-]')['inchikey']==structure('CC(=O)O')['inchikey']

def test_units_and_nonpositive_endpoint():
    n,_=endpoint('AUC','ng.hr.mL-1',1000)
    assert n['value']==1 and n['unit']=='ug*h/mL'
    n,_=endpoint('Cmax','nM',1000,200)
    assert n['value']==pytest.approx(.2)
    assert endpoint('AUC','mmHg hr-1',100)[0] is None
    assert endpoint('CL','mL.min-1.kg-1',0)[0] is None
    assert endpoint('T1/2','hr',float('inf'))[0] is None

def test_pksmart_measured_not_predicted_or_log_scaled():
    rows,_=pksmart_records(ROOT,True)
    r=[r for r in rows if r['compound_name']=='5-Fluorouracil' and r['species']=='rat' and r['endpoint']=='clearance'][0]
    assert r['value']==33 and r['dose_mg_kg'] is None
    assert all(r['source_record_id'].startswith('pksmart:S2:') and not r['strict_B_eligible'] for r in rows)
    assert not any(r['species']=='human' for r in rows)

def test_publication_compound_linkage_is_transitive():
    x=pd.DataFrame([dict(connectivity_key='a',publication_id='1'),dict(connectivity_key='b',publication_id='1'),
      dict(connectivity_key='b',publication_id='2'),dict(connectivity_key='c',publication_id='2'),dict(connectivity_key='d',publication_id='3')])
    m=linked_components(x);assert m['a']==m['b']==m['c'] and m['a']!=m['d']

def test_real_dataset_a_and_b_boundary():
    path=ROOT/'data/processed/molecular/endpoint_records.csv'
    if not path.exists():pytest.skip('Run molecular collection first')
    x=pd.read_csv(path,keep_default_na=False);good=x[x.eligible_A]
    assert not x.source_record_id.duplicated().any()
    assert good.inchikey.ne('').all()
    assert not x.loc[x.source.ne('strict_b'),'strict_B_eligible'].any()
    c=good[good.source.eq('chembl')]
    assert c.species.eq('rat').all() and c.in_vivo_evidence.astype(str).eq('True').all()
    assert pd.to_numeric(good.value).gt(0).all()
    assert c[c.dose_normalized.eq('yes')].dose_mg_kg.eq('').all()
    b=json.loads((ROOT/'reports/summary.json').read_text())
    assert b['primary_comparable_series']==59 and b['primary_comparable_compound_identifiers']==21
