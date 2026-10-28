import json
from pathlib import Path
import pytest
from pktk.collectors.cebs_reviewed import collect,numeric
from pktk.collectors.http import Collector
from pktk.expansion import build_series

ROOT=Path(__file__).resolve().parents[1]

def test_reviewed_cebs_protocol_separation():
    x=collect(Collector(ROOT,offline=True));x['molecular_group_id']=x.compound_id
    s=build_series(x)
    assert len(s)==33 and s.primary_eligible.all()
    assert x.casrn.str.strip().eq(x.casrn).all()
    assert set(x.loc[x.study_id=='CEBS:S0558','casrn'])=={'84-65-1'}
    b=x[x.study_id=='CEBS:K05097']
    assert len(build_series(b))==8
    assert b.vehicle.nunique()==2
    assert b.high_confidence_auc.all() and b.high_confidence_cmax.all()
    ftoh=x[x.study_id=='CEBS:K04003B']
    assert len(ftoh)==6 and ftoh.high_confidence_cmax.all() and not ftoh.high_confidence_auc.any()
    assert ftoh.compound_id.eq('DTXSID7029904').all()
    assert ftoh.cmax.iloc[0]==pytest.approx(.320)
    assert all(json.loads(v)['ANALYTE']=='Fluorotelomer alcohol 8+2' for v in ftoh.original_json)
    # First rat oral NCA group: source ng*min/mL -> ug*h/mL, ng/mL -> ug/mL.
    r=b[b.source_record_id=='CEBS:K05097:Sheet1:row6'].iloc[0]
    assert r.auc==pytest.approx(28700/60000)
    assert r.cmax==pytest.approx(.9132)
    # Never relabel whole blood as plasma or animal PK means as composites.
    n=x[x.study_id=='CEBS:S0571']
    assert set(n.matrix)=={'blood'}
    assert set(n[n.species=='rat'].endpoint_origin)=={'reported_individual_parameter_mean'}
    assert set(n[n.species=='mouse'].endpoint_origin)=={'reported_composite_profile'}
    # Only the two rat doses share the observed eight-hour AUC window.
    m=x[x.study_id=='CEBS:K02922']
    assert set(m.dose_mg_kg)=={10,50} and set(m.auc_definition)=={'AUC0-8h'}
    assert m.cmax.isna().all()
    # Source mismatch at 400 mg/kg in rat protocol must not enter the subset.
    f=x[(x.study_id=='CEBS:K07076') & (x.species=='rat')]
    assert 400 not in set(f.dose_mg_kg)
    assert not f.high_confidence_auc.any()
    # Actual dose, parent analyte only; the same workbook includes a metabolite.
    o=x[x.study_id=='CEBS:S0593']
    assert 107.88 in set(o.dose_mg_kg)
    assert all(json.loads(v)['Is_Parent_Compound']=='T' for v in o.original_json)

def test_reject_censored_or_ambiguous_numeric_cells():
    assert numeric('1,570')==1570
    assert numeric('') is None
    for value in ['1,57','<0.5','ND']:
        with pytest.raises(ValueError):numeric(value)
