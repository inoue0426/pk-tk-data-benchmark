import json
from pathlib import Path
import pandas as pd
import pytest
from pktk.collectors.http import Collector
from pktk.collectors.pmc_tables import collect,mean_sd
from pktk.standardization.resolve import valid_cas,choose,detail_candidate
from pktk.expansion import build_series


def test_cas_checksum_and_ambiguous_identities():
    assert valid_cas('477600-75-2')
    assert not valid_cas('477600-75-3')
    assert not valid_cas('tofacitinib')
    a={'eligible':True,'inchikey':'one'};b={'eligible':True,'inchikey':'two'}
    assert choose([a,b])[0] is None
    assert choose([a],complete=False)[0] is None
    assert choose([a])[0]==a


def test_conflicting_cas_or_name_not_forced():
    d={'chebi_accession':'CHEBI:test','name':'ethanol','default_structure':{
        'smiles':'CCO','standard_inchi_key':'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'},
       'database_accessions':{'CAS':[{'type':'CAS','accession_number':'64-17-5'}]}}
    assert detail_candidate(d,'64-17-5',['ethanol'])['eligible']
    assert not detail_candidate(d,'67-56-1',['ethanol'])['eligible']
    assert not detail_candidate(d,'64-17-5',['methanol'])['eligible']


def test_mean_sd_refuses_ranges_and_censoring():
    assert mean_sd('13,000 ± 1910')==(13000,1910)
    assert mean_sd('1840 ± 349 b')==(1840,349)
    with pytest.raises(ValueError):mean_sd('<3')
    with pytest.raises(ValueError):mean_sd('3 (1-5)')


@pytest.mark.skipif(not Path('data/raw/PMC6681021.xml').exists(),reason='Requires public XML cache')
def test_actual_tables_and_units():
    specs=json.loads(Path('configs/curation/pmc_expansion.json').read_text())
    x,_,failures=collect(Collector('.',offline=True),specs)
    assert not failures
    assert len(x)==41
    tof=x[(x.compound_name=='tofacitinib') & (x.route=='intravenous')]
    assert tof.auc.iloc[0]==pytest.approx(141/60)
    da=x[(x.compound_name=='DA-6886') & (x.route=='oral')]
    assert da.cmax.tolist()==pytest.approx([.017,.393,1.44])
    herb=x[(x.publication_doi=='10.3390/pharmaceutics10030133') & x.route.eq('oral')]
    assert len(herb)==8 and not herb.primary_eligible.any()
    # No dose-normalized AUC rows have been mistaken for measured AUC.
    assert tof.auc.iloc[1]==pytest.approx(342/60)


@pytest.mark.skipif(not Path('data/raw/chebi_71200.json').exists(),reason='Requires public identity cache')
def test_invalid_smiles_reconstructed_only_from_source_inchi():
    d=json.loads(Path('data/raw/chebi_71200.json').read_text())
    r=detail_candidate(d,'477600-75-2',['tofacitinib'])
    assert r['eligible']
    assert r['inchikey']=='UJLAWZDWDVHWOW-YPMHNXCESA-N'
    assert r['structure_input_kind']=='source_inchi_fallback_invalid_smiles'


def test_post_cache_request_body_is_part_of_identity(tmp_path):
    import hashlib
    c=Collector(tmp_path,offline=True)
    p=c.raw/'x.json';p.write_text('{}')
    meta=dict(url='https://example.org/query',version='v',request_body={'compound':'a'},sha256=hashlib.sha256(b'{}').hexdigest())
    p.with_suffix('.json.json').write_text(json.dumps(meta))
    assert c.fetch('x.json',meta['url'],'v',json_body={'compound':'a'})[0]==p
    with pytest.raises(ValueError):c.fetch('x.json',meta['url'],'v',json_body={'compound':'b'})


@pytest.mark.skipif(not Path('data/processed/reported_endpoints_expanded.csv').exists(),reason='Run expanded pipeline')
def test_expanded_grouping_no_publication_or_endpoint_leakage():
    x=pd.read_csv('data/processed/reported_endpoints_expanded.csv',keep_default_na=False)
    x=x[x.compound_name.isin(['HY-071085','tofacitinib','DA-6886','lumefantrine','saikosaponin A','imperatorin','paeonol','memantine'])].copy()
    for c in ['auc','cmax']:x[c]=pd.to_numeric(x[c],errors='coerce')
    s=build_series(x)
    assert len(s)==17
    assert s.primary_eligible.sum()==9
    assert s[s.primary_eligible].molecular_group_id.nunique()==7
    assert s[s.primary_eligible].publication_doi.nunique()==5
    # Same chemistry in a different publication must create a separate series.
    g=x[(x.protocol_id==x.protocol_id.iloc[0]) & (x.sex==x.sex.iloc[0])].copy()
    g['study_id']='doi:other';g['publication_doi']='other'
    g['source_record_id']='other:'+g.source_record_id
    assert len(build_series(pd.concat([x,g],ignore_index=True)))==len(s)+1
    assert not s[s.compound_name.eq('memantine')].primary_eligible.any()


def test_registry_merges_source_identifiers_without_losing_conflicts():
    from pktk.standardization.resolve import registry
    obs=pd.DataFrame([dict(compound_id='DTX:one',species='rat',parent_analyte=True,
        compound_name='ethanol',casrn='64-17-5',original_json='{}')])
    endpoints=pd.DataFrame([dict(compound_id='DTX:one',compound_name='ethyl alcohol',casrn='67-56-1')])
    r=registry(obs,endpoints)
    assert len(r)==1
    assert r[0]['cas_values']==['64-17-5','67-56-1']
    assert set(r[0]['names'])=={'ethanol','ethyl alcohol'}


@pytest.mark.skipif(not Path('data/raw/TKS_Formamide_1998_S0613_Group_TK_Parameters_2_copy.xlsx').exists(),reason='Requires CEBS cache')
def test_cebs_endpoints_are_not_predicted_peaks_or_nominal_doses():
    from pktk.collectors.cebs import collect as collect_cebs
    from pktk.qc.composite import audit
    x,_=collect_cebs(Collector('.',offline=True))
    assert len(x)==48
    form=x[x.compound_name=='Formamide']
    assert len(form)==24 and form.cmax.isna().all()
    assert form.high_confidence_auc.all()
    b=x[(x.compound_name=='Benzophenone') & (x.species=='rat') & (x.sex=='male')]
    assert b.dose_mg_kg.tolist()==pytest.approx([1.88,3.84,7.78])
    assert b.auc.iloc[0]==pytest.approx(32.7/60)
    iso=x[x.compound_name=='Isoeugenol']
    assert not iso[iso.species=='rat'].primary_eligible.any()
    assert iso[iso.species=='mouse'].high_confidence_cmax.all()
    assert not iso.high_confidence_auc.any()
    x['molecular_group_id']=x.compound_id
    s=build_series(x)
    assert len(s)==16 and s.primary_eligible.sum()==14
    a=audit(Path('.'),x)
    assert a['reviewed_candidates']==18
    assert a['status_counts']['validated_group_endpoints']==14


@pytest.mark.skipif(not Path('data/raw/PMC11506582.xml').exists(),reason='Requires round-two XML cache')
def test_round_two_cells_days_and_license():
    specs=json.loads(Path('configs/curation/pmc_round2.json').read_text())
    x,_,f=collect(Collector('.',offline=True),specs)
    assert not f
    me=x[x.compound_name=='methyl eugenol']
    assert len(me)==6 and me.protocol_id.nunique()==2
    assert me[me.regimen=='once'].cmax.iloc[0]==pytest.approx(.031258)
    assert me[me.regimen!='once'].cmax.iloc[0]==pytest.approx(.050602)
    assert not me[me.regimen!='once'].high_confidence_auc.any()
    lor=x[x.compound_name=='lorlatinib']
    assert lor.auc_definition.eq('AUC0-24h').all()
    assert lor.license_url.str.contains('by-nc').all()
    # An unreviewed restricted license cannot pass the default open-license adapter.
    bad=dict(specs[0]);bad.pop('license_url')
    _,_,f=collect(Collector('.',offline=True),[bad]);assert f
