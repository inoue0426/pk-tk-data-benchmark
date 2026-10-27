import hashlib
import json
from pathlib import Path
import pytest
from pktk.collectors.http import Collector
from pktk.collectors.pmc import parse,mean_sd

def test_cache_integrity_and_offline(tmp_path):
    c=Collector(tmp_path,offline=True)
    p=c.raw/'x';p.write_text('hello')
    meta={'url':'https://example.org/x','version':'1','sha256':hashlib.sha256(b'hello').hexdigest()}
    (c.raw/'x.json').write_text(json.dumps(meta))
    assert c.fetch('x',meta['url'],'1')[0]==p
    p.write_text('corrupted')
    with pytest.raises(ValueError): c.fetch('x',meta['url'],'1')
    with pytest.raises(FileNotFoundError): c.fetch('missing',meta['url'],'1')

def test_mean_sd():
    assert mean_sd('20.82 ± 3.01')==(20.82,3.01)
    with pytest.raises(ValueError): mean_sd('20.82 (3.01)')

@pytest.mark.skipif(not Path('data/raw/PMC9316000.xml').exists(),reason='Requires downloaded public XML')
def test_actual_pmc_adapter():
    meta=dict(url='public',retrieved_at='test',version='test',sha256='test')
    x=parse('data/raw/PMC9316000.xml',meta)
    assert len(x)==6
    assert x.cmax.iloc[0]==pytest.approx(.01133)
    assert x.auc.iloc[5]==pytest.approx(2.07961)
    assert x.n_animals.tolist()==[5]*6
    assert not x.high_confidence_auc.any()
    assert x.high_confidence_cmax.all()

def test_robots_denial_and_pagination_loop(tmp_path,monkeypatch):
    c=Collector(tmp_path)
    class Reply:
        status_code=200
        text='User-agent: *\nDisallow: /private\n'
        content=text.encode()
        def raise_for_status(self): pass
    monkeypatch.setattr(c.session,'get',lambda *a,**k:Reply())
    with pytest.raises(PermissionError): c.allowed('https://example.org/private/data')
    page=tmp_path/'page.json'
    page.write_text(json.dumps({'results':[1],'next':'https://example.org/api'}))
    monkeypatch.setattr(c,'fetch',lambda *a,**k:(page,{}))
    with pytest.raises(ValueError,match='loop'): list(c.pages('https://example.org/api'))
