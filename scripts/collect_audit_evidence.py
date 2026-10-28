"""Cache protocol evidence without interpreting it as validated endpoint data."""
import argparse,json
from pathlib import Path
from urllib.parse import urlencode,unquote
from pktk.collectors.http import Collector
p=argparse.ArgumentParser();p.add_argument('--offline',action='store_true');args=p.parse_args()
root=Path(__file__).resolve().parents[1];c=Collector(root,args.offline);out=[]
for item in json.loads((root/'reports/composite_source_discovery.json').read_text()):
 urls=[]
 if item.get('pmid') and item['document_id'] not in ['CvTdb:document:25297']:
  urls.append(('audit_pubmed_'+item['pmid']+'.json','https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urlencode({'query':'EXT_ID:'+item['pmid']+' AND SRC:MED','format':'json','resultType':'core'})))
 if item['document_id'] in ['CvTdb:document:'+v for v in ['159','160','25409','25460','25463','25474','25470']]:
  urls.extend((unquote(u.split('/')[-1]),u) for u in item['files'] if u.endswith('.xlsx') and ('Individual' in u or 'Parameters' in u))
 if item['document_id']=='CvTdb:document:131':urls.append(('cebs_678-39-7.html','https://cebs.niehs.nih.gov/cebs/test_article/678-39-7'))
 for name,url in urls:
  try:
   path,meta=c.fetch(name,url,'composite-evidence-2026-10-28')
   out.append(dict(document_id=item['document_id'],file=name,manifest=meta));print(item['document_id'],name,flush=True)
  except Exception as e:out.append(dict(document_id=item['document_id'],file=name,url=url,error=str(e)))
(root/'reports/composite_evidence_manifest.json').write_text(json.dumps(out,indent=2)+'\n')
