"""Cache discovered annotated CEBS group workbooks; discovery does not imply validation."""
import argparse,json
from pathlib import Path
from urllib.parse import unquote
from pktk.collectors.http import Collector
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--offline',action='store_true');args=p.parse_args()
c=Collector(root,args.offline)
items=json.loads((root/'reports/composite_source_discovery.json').read_text());out=[]
for item in items:
 for url in item['files']:
  if 'group' not in url.lower() or not url.lower().endswith('.xlsx'):continue
  name=unquote(url.split('/')[-1])
  try:
   version=json.loads((c.raw/(name+'.json')).read_text())['version'] if (c.raw/(name+'.json')).exists() else 'CEBS-expansion-2026-10-28'
   path,meta=c.fetch(name,url,version)
   out.append(dict(document_id=item['document_id'],file=name,manifest=meta));print(item['document_id'],name,flush=True)
  except Exception as e:out.append(dict(document_id=item['document_id'],file=name,error=str(e)))
(root/'reports/composite_workbooks.json').write_text(json.dumps(out,indent=2)+'\n')
