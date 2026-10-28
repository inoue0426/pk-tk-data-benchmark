"""Discover original CEBS TK files for the outstanding CvTdb audit queue."""
import argparse,json,re
from pathlib import Path
from urllib.parse import urljoin
import pandas as pd
from pktk.collectors.http import Collector
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--offline',action='store_true');args=p.parse_args()
c=Collector(root,args.offline)
a=pd.read_csv(root/'reports/composite_audit.csv',keep_default_na=False)
o=pd.read_csv(root/'data/processed/observations.csv',keep_default_na=False,low_memory=False)
results=[]
for doc,g in a.groupby('document_id'):
 r=o[o.document_id==doc].iloc[0];original=json.loads(r.original_json)
 item=dict(document_id=doc,n_candidates=len(g),compound_id=r.compound_id,compound_name=r.compound_name,casrn=r.casrn,doi=original.get('doi'),pmid=original.get('pmid'),source_experiment=original.get('other_study_identifier'),files=[])
 if 'DTXSID' in str(r.compound_id) and ('ntp' in str(original.get('doi','')).lower() or str(original.get('other_study_identifier','')).startswith(('S0','C2'))):
  url=f'https://cebs.niehs.nih.gov/cebs/test_article/{r.casrn}'
  try:
   path,meta=c.fetch(f'cebs_{r.casrn}.html',url,json.loads((c.raw/f'cebs_{r.casrn}.html.json').read_text())['version'] if (c.raw/f'cebs_{r.casrn}.html.json').exists() else 'CEBS-discovery-2026-10-28')
   item['landing_manifest']=meta
   links=sorted(set(re.findall(r'''(?:href|value)=["']([^"']+)["']''',path.read_text())))
   item['files']=[urljoin(url,x) for x in links if '/get_file/' in x and ('TK' in x or 'Toxicokinetic' in x or re.search(r'(?:S|K|C)\d+.*(?:Summary|summary)',x))]
   item['status']='landing_retrieved'
  except Exception as e:item.update(status='collection_failed',error=str(e))
 else:item['status']='literature_protocol_required'
 results.append(item);print(doc,item['status'],len(item['files']),flush=True)
(root/'reports/composite_source_discovery.json').write_text(json.dumps(results,indent=2)+'\n')
