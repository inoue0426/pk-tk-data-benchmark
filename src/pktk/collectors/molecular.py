"""Bounded, cached public molecular-PK source collection, separate from Dataset B."""
import json
from datetime import datetime,timezone
from urllib.parse import urljoin, urlencode
from pktk.collectors.http import Collector, sha256

BASE='https://www.ebi.ac.uk/chembl/api/data/'

def collect(root,offline=False):
    cfg=json.loads((root/'configs/molecular/sources.json').read_text())
    c=Collector(root,offline); manifest=[];live_checked=False
    def get(name,url,version):
        nonlocal live_checked
        if not offline and not live_checked and not (c.raw/name).exists() and name!='molecular_chembl_status.json':
            day=datetime.now(timezone.utc).strftime('%Y-%m-%d')
            check,_=c.fetch('molecular_release_check_'+day+'.json',BASE+'status.json','release-check-'+day)
            if json.loads(check.read_text())['chembl_db_version']!=cfg['chembl_release']:
                raise ValueError('Live ChEMBL release differs; never mix new pages with an old snapshot')
            live_checked=True
        p,m=c.fetch(name,url,version)
        manifest.append(dict(file=name,**m))
        return json.loads(p.read_text())
    status=get('molecular_chembl_status.json',BASE+'status.json',cfg['chembl_release'])
    if status['chembl_db_version']!=cfg['chembl_release']:raise ValueError('ChEMBL release changed; review and version a new snapshot')
    url=cfg['chembl_url'];rows=[];seen=set();page=0;total=None
    while url and len(rows)<cfg['maximum_activity_records']:
        if url in seen:raise ValueError('Pagination loop')
        seen.add(url)
        data=get(f'molecular_chembl_activity_{page:03d}.json',url,cfg['chembl_release'])
        total=data['page_meta']['total_count'];rows.extend(data['activities'])
        next_url=data['page_meta']['next'];url=urljoin('https://www.ebi.ac.uk',next_url) if next_url else None
        if url and not url.startswith(BASE):raise ValueError('Unexpected pagination origin')
        print('ChEMBL activity rows',len(rows),flush=True);page+=1
    rows=rows[:cfg['maximum_activity_records']]
    if len({r['activity_id'] for r in rows})!=len(rows):raise ValueError('Duplicate activity IDs across pages')
    # Original assay and publication metadata, in bounded API batches.
    metadata={}
    for resource,field,plural in [('assay','assay_chembl_id','assays'),('document','document_chembl_id','documents')]:
        ids=sorted({r[field] for r in rows if r.get(field)})
        records=[]
        for offset in range(0,len(ids),100):
            query=urlencode({field+'__in':','.join(ids[offset:offset+100]),'limit':100})
            data=get(f'molecular_chembl_{resource}_{offset//100:03d}.json',BASE+resource+'.json?'+query,cfg['chembl_release'])
            if data['page_meta']['next']:raise ValueError('Metadata batch unexpectedly paginated')
            records.extend(data[plural])
            print(resource,len(records),'/',len(ids),flush=True)
        metadata[resource]={r[field]:r for r in records}
        if set(metadata[resource])!=set(ids):raise ValueError('Missing/duplicate source metadata')
    report=dict(config=cfg,query_total_records=total,collected_records=len(rows),truncated=bool(url),downloads=manifest)
    (root/'reports/molecular_collection_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    return rows,metadata,report
