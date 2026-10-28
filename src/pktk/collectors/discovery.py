"""Bounded anonymous discovery; closed/licensed records never become benchmark data."""
import io,json,zipfile
from pathlib import Path
from pypdf import PdfReader

PROBES=[
 ('pkdb_public_statistics.json','https://pk-db.com/api/v1/statistics/'),
 ('pkdb_public_studies_small.json','https://pk-db.com/api/v1/studies/?page=1&page_size=1'),
]

def collect(collector):
    events=[]
    for name,url in PROBES:
        try:
            path,meta=collector.fetch(name,url,'discovery-2026-10-28')
            data=json.loads(path.read_text())
            if name.startswith('pkdb_public_studies'):
                studies=data['data']['data']
                events.append(dict(source='PK-DB',url=url,status='public_metadata_only',
                    studies=[{k:r.get(k) for k in ['sid','name','access','licence']} for r in studies],
                    decision='Do not import closed-licence study measurements; reviewed record is human, not preclinical.',manifest=meta))
            else:
                events.append(dict(source='PK-DB',url=url,status='public_statistics',counts=data,manifest=meta))
        except Exception as e:
            events.append(dict(source='PK-DB',url=url,status='failed',error=str(e)))
    try:
        path,meta=collector.fetch('ebi_da6886_supplement.zip',
            'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9024849/supplementaryFiles','discovery-2026-10-28')
        with zipfile.ZipFile(path) as outer:
            with zipfile.ZipFile(io.BytesIO(outer.read('pharmaceutics-14-00702-s001.zip'))) as inner:
                content=inner.read('pharmaceutics-1629108-SM.pdf')
        reader=PdfReader(io.BytesIO(content))
        txt='\n'.join(p.extract_text() for p in reader.pages)
        (collector.root/'data/interim/DA6886_supplement.txt').write_text(txt)
        for phrase in ['Table S1','Table S2','Stability']:
            if phrase not in txt:raise ValueError('Supplement classification requires renewed review')
        events.append(dict(source='PMC9024849_supplement',status='downloaded_and_inspected',manifest=meta,
            pages=len(reader.pages),new_systemic_pk_groups=0,
            contents='Tissue calibration standards, MRM chromatograms, plasma stability controls; not additional in vivo dose series.'))
    except Exception as e:
        events.append(dict(source='PMC9024849_supplement',status='failed',error=str(e)))
    (collector.root/'reports/expansion_discovery.json').write_text(json.dumps(events,indent=2)+'\n')
    return events
