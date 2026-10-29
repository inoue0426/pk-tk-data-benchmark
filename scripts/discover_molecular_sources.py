"""Replay the negative-source and original-publication availability audit."""
import argparse,json,zipfile
from pathlib import Path
from urllib.parse import urlencode
import pandas as pd
from pktk.collectors.http import Collector

def main():
    p=argparse.ArgumentParser();p.add_argument('--offline',action='store_true');a=p.parse_args()
    root=Path(__file__).resolve().parents[1];c=Collector(root,a.offline)
    archive,m=c.fetch('frdb-v2024-12-30.zip','https://drugs.ncats.io/downloads-public/frdb-v2024-12-30.zip','FRDB-2024-12-30')
    with zipfile.ZipFile(archive) as z:x=pd.read_csv(z.open('frdb/frdb-pk.tsv'),sep='\t',keep_default_na=False)
    pub,meta=c.fetch('pksmart_supp_metadata.json','https://api.figshare.com/v2/articles/30224986','discovery-2026-10-29')
    pk=json.loads(pub.read_text());records=[]
    candidates=pd.read_csv(root/'data/processed/molecular/dose_candidates_unverified.csv')
    for i,doi in enumerate(sorted(candidates.publication_id.unique())):
        query='DOI:"'+doi.removeprefix('doi:')+'"'
        url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urlencode({'query':query,'format':'json','resultType':'core'})
        path,provenance=c.fetch(f'molecular_candidate_lookup_{i}.json',url,'dose-candidate-discovery-2026-10-29')
        data=json.loads(path.read_text())
        found=[dict(pmid=v.get('id'),pmcid=v.get('pmcid'),is_open_access=v.get('isOpenAccess'),in_pmc=v.get('inPMC'),title=v.get('title')) for v in data['resultList']['result']]
        records.append(dict(publication=doi,source_url=url,raw_sha256=provenance['sha256'],results=found,decision='Original protocol not obtained; no strict-series promotion'))
    (root/'reports/molecular_candidate_access.json').write_text(json.dumps(records,indent=2)+'\n')
    result=dict(ncats_pk_rows=len(x),ncats_species_counts=x.pk_species.value_counts().to_dict(),ncats_archive=m,
                pksmart_publisher_license=pk['license'],pksmart_supplement_metadata=meta,candidate_publications=len(records),
                original_full_protocols_obtained=0)
    (root/'reports/molecular_discovery_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['ncats_pk_rows','ncats_species_counts','candidate_publications','original_full_protocols_obtained']},indent=2))
if __name__=='__main__':main()
