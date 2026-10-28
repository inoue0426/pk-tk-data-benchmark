"""Supplement names only when EPA httk CAS and DTXSID agree exactly."""
import json,tarfile
import rdata


def enrich(collector,entries):
    path,meta=collector.fetch('httk_2.7.4.tar.gz','https://cran.r-project.org/src/contrib/httk_2.7.4.tar.gz',
        'httk-2.7.4','e186171a90d97fe456e7836ec3096719feb49354afd71860df2d2a5c8c4c51b2')
    # Extract one known data member without unpacking arbitrary archive paths.
    with tarfile.open(path) as archive:
        raw=archive.extractfile('httk/data/Tables.RData').read()
    temp=collector.root/'data/interim/httk_tables.RData';temp.write_bytes(raw)
    table=rdata.read_rda(temp)['chem.physical_and_invitro.data']
    events=[]
    for entry in entries:
        if len(entry['cas_values'])!=1:continue
        cas=entry['cas_values'][0]
        hit=table[(table['CAS']==cas)&(table['DTXSID']==entry['compound_id'])]
        if len(hit)!=1:continue
        names=set(str(hit.iloc[0]['All.Compound.Names']).split('|'))|{str(hit.iloc[0]['Compound'])}
        names={n for n in names if n and n!=cas and not n.upper().startswith('DTXSID')}
        added=sorted(names-set(entry['names']))
        if not added:continue
        entry['names']=sorted(set(entry['names'])|names)
        entry['alias_source_url']=meta['url'];entry['alias_raw_sha256']=meta['sha256']
        entry['alias_match_basis']='exact_CAS_and_DTXSID';entry['added_aliases_json']=json.dumps(added)
        events.append(dict(compound_id=entry['compound_id'],cas=cas,added_names=added,source_url=meta['url'],raw_sha256=meta['sha256']))
    (collector.root/'data/processed/identity_alias_evidence.json').write_text(json.dumps(events,indent=2)+'\n')
    return entries
