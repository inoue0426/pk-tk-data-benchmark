"""Conservative ChEBI resolution with explicit alternatives and exact InChIKey checks."""
import hashlib,html,json,re
from urllib.parse import urlencode
from pathlib import Path
import pandas as pd
from rdkit import Chem
from pktk.standardization.core import clean,structure

BASE='https://www.ebi.ac.uk/chebi/backend/api/public/'

def norm(value):
    return ' '.join(html.unescape(re.sub('<[^>]+>','',str(value))).casefold().split())

def valid_cas(value):
    if not re.fullmatch(r'\d{2,7}-\d{2}-\d',value):return False
    digits=value.replace('-','')
    return sum((i+1)*int(n) for i,n in enumerate(reversed(digits[:-1])))%10==int(digits[-1])

def detail_candidate(detail,cas,names):
    synonyms={norm(detail.get('name','')),norm(detail.get('ascii_name',''))}
    for group in detail.get('names',{}).values():
        for r in group:
            synonyms.update(norm(r.get(k,'')) for k in ['name','ascii_name'])
    xrefs=[r for group in detail.get('database_accessions',{}).values() for r in group]
    cas_values=sorted({r['accession_number'] for r in xrefs if r['type']=='CAS' or r.get('source_name')=='CAS'})
    cids=sorted({r['accession_number'] for r in xrefs if r.get('prefix')=='pubchem.compound' or r.get('source_name')=='PubChem'})
    st=detail.get('default_structure') or {}
    standardized=structure(st.get('smiles',''))
    input_kind='source_smiles'
    if standardized['structure_flag']=='invalid_structure' and st.get('standard_inchi'):
        mol=Chem.MolFromInchi(st['standard_inchi'])
        if mol is not None:
            standardized=structure(Chem.MolToSmiles(mol,isomericSmiles=True))
            input_kind='source_inchi_fallback_invalid_smiles'
    key=st.get('standard_inchi_key','')
    name_match=bool({norm(n) for n in names if n}&synonyms)
    cas_match=bool(cas and cas in cas_values)
    eligible=bool(standardized['inchikey'] and standardized['inchikey']==key and name_match and (not cas or cas_match))
    return dict(chebi_id=detail['chebi_accession'],name=detail['name'],cas_values=cas_values,
                pubchem_cids=cids,source_inchikey=key,cas_match=cas_match,name_match=name_match,
                eligible=eligible,structure_input_kind=input_kind,source_smiles=st.get('smiles',''),**standardized)


def choose(candidates,complete=True):
    eligible=[c for c in candidates if c['eligible']]
    keys={c['inchikey'] for c in eligible}
    if not complete:return None,'search_incomplete'
    if len(keys)>1:return None,'ambiguous_multiple_structures'
    if len(keys)==1:return eligible[0],'resolved_exact'
    return None,'unresolved_or_identifier_conflict'


def registry(observations,endpoints):
    rows=[]
    pre=observations[(observations.species!='human') & observations.parent_analyte]
    for cid,g in pre.groupby('compound_id'):
        originals=[json.loads(v) for v in g.original_json.drop_duplicates()]
        names=set(g.compound_name)
        for r in originals:
            names.update(r.get(k,'') for k in ['analyzed_chem_name','analyzed_chem_name_original','analyzed_chem_name_secondary_original'])
        rows.append(dict(compound_id=cid,names=sorted(n for n in names if n),cas_values=sorted(set(g.casrn)-{''}),source='cvtdb'))
    for cid,g in endpoints.groupby('compound_id'):
        rows.append(dict(compound_id=cid,names=sorted(set(g.compound_name)),cas_values=sorted(set(g.get('casrn',pd.Series(dtype=str)))-{''}),source='pmc'))
    # Multiple sources may already share a curated identifier (e.g. CEBS/DTXSID).
    # Merge evidence, never manufacture another registry row for the same identifier.
    merged={}
    for r in rows:
        if r['compound_id'] not in merged:merged[r['compound_id']]=r.copy()
        else:
            old=merged[r['compound_id']]
            for k in ['names','cas_values']:old[k]=sorted(set(old[k])|set(r[k]))
    return list(merged.values())


def resolve(collector,entries,priority_ids=(),limit=None,bulk_candidates=None):
    bulk_candidates=bulk_candidates or {}
    entries=sorted(entries,key=lambda r:(r['source']!='pmc',r['compound_id'] not in priority_ids,r['compound_id']))
    output=[];candidates_out=[];manifests=[]
    for i,entry in enumerate(entries):
        record={**entry,'names_json':json.dumps(entry['names']),'cas_values_json':json.dumps(entry['cas_values']),
                'identity_status':'not_queried','canonical_smiles':'','inchikey':'','pubchem_cid':'','pubchem_cids_json':'[]',
                'chebi_id':'','resolved_cas_json':'[]','scaffold':'','structure_flag':'missing_structure','molecular_weight':None}
        del record['names'];del record['cas_values']
        if limit is not None and i>=limit:output.append(record);continue
        if len(entry['cas_values'])>1:
            record['identity_status']='ambiguous_source_cas';output.append(record);continue
        cas=entry['cas_values'][0] if entry['cas_values'] else ''
        if cas and not valid_cas(cas):
            record['identity_status']='invalid_cas_checksum';output.append(record);continue
        query=cas or (entry['names'][0] if entry['names'] else '')
        if not query:record['identity_status']='missing_identifiers';output.append(record);continue
        record['identity_query']=query
        candidates=list(bulk_candidates.get(entry["compound_id"],[]));complete=True
        try:
            url=BASE+'es_search/?'+urlencode(dict(term=query,page=1,size=20))
            name='chebi_search_'+hashlib.sha256(query.encode()).hexdigest()[:20]+'.json'
            # Consult cached ChEBI evidence even when EPA resolves the identifier,
            # so source conflicts cannot be hidden by changing provider priority.
            use_chebi=not any(c['eligible'] for c in candidates) or (collector.raw/name).exists()
            hits=[];result={'number_pages':1,'total':0}
            if use_chebi:
                p,m=collector.fetch(name,url,'identity-2026-10-28',allowed_statuses=(404,));manifests.append(m)
                record['search_source_url']=m['url'];record['search_raw_sha256']=m['sha256'];record['search_retrieved_at']=m['retrieved_at']
                if m.get('http_status')==404:
                    record['search_status']='not_found'
                else:
                    result=json.loads(p.read_text());complete=result['number_pages']<=1;hits=result['results']
            else:
                record['search_status']='EPA_exact_identifiers; additional_ChEBI_query_not_needed'
            # Avoid broad fuzzy lookups: names require exact preferred name before detail retrieval.
            # CAS search hits are verified against the full record, never accepted on search score.
            if not cas:
                exact=[h for h in hits if norm(h['_source']['name']) in {norm(n) for n in entry['names']}]
                hits=exact if exact else hits[:5]
                if not exact and result['total']>5:complete=False
            if len(hits)>20:complete=False
            for hit in hits[:20]:
                chebi=hit['_id']
                p,m=collector.fetch('chebi_'+chebi+'.json',BASE+f'compound/{chebi}/','identity-2026-10-28');manifests.append(m)
                candidate=detail_candidate(json.loads(p.read_text()),cas,entry['names'])
                candidate.update(provider='ChEBI',compound_id=entry['compound_id'],query=query,source_url=m['url'],raw_sha256=m['sha256'],retrieved_at=m['retrieved_at'])
                candidates.append(candidate)
            selected,status=choose(candidates,complete)
            if not candidates and record.get('search_status')=='not_found':status='not_found'
            record['identity_qc_flags']=';'.join(sorted({flag for c in candidates for flag in c.get('identity_qc_flags',[])}))
            record['identity_status']=status
            if selected:
                for k in ['canonical_smiles','inchikey','chebi_id','scaffold','structure_flag','molecular_weight','structure_input_kind']:
                    record[k]=selected[k]
                record['identity_provider']=selected.get('provider','ChEBI')
                record['identity_source_member']=selected.get('source_member','')
                record['identity_source_row']=selected.get('source_row','')
                record['resolved_cas_json']=json.dumps(selected['cas_values'])
                record['identity_source_url']=selected['source_url'];record['identity_raw_sha256']=selected['raw_sha256']
                cids=selected['pubchem_cids']
                # UniChem exact structure lookup is a cross-reference, not a parent/salt transformation.
                # Resolve optional CIDs for primary/candidate compounds and reuse every
                # existing cross-reference cache. Chemical identity does not require a CID.
                lookup_cid=entry['source']=='pmc' or entry['compound_id'] in priority_ids or (collector.raw/('unichem_'+record['inchikey']+'.json')).exists()
                if lookup_cid:
                    try:
                        url='https://www.ebi.ac.uk/unichem/api/v1/compounds'
                        p,m=collector.fetch('unichem_'+record['inchikey']+'.json',url,'identity-2026-10-28',json_body={'compound':record['inchikey'],'type':'inchikey'});manifests.append(m)
                        links=json.loads(p.read_text())
                        cids=sorted(set(cids)|{str(r['compoundId']) for compound in links.get('compounds',[]) if Chem.InchiToInchiKey(compound['inchi']['inchi'])==record['inchikey'] for r in compound['sources'] if r['id']==22})
                    except Exception as e:
                        if collector.offline and isinstance(e,FileNotFoundError):raise
                        record['pubchem_lookup_status']='unavailable: '+str(e)[:160]
                elif not cids:
                    record['pubchem_lookup_status']='not_queried_optional_crossref'
                record['pubchem_cids_json']=json.dumps(cids)
                record['pubchem_cid']=cids[0] if len(cids)==1 else ''
                record['pubchem_lookup_status']=record.get('pubchem_lookup_status','unique_crossref' if len(cids)==1 else 'ambiguous_or_absent_crossref')
        except Exception as e:
            if collector.offline and isinstance(e,FileNotFoundError):raise
            record['identity_status']='lookup_failed';record['error']=str(e)
        candidates_out.extend(candidates);output.append(record)
        print(f'identity {i+1}/{len(entries)} {entry["compound_id"]}: {record["identity_status"]}',flush=True)
    registry=pd.DataFrame(output)
    registry.to_csv(collector.root/'data/processed/compound_identity_registry.csv',index=False)
    (collector.root/'data/processed/compound_identity_candidates.json').write_text(json.dumps(candidates_out,indent=2,allow_nan=False)+'\n')
    # POST lookups share a URL; request bodies are part of provenance identity.
    manifest_keys={json.dumps([m['url'],m.get('request_body')],sort_keys=True):m for m in manifests}
    (collector.root/'reports/identity_source_manifest.json').write_text(json.dumps(list(manifest_keys.values()),indent=2)+'\n')
    return registry
