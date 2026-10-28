"""Exact DTXSID/CAS matching against EPA's pinned public CC0 structure release."""
import json,zipfile
import pandas as pd
from pktk.standardization.core import structure

FILE='DSSTox_CCD_dump_12092025_CSVs.zip'
URL='https://clowder.edap-cluster.com/api/files/69529775e4b0731a616efc4b'
SHA='66fc9d4d3bda053ab60ee01e2ff0ed4ba030339588c40711b6939faede6cf4c3'


def candidate(row,entry,meta,source_row):
    st=structure(row['SMILES'])
    cas=entry['cas_values']
    identifier_match=row['DTXSID']==entry['compound_id']
    cas_match=len(cas)==1 and cas[0]==row['CASRN']
    flags=[]
    if not cas_match:flags.append('source_CAS_disagrees')
    if not st['inchikey']:flags.append('missing_or_invalid_source_structure')
    elif st['inchikey']!=row['INCHIKEY']:
        flags.append('nonstandard_source_InChIKey_options_unverified' if row['INCHIKEY'][-4:-3]=='N' else 'source_smiles_InChIKey_disagree')
    eligible=bool(identifier_match and cas_match and st['inchikey'] and st['inchikey']==row['INCHIKEY'])
    return dict(compound_id=entry['compound_id'],provider='EPA_DSSTox',chebi_id='',name=row['PREFERRED_NAME'],
        cas_values=[row['CASRN']],pubchem_cids=[],source_inchikey=row['INCHIKEY'],cas_match=cas_match,
        identifier_match=identifier_match,name_match=None,eligible=eligible,structure_input_kind='source_smiles',
        source_smiles=row['SMILES'],source_url=meta['url'],raw_sha256=meta['sha256'],retrieved_at=meta['retrieved_at'],
        source_member='DSSTox_CCD_dump_12092025/DSSToxCCDdump.csv',source_row=source_row,
        source_dtxcid=row['DTXCID'],query=entry['compound_id'],identity_qc_flags=flags,**st)


def collect(collector,entries):
    p,m=collector.fetch(FILE,URL,'DSSTox-Dec2025',SHA)
    lookup={e['compound_id']:e for e in entries};found={};total=0
    with zipfile.ZipFile(p) as z:
        # The other CSVs are split copies: read the full master exactly once.
        with z.open('DSSTox_CCD_dump_12092025/DSSToxCCDdump.csv') as f:
            for chunk in pd.read_csv(f,chunksize=50000,dtype=str,keep_default_na=False):
                for index,row in chunk[chunk.DTXSID.isin(lookup)].iterrows():
                    found.setdefault(row.DTXSID,[]).append(candidate(row,lookup[row.DTXSID],m,int(index)+2))
                total+=len(chunk)
    report=dict(source_manifest=m,total_release_rows=total,matched_source_identifiers=len(found),
                eligible_structure_candidates=sum(c['eligible'] for cs in found.values() for c in cs),
                acceptance='Exact DTXSID and primary CAS; source SMILES must regenerate full source InChIKey; conflicts retained; no QSAR-ready/MS-ready substitution',
                license='CC0',citation='EPA DSSTox, December 2025 CSV release; doi:10.23645/epacomptox.5588566.v8')
    (collector.root/'reports/dsstox_identity_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    return found
