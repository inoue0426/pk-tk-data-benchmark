"""Pinned, reviewed XML table adapters. No generic inferred dose or unit mappings."""
import hashlib
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path
import pandas as pd
from pktk.standardization.core import identity


def text(node):
    return ' '.join(''.join(node.itertext()).split())


def mean_sd(cell):
    # Restrict suffixes to significance footnotes; never parse ranges or censoring as means.
    m=re.fullmatch(r'([0-9.,]+)(?:\s*±\s*([0-9.,]+))?\s*[a-z*#\s]*',cell)
    if not m: raise ValueError('Unsupported numeric endpoint cell: '+cell)
    return float(m[1].replace(',','')),float(m[2].replace(',','')) if m[2] else None


def parse(path,meta,spec):
    if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=spec['source_sha256']:
        raise ValueError('Source snapshot differs from reviewed curation')
    root=ET.parse(path).getroot()
    permissions=ET.tostring(root.find('.//permissions'),encoding='unicode')
    if not re.search(r'creativecommons.org/licenses/by/[24].0',permissions):
        raise ValueError('Unreviewed article license')
    full=text(root)
    for evidence in spec['evidence']:
        if evidence not in full: raise ValueError('Missing protocol evidence: '+evidence)
    table=root.find(f'.//table-wrap[@id="{spec["table_id"]}"]')
    rows=[[text(c) for c in tr] for tr in table.findall('.//tr')]
    if hashlib.sha256(json.dumps(rows,ensure_ascii=False).encode()).hexdigest()!=spec['table_sha256']:
        raise ValueError('Reviewed table changed')
    doi=next(e.text for e in root.findall('.//article-id') if e.get('pub-id-type')=='doi').lower()
    pubmed=next((e.text for e in root.findall('.//article-id') if e.get('pub-id-type')=='pmid'),'')
    pmc=spec['pmc'];records=[]
    for dose,column,n in zip(spec['doses'],spec['columns'],spec['n_animals'],strict=True):
        values={};original={}
        for endpoint in ['auc','cmax']:
            label=spec[endpoint+'_row']
            # Repeated matrix rows (memantine) use first explicitly selected plasma block only.
            cells=next((r for r in rows if r and r[0]==label),None) if label else None
            value,sd=mean_sd(cells[column]) if cells else (None,None)
            factor=spec[endpoint+'_factor']
            values[endpoint]=value*factor if value is not None else None
            values[endpoint+'_sd']=sd*factor if sd is not None else None
            original[endpoint]={'row':label,'column':column,'cell':cells[column] if cells else None,'factor':factor}
        rid=f'{pmc}:{spec["table_id"]}:column{column}'
        protocol=f'{pmc}:{spec["compound_name"]}:{spec["species"]}:{spec["route"]}:single'
        flags=spec['flags']+['group_summary_only','missing_structure']
        eligible=spec['primary_eligible']
        records.append(dict(source='pmc',source_record_id=rid,experiment_id=rid,compound_id='publication:'+spec['compound_name'],
            compound_name=spec['compound_name'],casrn='',pubchem_cid='',canonical_smiles='',inchikey='',
            species=spec['species'],strain=spec['strain'],sex=spec['sex'],route=spec['route'],matrix=spec['matrix'],
            dose_mg_kg=dose,dose_unit='mg/kg',dose_basis='reported; see source methods',**values,
            auc_definition=spec['auc_definition'],auc_unit=spec['auc_unit'],cmax_unit=spec['cmax_unit'],
            n_animals=n,regimen=spec['regimen'],vehicle=spec['vehicle'],fasting=spec['fasting'],
            infusion_duration_min=spec['infusion_duration_min'],protocol_id=protocol,
            study_id='doi:'+doi,document_id='doi:'+doi,publication_doi=doi,pmid=pubmed,
            publication_url=f'https://pmc.ncbi.nlm.nih.gov/articles/{pmc}/',source_url=meta['url'],
            retrieved_at=meta['retrieved_at'],dataset_version=meta['version'],raw_sha256=meta['sha256'],
            original_json=json.dumps(original,ensure_ascii=False),protocol_evidence_json=json.dumps(spec['evidence'],ensure_ascii=False),
            endpoint_origin=spec['endpoint_origin'],protocol_verified=True,primary_eligible=eligible,
            high_confidence_cmax=bool(eligible and values['cmax'] is not None),
            high_confidence_auc=bool(eligible and values['auc'] is not None),has_numeric_curve=False,
            qc_flags=';'.join(flags),curation_spec_id=identity([pmc,spec['table_id'],protocol])))
    return pd.DataFrame(records)


def collect(collector,config):
    frames=[];manifests={};failures=[]
    for spec in config:
        pmc=spec['pmc']
        try:
            path,meta=collector.fetch(pmc+'.xml',f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML',
                                      pmc+'-expansion-2026-10-28',spec['source_sha256'])
            frames.append(parse(path,meta,spec));manifests[pmc]=meta
        except Exception as e:
            failures.append({'pmc':pmc,'table':spec['table_id'],'compound':spec['compound_name'],'error':str(e)})
    return pd.concat(frames,ignore_index=True) if frames else pd.DataFrame(),manifests,failures
