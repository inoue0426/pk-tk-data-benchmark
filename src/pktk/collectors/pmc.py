"""Validated article-specific XML adapter; table metadata changes fail closed."""
import json
import re
import xml.etree.ElementTree as ET
import pandas as pd

URL='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9316000/fullTextXML'
PREFIX='pharmaceuticals-15-00801-'

def text(element):
    return ''.join(element.itertext()).strip()

def mean_sd(cell):
    m=re.fullmatch(r'\s*([0-9.]+)\s*±\s*([0-9.]+)\s*\**\s*',cell)
    if not m: raise ValueError('Unrecognized mean/SD: '+cell)
    return float(m[1]),float(m[2])

def parse(path, meta):
    root=ET.parse(path).getroot()
    if 'creativecommons.org/licenses/by/4.0' not in ET.tostring(root.find('.//permissions'),encoding='unicode'):
        raise ValueError('Article license changed')
    full=text(root)
    for evidence in ['Sprague Dawley','12 h','0.5% sodium carboxymethylcellulose']:
        if evidence not in full: raise ValueError('Missing protocol evidence: '+evidence)
    table=root.find(f'.//table-wrap[@id="{PREFIX}t001"]')
    headings=[text(x) for x in table.findall('./table/thead/tr')[1]]
    if headings[:3]!=['S-IG (3.25 mg/kg)','S-IG (6.5 mg/kg)','S-IG (13 mg/kg)']:
        raise ValueError('Dose headers changed')
    sexes=[text(x).lower() for x in table.findall('./table/thead/tr')[2]][:6]
    if sexes!=['male','female']*3: raise ValueError('Sex headers changed')
    values={text(tr[0]):[text(c) for c in list(tr)[1:]] for tr in table.findall('./table/tbody/tr')}
    design=root.find(f'.//table-wrap[@id="{PREFIX}t005"]')
    design_rows=design.findall('./table/tbody/tr')[:6]
    sampling=text(design_rows[0][-1])
    if sampling!='Blood: 5, 10, 15, 30, 45, 60, 90, 120, 360, 480 min':
        raise ValueError('Sampling schedule changed')
    rows=[]
    for i,(sex,dose) in enumerate(zip(sexes,[3.25,3.25,6.5,6.5,13.,13.])):
        d=[text(c) for c in design_rows[i]]
        if d[1]!=('M' if sex=='male' else 'F') or d[3]!='HY-071085' or d[4]!='S-IG' or float(d[6])!=dose:
            raise ValueError('Design/table mismatch')
        auc,auc_sd=mean_sd(values['AUC0-tn (h·ng/mL)'][i])
        cmax,cmax_sd=mean_sd(values['Cmax (ng/mL)'][i])
        rid=f'PMC9316000:Table1:column{i+1}'
        rows.append(dict(source='pmc',source_record_id=rid,compound_id='publication:HY-071085',compound_name='HY-071085',
            canonical_smiles='',inchikey='',species='rat',strain='Sprague Dawley',sex=sex,route='oral',
            matrix='plasma',dose_mg_kg=dose,dose_unit='mg/kg',auc=auc/1000,auc_sd=auc_sd/1000,
            auc_unit='ug*h/mL',auc_definition='AUC0-tn_reported',auc_terminal_time_h=None,
            cmax=cmax/1000,cmax_sd=cmax_sd/1000,cmax_unit='ug/mL',endpoint_origin='reported_group_mean',
            study_id='doi:10.3390/ph15070801',document_id='doi:10.3390/ph15070801',
            experiment_id=rid,protocol_id='PMC9316000:rat:single_oral:Table5',n_animals=int(d[5]),
            sampling_schedule_h=json.dumps([v/60 for v in [5,10,15,30,45,60,90,120,360,480]]),
            regimen='once',vehicle='0.5% sodium carboxymethylcellulose',fasting='12 h',age='6 weeks',
            publication_doi='10.3390/ph15070801',publication_url='https://pmc.ncbi.nlm.nih.gov/articles/PMC9316000/',
            source_url=meta['url'],retrieved_at=meta['retrieved_at'],dataset_version=meta['version'],raw_sha256=meta['sha256'],
            qc_flags='missing_structure;auc_terminal_time_unverified;group_mean_only',
            high_confidence_cmax=True,high_confidence_auc=False,has_numeric_curve=False,
            original_json=json.dumps({'dose':dose,'dose_unit':'mg/kg','cmax':values['Cmax (ng/mL)'][i],
                'cmax_unit':'ng/mL','auc':values['AUC0-tn (h·ng/mL)'][i],'auc_unit':'h*ng/mL',
                'sampling':sampling,'design_row':d,'table_id':PREFIX+'t001','column':i+1})))
    return pd.DataFrame(rows)

def collect(collector,checksum=None):
    path,meta=collector.fetch('PMC9316000.xml',URL,'PMC9316000-table1-2022',checksum)
    return parse(path,meta),meta
