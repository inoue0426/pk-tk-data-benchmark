"""Reviewed source-cell mappings for the second expansion; measurements stay in source XML."""
import json,hashlib,xml.etree.ElementTree as E
from pathlib import Path

def base(pmc,table,compound):
    p=Path('data/raw',pmc+'.xml');root=E.parse(p).getroot();t=root.find(f'.//table-wrap[@id="{table}"]')
    rows=[[' '.join(''.join(c.itertext()).split()) for c in r] for r in t.findall('.//tr')]
    return dict(pmc=pmc,table_id=table,compound_name=compound,source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
        table_sha256=hashlib.sha256(json.dumps(rows,ensure_ascii=False).encode()).hexdigest(),species='rat',strain='Sprague-Dawley',
        sex='male',route='oral',matrix='plasma',regimen='once',endpoint_origin='reported_group_mean',
        auc_unit='ug*h/mL',cmax_unit='ug/mL',infusion_duration_min=0,protocol_verified=True,primary_eligible=True,
        auc_factor=.001,cmax_factor=.001,n_animals=[6,6,6],license_url='https://creativecommons.org/licenses/by/4.0/')

x=base('PMC12570999','t0004','lorlatinib')
x.update(doses=[5,10,50],columns=[1,2,3],auc_row='AUC(0–24h) (ng·h/mL)',cmax_row='Cmax (ng/mL)',auc_definition='AUC0-24h',
    vehicle='0.5% carboxymethylcellulose sodium',fasting='12 h',license_url='https://creativecommons.org/licenses/by-nc/4.0/',
    evidence=['Rats were fasted for 12 h','All drug solutions were made up in 0.5% carboxymethylcellulose sodium','low-dose group, single dose of lorlatinib (5 mg/kg)','high-dose group, single dose of lorlatinib (50 mg/kg)'],
    flags=['license_noncommercial','source_total_animal_count_inconsistent','AUCinf_large_tail_not_selected'])
specs=[x]
for day in [1,28]:
 x=base('PMC11506582','cimb-46-00673-t001','methyl eugenol')
 start=1 if day==1 else 4
 cells={ep:[[start,col],[start+1,col-1],[start+2,col-1]] for ep,col in [('auc',5 if day==1 else 6),('cmax',4)]}
 x.update(species='mouse',strain='Kunming',compound_id='DTXSID2061045',casrn='93-15-2',doses=[20,40,60],columns=list(range(start,start+3)),
    auc_row=None,cmax_row=None,endpoint_cells=cells,auc_definition='AUC0-inf' if day==1 else 'AUC0-t_unspecified',auc_eligible=day==1,
    protocol_suffix='day'+str(day),regimen='once' if day==1 else 'daily; day 28',vehicle='0.5% CMC-Na',fasting='not reported',
    evidence=['Twenty-four male Kunming (KM) mice','ME dissolved in 0.5% CMC-Na','on Days 1 and 28'],
    flags=['fasting_not_reported','shared_animals_across_days']+([] if day==1 else ['auc_terminal_time_unresolved','repeat_dose_AUCinf_not_selected']))
 specs.append(x)
# DTXSID comes from the existing CvTdb CAS registry, never guessed from name.
import pandas as pd
obs=pd.read_csv('data/processed/observations.csv',keep_default_na=False,low_memory=False)
ids=obs.loc[obs.casrn.eq('93-15-2') & obs.parent_analyte,'compound_id'].unique()
assert len(ids)==1
for x in specs[1:]:x['compound_id']=ids[0]
Path('configs/curation/pmc_round2.json').write_text(json.dumps(specs,indent=2,ensure_ascii=False)+'\n')
