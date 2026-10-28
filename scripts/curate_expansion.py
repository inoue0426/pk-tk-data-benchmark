"""Maintainer utility: encode reviewed table selections without transcribing measurements."""
import hashlib,json,xml.etree.ElementTree as E
from pathlib import Path

specs=[]
def add(pmc,table,compound,doses,columns,auc_row,cmax_row,auc_factor=1,cmax_factor=1,**meta):
    path=Path('data/raw')/(pmc+'.xml');root=E.parse(path).getroot()
    t=root.find(f'.//table-wrap[@id="{table}"]')
    rows=[[' '.join(''.join(c.itertext()).split()) for c in r] for r in t.findall('.//tr')]
    # Exact table checksum and source checksum reject silent layout/content changes.
    for label in [auc_row,cmax_row]:
        if label is not None and not any(row and row[0]==label for row in rows): raise ValueError((pmc,table,label,rows))
    specs.append(dict(pmc=pmc,table_id=table,compound_name=compound,source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        table_sha256=hashlib.sha256(json.dumps(rows,ensure_ascii=False).encode()).hexdigest(),doses=doses,columns=columns,
        auc_row=auc_row,cmax_row=cmax_row,auc_factor=auc_factor,cmax_factor=cmax_factor,
        species='rat',strain='Sprague-Dawley',sex='male',matrix='plasma',regimen='once',
        endpoint_origin='reported_group_mean',auc_definition='AUC0-inf',auc_unit='ug*h/mL',cmax_unit='ug/mL',
        protocol_verified=True,primary_eligible=True,**meta))

add('PMC6681021','pharmaceutics-11-00318-t001','tofacitinib',[5,10,20,50],[1,2,3,4],
    'AUC (μg∙min/mL) a',None,1/60,route='intravenous',n_animals=[9,8,7,7],
    vehicle='0.9% NaCl with 0.5% beta-cyclodextrin',fasting='food available after recovery',infusion_duration_min=1,
    evidence=['Male Sprague-Dawley rats','time zero to infinity','0.5% β-cyclodextrin','injected via the jugular vein for 1 min'],flags=[])
add('PMC6681021','pharmaceutics-11-00318-t002','tofacitinib',[10,20,50,100],[1,2,3,4],
    'AUC (μg∙min/mL) a','Cmax (μg/mL)',1/60,route='oral',n_animals=[7,8,9,7],
    vehicle='not explicitly specified in oral subsection',fasting='overnight',infusion_duration_min=0,
    evidence=['For oral administration, approximately 1.0 mL tofacitinib','rats were fasted overnight','time zero to infinity'],flags=['oral_vehicle_not_explicit'])
for table,route,cmax in [('t002','intravenous',None),('t003','oral','Cmax (ng/mL) 1, 2')]:
 add('PMC9024849','pharmaceutics-14-00702-'+table,'DA-6886',[2,10,20],[1,2,3],
     'AUC0–∞ (ng∙h/mL) 1, 2',cmax,.001,.001,route=route,n_animals=[4,4,4],vehicle='water for injection',fasting='18 h',
     infusion_duration_min=1 if route=='intravenous' else 0,
     evidence=['Male Sprague-Dawley rats','animals were fasted for 18 h','mg as free form per kilogram','dissolved in water for injection'],flags=[])
add('PMC3305926','T2','lumefantrine',[10,20,40],[2,3,4],'AUC0-∞ (hr*ng/mL)','Cmax (ng/mL)',.001,.001,
    route='oral',n_animals=[5,5,5],vehicle='0.25% CMC suspension',fasting='12-14 h',infusion_duration_min=0,
    evidence=['Male SD rats weighing 200-220 g were fasted overnight','groups of five rats at each dose level','0.25% CMC suspension'],flags=[])
for compound,doses,cols,animals,vehicle in [
 ('saikosaponin A',[2.3,4.6],[1,2],[5,4],'DMSO:Tween80:water 2:1:7'),
 ('paeonol',[14.8,29.6,59.2],[3,4,5],[5,5,6],'DMSO:Tween80:20% Solutol in water 2:2:6'),
 ('imperatorin',[2.9,5.8,11.5],[6,7,8],[5,6,5],'DMSO:Tween80:water 2:2:6')]:
 add('PMC6160927','pharmaceutics-10-00133-t001',compound,doses,cols,'AUC0–inf (μg∙min/mL) a',None,1/60,
     route='intravenous',n_animals=animals,vehicle=vehicle,fasting='not reported',infusion_duration_min=1,
     evidence=['Male Sprague-Dawley rats','each pure active component','infused for 1 min via the jugular vein'],flags=['fasting_not_reported'])
# Confounded oral comparisons are retained as explicitly ineligible, not silently removed.
for table,compound,doses,cols,animals,auclabel in [
 ('t002','saikosaponin A',[2.3,4.6],[2,4],[4,5],'AUC0–8h (μg∙min/mL) a'),
 ('t003','paeonol',[14.8,29.6,59.2],[1,3,5],[4,5,4],'AUC0–inf (μg∙min/mL) a'),
 ('t004','imperatorin',[1.4,2.9,5.8],[1,3,5],[6,6,5],'AUC0–inf (μg∙min/mL) a')]:
 add('PMC6160927','pharmaceutics-10-00133-'+table,compound,doses,cols,auclabel,'Cmax (μg/mL) a',1/60,
     route='oral',n_animals=animals,vehicle='ethanol proportion changes with dose',fasting='not reported',infusion_duration_min=0,
     evidence=['ethanol with dose elevation','Equivalent doses of each pure'],flags=['vehicle_composition_changes_with_dose']+(['suspect_cmax_unit'] if compound!='paeonol' else []))
 specs[-1]['primary_eligible']=False
 if compound=='saikosaponin A':specs[-1]['auc_definition']='AUC0-8h'
# Preserve molar units. Do not normalize suspicious AUC nM hr/kg labels or invent per-timepoint n.
for pmc,table,species,strain,low in [('PMC3269340','d20e593','rat','Sprague-Dawley',1.22),('PMC3269340','d20e1470','mouse','C57Bl/6',1)]:
 for route,columns,doses in [('oral',[3,4],[low,10]),('subcutaneous',[5,6],[1,10])]:
  add(pmc,table,'memantine',doses,columns,None,'Observed Cmax',route=route,n_animals=[None,None],
      vehicle='physiological saline; dose volume changes',fasting='fasted; duration not reported',infusion_duration_min=0,
      evidence=['Since all time points were terminal','solution in physiological saline','2 mL/kg for both rats and mice'],
      flags=['dose_volume_changes','composite_PK','auc_unit_ambiguous','per_timepoint_n_unavailable'])
  specs[-1].update(species=species,strain=strain,cmax_unit='nM',primary_eligible=False,endpoint_origin='reported_composite_profile')
for spec in specs:
 if 'oral_vehicle_not_explicit' in spec['flags']:spec['primary_eligible']=False
Path('configs/curation/pmc_expansion.json').write_text(json.dumps(specs,indent=2,ensure_ascii=False)+'\n')
print(len(specs),'reviewed protocols')
