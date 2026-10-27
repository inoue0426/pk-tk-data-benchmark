"""Conservative mass-based normalization; no implicit molar or body-weight conversion."""
import hashlib
import json
import math
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

MISSING = {'', 'nr', 'na', 'n/a', 'nan', 'none', 'not reported', '<na>'}

def clean(x):
    if x is None or pd.isna(x):
        return ''
    s = str(x).strip()
    return '' if s.lower() in MISSING else s

def number(x):
    try:
        n = float(x)
        return n if math.isfinite(n) else float('nan')
    except (ValueError, TypeError):
        return float('nan')

def unit(x):
    return clean(x).lower().replace('μ','u').replace('µ','u').replace(' ','')

def convert(value, units, kind):
    tables = {
        'dose': {'mg/kg':1., 'mg/kgbw':1., 'ug/kg':.001, 'g/kg':1000.},
        'concentration': {'ug/ml':1., 'ng/ml':.001, 'mg/l':1., 'ug/l':.001, 'mg/ml':1000., 'ng/ul':1.},
        'time': {'h':1., 'hr':1., 'hours':1., 'hour':1., 'min':1/60, 'minutes':1/60, 'minute':1/60, 's':1/3600, 'sec':1/3600, 'day':24., 'days':24.},
    }
    return number(value) * tables[kind].get(unit(units), float('nan'))

def identity(parts):
    return hashlib.sha256(json.dumps(list(parts), default=str, ensure_ascii=True).encode()).hexdigest()[:20]

def structure(smiles):
    """Retain stereochemistry and fragments; do not merge salts/tautomers automatically."""
    if not clean(smiles):
        return dict(canonical_smiles='', inchikey='', scaffold='', molecular_weight=None, structure_flag='missing_structure')
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return dict(canonical_smiles='', inchikey='', scaffold='', molecular_weight=None, structure_flag='invalid_structure')
    return dict(canonical_smiles=Chem.MolToSmiles(mol, isomericSmiles=True), inchikey=Chem.MolToInchiKey(mol),
                scaffold=MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False),
                molecular_weight=Descriptors.MolWt(mol), structure_flag='')

def harmonize_cvtdb(raw, provenance):
    rows = []
    for original in raw.to_dict('records'):
        r = {k: clean(v) for k,v in original.items()}
        cid = r['analyzed_chem_dtxsid'] or 'CvTdb:chemical:'+r['fk_analyzed_chemical_id']
        dose = convert(r['dose_level_normalized'], r['dose_level_units_normalized'], 'dose')
        conc = convert(r['conc'], r['conc_units_normalized'], 'concentration')
        flags = ['missing_strain','unverified_protocol','missing_structure']
        if not r['analyzed_chem_dtxsid']:
            flags.append('missing_external_identifier')
        for normalized, value, units, kind in [(dose,r['dose_level_original'],r['dose_level_units_original'],'dose'),(conc,r['conc_original'],r['conc_units_original'],'concentration'),(number(r['time_hr']),r['time_original'],r['time_units_original'],'time')]:
            converted = convert(value, units, kind)
            if math.isfinite(converted) and math.isfinite(normalized) and not math.isclose(converted, normalized, rel_tol=.01, abs_tol=1e-5):
                flags.append('conflicting_'+kind+'_normalization')
        sex = r['sex'].lower().replace('male and female','mixed')
        record = dict(source='cvtdb', source_record_id=r['conc_time_id'],
            compound_id=cid, compound_name=r['analyzed_chem_name'] or r['analyzed_chem_name_original'],
            casrn=r['analyzed_chem_casrn'], canonical_smiles='', inchikey='', species=r['species'].lower(),
            strain='', sex=sex, route={'oral':'oral','iv':'intravenous'}.get(r['administration_route_normalized'],''),
            dose_mg_kg=dose, dose_unit='mg/kg', time_h=number(r['time_hr']), time_unit='h',
            concentration_ug_ml=conc, concentration_unit='ug/mL', matrix=r['conc_medium_normalized'].lower(),
            study_id='CvTdb:study:'+r['fk_study_id'], document_id='CvTdb:document:'+r['fk_extraction_document_id'],
            experiment_id='CvTdb:series:'+r['fk_series_id'], subject_id=r['fk_subject_id'],
            n_animals=number(r['n_subjects_in_series']), publication_doi=r['doi'].lower(), publication_url=r['url'],
            parent_analyte=r['fk_dosed_chemical_id']==r['fk_analyzed_chemical_id'],
            regimen=r['dose_frequency_normalized'].lower(), vehicle=r['dose_vehicle'].lower(),
            fasting=r['fasting_period'].lower(), age=r['age'], age_unit=r['age_units'],
            dose_duration=r['dose_duration'], dose_duration_unit=r['dose_duration_units'],
            radiolabeled=r['radiolabeled'], source_url=provenance['url'], retrieved_at=provenance['retrieved_at'],
            dataset_version=provenance['version'], raw_sha256=provenance['sha256'],
            original_json=json.dumps(r, sort_keys=True), qc_flags=';'.join(flags))
        rows.append(record)
    return pd.DataFrame(rows)
