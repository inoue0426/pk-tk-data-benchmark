"""Conservative molecular PK harmonization. Unknown context is never imputed."""
import math,re
from functools import lru_cache
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

@lru_cache(maxsize=20000)
def structure(smiles):
    m=Chem.MolFromSmiles(smiles) if smiles else None
    if m is None:return dict(structure_status='invalid_or_missing_structure')
    result=dict(original_smiles=smiles,source_inchikey=Chem.MolToInchiKey(m),molecular_weight=Descriptors.MolWt(m))
    if len(Chem.GetMolFrags(m))!=1:return result|dict(structure_status='multicomponent_structure_unresolved')
    # Preserve stereochemistry; remove reversible protonation only. No largest-fragment guess.
    n=rdMolStandardize.Uncharger().uncharge(rdMolStandardize.Cleanup(m))
    key=Chem.MolToInchiKey(n)
    if not key:return result|dict(structure_status='invalid_inchikey')
    return result|dict(structure_status='resolved',canonical_smiles=Chem.MolToSmiles(n,isomericSmiles=True),
        inchikey=key,connectivity_key=key.split('-')[0],scaffold=MurckoScaffold.MurckoScaffoldSmiles(mol=n,includeChirality=False))

def context(description):
    s=description.lower();flags=[]
    routes=[]
    patterns={'intravenous':r'\bintravenous\w*|\bi\.?v\.?\b','oral':r'\boral\w*|\bperoral\w*|\bp\.?o\.?\b|\bgavage\b',
              'subcutaneous':r'\bsubcutaneous\w*|\bs\.?c\.?\b','intraperitoneal':r'\bintraperitoneal\w*|\bi\.?p\.?\b',
              'intraduodenal':r'\bintraduodenal\w*','intramuscular':r'\bintramuscular\w*'}
    for route,pattern in patterns.items():
        if re.search(pattern,s):routes.append(route)
    route=routes[0] if len(routes)==1 else ''
    if len(routes)>1:flags.append('multiple_routes')
    if not routes:flags.append('missing_route')
    mats=[m for m in ['plasma','blood','serum','brain','csf','liver','urine'] if re.search(r'\b'+m+r'\b',s)]
    matrix=mats[0] if len(mats)==1 else ''
    if len(mats)>1:flags.append('multiple_matrices')
    if not mats:flags.append('missing_matrix')
    if 'portal vein' in s or 'portal blood' in s:
        matrix='portal_blood'
        flags.append('non_systemic_sampling_site')
    auc_definition='AUC_unspecified'
    if re.search(r'(?:auc|area under).{0,100}(?:infinity|infinite|\binf\b)',s):auc_definition='AUC0-inf'
    windows=re.findall(r'\b0\s*(?:-|–|to)\s*(\d+(?:\.\d+)?)\s*(?:hr?\b|hours?\b)',s)
    if len(set(windows))==1 and ('auc' in s or 'area under' in s):auc_definition='AUC0-'+str(float(windows[0]))+'h'
    normalized=bool(re.search(r'dose[- ]normali[sz]|normali[sz]\w*.{0,35}dose|(?:auc|cmax)\s*/\s*dose',s))
    matches=re.findall(r'(?<![\w.])(\d+(?:\.\d+)?)\s*\(?\s*(?:mg\s*/\s*kg|mg\s*kg\s*[-−]1)',s)
    doses=sorted(set(float(d) for d in matches))
    # A single numeric match in a range is not a single administered dose.
    ambiguous=bool(re.search(r'\d\s*(?:-|–|to|and|or|,)\s*\d+(?:\.\d+)?\s*\(?\s*mg\s*(?:/\s*kg|kg\s*[-−]1)',s))
    dose=doses[0] if len(doses)==1 and not ambiguous and not normalized else None
    if dose is None:flags.append('missing_or_ambiguous_administered_dose')
    if normalized:flags.append('dose_normalized_not_dose_specific')
    invivo=bool(re.search(r'\bin vivo\b',s) or routes)
    non_pk=bool(re.search(r'in vitro|ex vivo|effect[ -]time|effect.{0,10}curve|response[ -]time|metabolite|radioactiv|microsom|hepatocyt|intrinsic|homogenate|\bs9\b|cell line|liver slice|perfus|predicted|simulated|calculated from|\bdiuresis\b|pressor|blood pressure|camp levels|enzyme|receptor occupancy',s))
    return dict(route=route,matrix=matrix,dose_mg_kg=dose,dose_normalized='yes' if normalized else 'not_reported',
        in_vivo_evidence=invivo and not non_pk,context_flags=flags,
        auc_definition=auc_definition,
        sex='female' if 'female' in s else 'male' if re.search(r'\bmale\b',s) else '')

# ChEMBL normalized units; no guesses for malformed dimensions.
UNITS={('CL','mL.min-1.kg-1'):('clearance','mL/min/kg',1.),('Vd','L.kg-1'):('volume_distribution','L/kg',1.),
 ('T1/2','hr'):('half_life','h',1.),('AUC','ng.hr.mL-1'):('auc','ug*h/mL',.001),
 ('Cmax','ug.mL-1'):('cmax','ug/mL',1.)}

def endpoint(kind,unit,value,mw=None):
    try:v=float(value)
    except (TypeError,ValueError):return None,'missing_numeric_endpoint'
    if not math.isfinite(v) or v<=0:return None,'nonpositive_or_nonfinite_endpoint'
    if kind=='Cmax' and unit=='nM' and mw:
        return dict(endpoint='cmax',value=v*mw/1e6,unit='ug/mL',conversion='nM_to_mass_using_source_structure_MW'),'molar_mass_conversion'
    if (kind,unit) not in UNITS:return None,'unsupported_or_conflicting_unit'
    ep,u,f=UNITS[kind,unit]
    return dict(endpoint=ep,value=v*f,unit=u,conversion=f'multiply_by_{f:g}'),''
