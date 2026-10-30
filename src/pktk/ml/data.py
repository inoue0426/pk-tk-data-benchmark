"""Conservative conditional PK task definitions and structure/study split audits."""
from pathlib import Path
import hashlib,json,re
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'data/processed/ml'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dump(path,obj):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(obj,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x),allow_nan=False)+'\n')


def structure(smiles):
    m=Chem.MolFromSmiles(smiles)
    if m is None or len(Chem.GetMolFrags(m))!=1:raise ValueError('Invalid/multicomponent structure')
    # A conservative split identity, not a replacement for source molecular identity.
    p=rdMolStandardize.Uncharger().uncharge(Chem.Mol(m))
    enumerator=rdMolStandardize.TautomerEnumerator()
    tautomers=enumerator.Enumerate(p)
    if str(tautomers.status)!='Completed':raise ValueError('incomplete_tautomer_enumeration')
    p=enumerator.PickCanonical(tautomers)
    Chem.RemoveStereochemistry(p)
    family=Chem.MolToInchiKey(p).split('-')[0]
    scaffold=MurckoScaffold.MurckoScaffoldSmiles(mol=p,includeChirality=False) or '__ACYCLIC__'
    return family,scaffold,len(Chem.FindMolChiralCenters(m,includeUnassigned=True)),sum(x[1]=='?' for x in Chem.FindMolChiralCenters(m,includeUnassigned=True))


def exclusion(row):
    reasons=[]
    if row.source!='chembl':reasons.append('not_original_publication_annotated_chembl')
    if row.species!='rat':reasons.append('species_not_rat_task')
    if row.matrix!='plasma':reasons.append('matrix_not_explicit_plasma')
    if not row.publication_id:reasons.append('missing_publication')
    if row.endpoint not in ['clearance','cmax']:reasons.append('no_compatible_sufficient_task_definition')
    wanted='intravenous' if row.endpoint=='clearance' else 'oral'
    if row.route!=wanted:reasons.append('route_not_task_route')
    dose=pd.to_numeric(row.dose_mg_kg,errors='coerce')
    if not np.isfinite(dose) or dose<=0:reasons.append('missing_positive_administered_dose')
    doses=re.findall(r'(\d+(?:\.\d+)?)\s*mg\s*/\s*kg(?!\s*/)',row.description,re.I)
    if not doses or len(set(doses))!=1 or not np.isclose(float(doses[0]),dose):reasons.append('administered_mg_kg_not_independently_confirmed')
    if row.dose_normalized=='yes':reasons.append('dose_normalized_not_absolute_endpoint')
    desc=row.description.lower()
    if re.search(r'unbound|uncharged|intrinsic|renal clearance|hepatic clearance|steady.state|infus|diabet|diseas|knock.?out|pretreat|co.?administ|repeated|multiple dos',desc):
        reasons.append('special_analyte_regimen_or_organ_clearance')
    if row.endpoint=='cmax' and re.search(r'after\s+(?:1\s*h|60\s*min)',desc):
        reasons.append('fixed_time_concentration_ambiguity')
    if row.endpoint=='clearance' and not re.search(r'plasma.{0,30}clearance|clearance.{0,40}plasma|clearance clp',desc):
        reasons.append('systemic_plasma_clearance_not_explicit')
    expected='mL/min/kg' if row.endpoint=='clearance' else 'ug/mL'
    if row.unit!=expected:reasons.append('incompatible_unit')
    return sorted(set(reasons))


def connected_groups(d,kind):
    parents={i:i for i in d.index}
    def root(i):
        while parents[i]!=i:parents[i]=parents[parents[i]];i=parents[i]
        return i
    def union(a,b):
        a,b=root(a),root(b)
        if a!=b:parents[max(a,b)]=min(a,b)
    cols=['family','connectivity_key']
    if 'scaffold' in kind:cols.append('ml_scaffold')
    if 'study' in kind:cols.append('publication_id')
    for col in cols:
        for _,g in d.groupby(col):
            for i in g.index[1:]:union(g.index[0],i)
    return pd.Series({i:str(root(i)) for i in d.index})


def allocate(d,kind,seed,config):
    groups=connected_groups(d,kind)
    # Randomized group ordering, greedy allocation to remaining row quotas. No outcome access.
    rng=np.random.default_rng(seed)
    ids=np.array(sorted(groups.unique()));rng.shuffle(ids)
    roles=list(config['fractions']);need={r:len(d)*config['fractions'][r] for r in roles}
    assigned={}
    for g in ids:
        n=int((groups==g).sum());role=max(roles,key=lambda r:need[r])
        assigned[g]=role;need[role]-=n
    out=d[['row_id','inchikey','family','connectivity_key','ml_scaffold','publication_id']].copy()
    out['split_group']=groups;out['partition']=groups.map(assigned)
    for col in ['family','connectivity_key']+(['ml_scaffold'] if 'scaffold' in kind else [])+(['publication_id'] if 'study' in kind else []):
        if out.groupby(col).partition.nunique().max()>1:raise ValueError('Split leakage: '+col)
    sizes={r:int(out[out.partition==r].family.nunique()) for r in roles}
    valid=all(sizes[r]>=config['minimum_partition_compounds'][r] for r in roles)
    meta=dict(task=d.task.iloc[0],split=kind,seed=seed,eligible=valid,partition_compounds=sizes,
              partition_records=out.groupby('partition').size().to_dict(),n_groups=groups.nunique(),largest_group_records=int(groups.value_counts().max()),
              studies_crossing_partitions=int((out.groupby('publication_id').partition.nunique()>1).sum()),
              scaffolds_crossing_partitions=int((out.groupby('ml_scaffold').partition.nunique()>1).sum()))
    return out,meta


def build():
    OUT.mkdir(parents=True,exist_ok=True)
    source=ROOT/'data/processed/molecular/eligible_endpoints.csv'
    cfg=json.loads((ROOT/'configs/ml/benchmark.json').read_text())
    raw=pd.read_csv(source,keep_default_na=False)
    rows=[];selected=[];chem={}
    for r in raw.itertuples(index=False):
        reasons=exclusion(r)
        if not reasons:
            try:
                if r.canonical_smiles not in chem:chem[r.canonical_smiles]=structure(r.canonical_smiles)
            except Exception as e:reasons.append('structure_failure:'+str(e))
        rows.append(dict(source_record_id=r.source_record_id,endpoint=r.endpoint,source=r.source,eligible_ml=not reasons,reasons=';'.join(reasons)))
        if reasons:continue
        family,scaffold,centers,unknown=chem[r.canonical_smiles]
        row=r._asdict();row.update(family=family,ml_scaffold=scaffold,stereocenters=centers,unspecified_stereocenters=unknown,
                                  task='rat_iv_plasma_clearance' if r.endpoint=='clearance' else 'rat_oral_plasma_cmax',
                                  label_kind='source_reported_PK_endpoint_not_individual_animal',log_dose=np.log10(float(r.dose_mg_kg)),y=np.log10(float(r.value)))
        selected.append(row)
    pd.DataFrame(rows).to_csv(OUT/'eligibility.csv',index=False)
    d=pd.DataFrame(selected)
    # Same full identity / publication / dose / endpoint: median in log space; retain every member ID.
    dedup=[]
    for _,g in d.groupby(['task','inchikey','publication_id','dose_mg_kg'],sort=True):
        r=g.iloc[0].to_dict();r['y']=float(g.y.median());r['value']=10**r['y']
        r['original_record_ids_json']=json.dumps(sorted(g.source_record_id));r['original_records']=len(g)
        r['within_cell_log_range']=float(g.y.max()-g.y.min())
        r['row_id']=hashlib.sha256('|'.join(sorted(g.source_record_id)).encode()).hexdigest()[:20]
        dedup.append(r)
    d=pd.DataFrame(dedup).sort_values('row_id').reset_index(drop=True)
    d.to_csv(OUT/'tasks.csv',index=False)
    # Freeze all original A/B artifacts, numerical benchmark outputs and source manifests.
    snapshot={str(p.relative_to(ROOT)):sha(p) for folder in ['data/processed/molecular','results/research'] for p in sorted((ROOT/folder).glob('*')) if p.is_file()}
    for name in ['data/processed/high_confidence_endpoints.csv','data/processed/crossdose_series_expanded.csv','reports/molecular_summary.json','reports/research_manifest.json']:
        snapshot[name]=sha(ROOT/name)
    preserved=ROOT/'reports/ml/original_artifacts.json'
    if preserved.exists():
        assert json.loads(preserved.read_text())==snapshot,'Original A/B artifacts changed'
    else:dump(preserved,snapshot)
    splits=[];metadata=[]
    for task,g in d.groupby('task'):
        for kind in cfg['splits']:
            for seed in cfg['seeds']:
                split,info=allocate(g,kind,seed,cfg);split['task']=task;split['split']=kind;split['seed']=seed
                splits.append(split);metadata.append(info)
    pd.concat(splits).to_csv(OUT/'splits.csv',index=False)
    dump(OUT/'split_audit.json',metadata)
    counts=d.groupby(['task','endpoint','species','route','matrix','unit']).agg(records=('row_id','size'),original_records=('original_records','sum'),compounds=('inchikey','nunique'),families=('family','nunique'),publications=('publication_id','nunique')).reset_index().to_dict('records')
    summary=dict(input_sha256=sha(source),config_sha256=sha(ROOT/'configs/ml/benchmark.json'),rows_considered=len(raw),records_eligible=sum(r['eligible_ml'] for r in rows),compounds=d.inchikey.nunique(),tasks=counts,eligible_split_runs=sum(m['eligible'] for m in metadata),
                 missing_stereochemistry_rows=int((d.unspecified_stereocenters>0).sum()),
                 exclusions=pd.DataFrame(rows).groupby('endpoint').eligible_ml.agg(['size','sum']).reset_index().to_dict('records'))
    dump(ROOT/'reports/ml/audit.json',summary);print(json.dumps(summary,indent=2));print(pd.DataFrame(metadata).to_string(index=False))

if __name__=='__main__':build()
