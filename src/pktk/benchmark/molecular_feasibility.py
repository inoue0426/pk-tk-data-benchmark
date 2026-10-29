"""Dataset A feasibility only: no model fitting or invented study independence."""
import json,hashlib
import pandas as pd

def linked_components(rows):
    parent={}
    def find(a):
        parent.setdefault(a,a)
        if parent[a]!=a:parent[a]=find(parent[a])
        return parent[a]
    def union(a,b):
        a,b=find(a),find(b)
        if a!=b:parent[max(a,b)]=min(a,b)
    for _,r in rows.iterrows():union('compound:'+r.connectivity_key,'publication:'+r.publication_id)
    return {key:hashlib.sha256(find('compound:'+key).encode()).hexdigest()[:16] for key in rows.connectivity_key.unique()}

def assess(good,out):
    # Unknown original publications cannot participate in a study-disjoint claim.
    known=good[good.publication_id.ne('')].copy();mapping=linked_components(known)
    known['dependency_component']=known.connectivity_key.map(mapping)
    known[['source_record_id','inchikey','connectivity_key','publication_id','dependency_component']].to_csv(out/'dependency_components.csv',index=False)
    rows=[]
    for keys,g in good.groupby(['endpoint','species','route','matrix','unit','endpoint_definition','dose_normalized'],sort=True):
        k=g[g.publication_id.ne('')]
        components=k.connectivity_key.map(mapping).value_counts()
        rows.append(dict(zip(['endpoint','species','route','matrix','unit','endpoint_definition','dose_normalized'],keys))|
          dict(compounds=g.inchikey.nunique(),connectivity_groups=g.connectivity_key.nunique(),records=len(g),
               source_publications=k.publication_id.nunique(),dependency_components=len(components),
               largest_component_record_fraction=float(components.max()/len(k)) if len(k) else None,
               missing_original_publication_records=int(g.publication_id.eq('').sum())))
    table=pd.DataFrame(rows);table.to_csv(out/'prediction_cohort_feasibility.csv',index=False)
    knowncounts=known.groupby('dependency_component').connectivity_key.nunique()
    result=dict(molecular_structures=good.inchikey.nunique(),connectivity_groups=good.connectivity_key.nunique(),
       scaffold_groups=good.scaffold.nunique(),acyclic_compounds=good[good.scaffold.eq('')].inchikey.nunique(),
       known_publication_dependency_components=len(knowncounts),largest_component_compounds=int(knowncounts.max()) if len(knowncounts) else 0,
       no_models_trained=True,no_evaluation_split_claimed=True,
       conclusion='Structure diversity exceeds the former cross-dose cohort; compound-disjoint endpoint modeling is structurally possible in populated, compatible cohorts. Original study independence is not established by publication IDs. No train/test split is frozen before source/protocol and duplicate review.',
       limitations=['ChEMBL scope is the first fixed activity-ID tranche, not random or exhaustive.',
        'PKSmart aggregate source rows lack original publication and administered dose; exclude from study-disjoint evaluation until mapped.',
        'Connectivity groups conservatively keep unspecified stereochemistry/charge variants together in future splits.',
        'Never pool AUC definitions, matrices, routes, dose-normalized and dose-specific values.',
        'Endpoint cohorts need dose, matrix, assay-method and species-appropriate review; AUC/Cmax at unknown doses cannot be unconditional molecular targets.',
        'Scaffold counts include one shared acyclic group; scaffold-disjoint performance has not been evaluated.'])
    (out/'feasibility.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
