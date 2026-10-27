import math
import pandas as pd

# Flags block derived endpoints; missing metadata also prevents high-confidence use separately.
BLOCKING = {'not_preclinical','nonparent_analyte','invalid_dose','invalid_concentration','invalid_time',
            'ambiguous_route','unsupported_matrix','duplicate_record','conflicting_dose_normalization',
            'conflicting_concentration_normalization','conflicting_time_normalization','not_single_dose',
            'radiolabeled_analyte','inconsistent_series_metadata'}

def check_observations(data):
    data = data.copy()
    data['eligible_for_curve'] = False
    duplicates = data.duplicated(['source','source_record_id'], keep=False)
    conflicts = set()
    metadata = ['compound_id','species','sex','matrix','route','dose_mg_kg','study_id','document_id','subject_id']
    for exp, g in data.groupby('experiment_id'):
        if any(g[c].nunique(dropna=False)>1 for c in metadata):
            conflicts.add(exp)
    for i,r in data.iterrows():
        flags = set(filter(None,r.qc_flags.split(';')))
        if r.species not in {'rat','mouse','hamster','dog','monkey','rabbit'}: flags.add('not_preclinical')
        if not r.parent_analyte: flags.add('nonparent_analyte')
        if not math.isfinite(r.dose_mg_kg) or r.dose_mg_kg<=0: flags.add('invalid_dose')
        if math.isfinite(r.dose_mg_kg) and r.dose_mg_kg>10000: flags.add('extreme_dose_review')
        if not math.isfinite(r.concentration_ug_ml) or r.concentration_ug_ml<0: flags.add('invalid_concentration')
        if r.concentration_ug_ml==0: flags.add('zero_concentration_censoring_unknown')
        if not math.isfinite(r.time_h) or r.time_h<0: flags.add('invalid_time')
        if r.route not in {'oral','intravenous'}: flags.add('ambiguous_route')
        if r.matrix not in {'blood','plasma','serum'}: flags.add('unsupported_matrix')
        if not r.sex: flags.add('missing_sex')
        if r.regimen!='once': flags.add('not_single_dose')
        if r.radiolabeled in {'1','1.0','True'}: flags.add('radiolabeled_analyte')
        if not r.radiolabeled: flags.add('radiolabel_status_unknown')
        if duplicates.loc[i]: flags.add('duplicate_record')
        if r.experiment_id in conflicts: flags.add('inconsistent_series_metadata')
        if not r.publication_doi: flags.add('missing_doi')
        data.at[i,'qc_flags'] = ';'.join(sorted(flags))
        data.at[i,'eligible_for_curve'] = not bool(flags & BLOCKING)
    return data

def flag_counts(data):
    values=data.qc_flags.str.split(';').explode()
    return values[values.ne('')].value_counts().to_dict()

def overlap_candidates(observations, endpoints):
    """Flag shared publication+name or molecular ID across sources, not automatic deduplication."""
    cols = ['source','compound_id','compound_name','publication_doi']
    x = pd.concat([observations[cols],endpoints[cols]],ignore_index=True).drop_duplicates()
    x['name_key'] = x.compound_name.str.lower().str.replace(r'[^a-z0-9]','',regex=True)
    out=[]
    for key in ['compound_id','name_key']:
        for (doi,identity),g in x.groupby(['publication_doi',key]):
            if doi and identity and g.source.nunique()>1:
                out.append({'doi':doi,'match_type':key,'identity':identity,'sources':sorted(g.source.unique())})
    return out
