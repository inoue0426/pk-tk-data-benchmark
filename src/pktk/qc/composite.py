"""Representative, protocol-grounded audit of the 112 exploratory candidate groups."""
import json
import pandas as pd
from pktk.standardization.core import identity

KEYS=['compound_id','species','sex','route','matrix','document_id','regimen']

def audit(root,endpoints):
    cells=pd.read_csv(root/'data/processed/cvtdb_composite_candidates.csv',keep_default_na=False)
    counts=cells.groupby(KEYS).dose_mg_kg.nunique()
    keys=counts[counts>=2].reset_index()[KEYS]
    selection={'CvTdb:document:171':'S0592','CvTdb:document:172':'S0613','CvTdb:document:25297':'K10262'}
    observations=pd.read_csv(root/'data/processed/observations.csv',keep_default_na=False,low_memory=False)
    rows=[]
    for _,key in keys.iterrows():
        d=key.to_dict();d['candidate_id']=identity(list(key));d['audit_status']='not_reviewed'
        d['decision_reason']='Outside initial stratified protocol sample'
        if key.document_id in selection:
            e=endpoints[(endpoints.document_id==key.document_id)&(endpoints.species==key.species)&(endpoints.sex==key.sex)&(endpoints.route==key.route)]
            o=observations[(observations.document_id==key.document_id)&(observations.species==key.species)&(observations.sex==key.sex)&(observations.route==key.route)&(observations.matrix==key.matrix)]
            d.update(study_id=selection[key.document_id],n_observations=len(o),n_subject_ids=o.subject_id.nunique(),
                     n_source_series=o.experiment_id.nunique(),source_doses_json=json.dumps(sorted(pd.to_numeric(o.dose_mg_kg).unique().tolist())))
            if key.matrix!='plasma':
                d.update(audit_status='withheld_matrix_conflict',decision_reason='CEBS source specifies plasma; CvTdb also labels rows blood. No blood-plasma pooling or silent correction.')
            elif e.empty:
                d.update(audit_status='withheld_no_matching_endpoint',decision_reason='No reviewed group endpoint for candidate')
            elif not e.primary_eligible.all():
                d.update(audit_status='withheld_dose_conflict',decision_reason='Isoeugenol rat middle-dose field says 35 mg/kg; group label and CvTdb say 70. Do not guess the dose.')
            else:
                d.update(audit_status='validated_group_endpoints',decision_reason='Source explicitly describes sparse sampling and group analysis. Use CEBS reported endpoints; no invented individual trajectories.',
                         n_doses=e.dose_mg_kg.nunique(),eligible_auc=bool(e.high_confidence_auc.all()),eligible_cmax=bool(e.high_confidence_cmax.all()),
                         matched_source_records_json=json.dumps(e.source_record_id.tolist()),protocol_evidence_json=e.protocol_evidence_json.iloc[0])
                # Every accepted dose must match a source dose in the existing candidate observations.
                if not set(e.dose_mg_kg).issubset(set(pd.to_numeric(o.dose_mg_kg))):
                    d.update(audit_status='withheld_dose_reconciliation',decision_reason='Actual endpoint doses do not all match CvTdb dose cells')
        rows.append(d)
    out=pd.DataFrame(rows);out.to_csv(root/'reports/composite_audit.csv',index=False)
    summary=dict(total_candidates=len(out),reviewed_candidates=int(out.audit_status.ne('not_reviewed').sum()),
                 status_counts=out.audit_status.value_counts().to_dict(),selection='Purposive stratified sample: three source studies, rat/mouse, sex, route and matrix; not a statistical yield estimate',
                 expansion_decision='Initial recovery supports targeted CEBS extraction. No blanket promotion; reconcile dose and matrix errors first. Audit stays at three source studies for this run.')
    (root/'reports/composite_audit.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary
