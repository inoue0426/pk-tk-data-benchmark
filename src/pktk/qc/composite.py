"""Representative, protocol-grounded audit of the 112 exploratory candidate groups."""
import json
from pktk.collectors.http import sha256
import pandas as pd
from pktk.standardization.core import identity

KEYS=['compound_id','species','sex','route','matrix','document_id','regimen']

def audit(root,endpoints,write_reports=True):
    cells=pd.read_csv(root/'data/processed/cvtdb_composite_candidates.csv',keep_default_na=False)
    counts=cells.groupby(KEYS).dose_mg_kg.nunique()
    keys=counts[counts>=2].reset_index()[KEYS]
    selection={'CvTdb:document:171':'S0592','CvTdb:document:172':'S0613','CvTdb:document:25297':'K10262'}
    specs=json.loads((root/'configs/curation/cebs_expansion.json').read_text())
    decisions={v['document_id']:v for v in json.loads((root/'configs/curation/composite_decisions.json').read_text())}
    for v in decisions.values():
        for ref in v['evidence']:
            if sha256(root/'data/raw'/ref['file'])!=ref['sha256']:raise ValueError('Audit evidence hash mismatch')
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
        if key.document_id in decisions:
            review=decisions[key.document_id]
            d.update(audit_status=review['status'],decision_reason=review['decision_reason'],
                     protocol_reviewed=review['protocol_reviewed'],evidence_json=json.dumps(review['evidence']))
            spec=next((v for v in specs if v['cvtdb_document_id']==key.document_id),None)
            if spec:
                e=endpoints[(endpoints.document_id==key.document_id)&(endpoints.compound_id==key.compound_id)&
                    (endpoints.species==key.species)&(endpoints.sex==key.sex)&(endpoints.route==key.route)&(endpoints.matrix==key.matrix)]
                o=cells[(cells.document_id==key.document_id)&(cells.compound_id==key.compound_id)&
                    (cells.species==key.species)&(cells.sex==key.sex)&(cells.route==key.route)&(cells.matrix==key.matrix)]
                if e.empty:
                    d.update(audit_status='withheld_no_reviewed_endpoint',decision_reason=spec['review_note']+' No reviewed endpoint stratum matches this candidate.')
                else:
                    # Endpoint rows are already protocol-separated; candidate itself may contain several formulations.
                    from pktk.expansion import build_series
                    trial=e.copy();trial['molecular_group_id']=trial.compound_id
                    series=build_series(trial);accepted=series[series.primary_eligible]
                    doses=set(pd.to_numeric(o.dose_mg_kg))
                    dose_match=all(any(abs(a-b)<1e-6 for b in doses) for a in e.dose_mg_kg)
                    d.update(audit_status='validated_group_endpoints' if dose_match and len(accepted) else 'withheld_dose_reconciliation',
                        matched_source_records_json=json.dumps(e.source_record_id.tolist()),
                        n_validated_endpoint_series=len(accepted),n_doses=e.dose_mg_kg.nunique(),
                        eligible_auc=bool(accepted.eligible_auc.any()),eligible_cmax=bool(accepted.eligible_cmax.any()),
                        source_doses_json=json.dumps(sorted(doses)),protocol_evidence_json=e.protocol_evidence_json.iloc[0],
                        study_id=spec['study'])
                    if not dose_match:d['decision_reason']+=' Source endpoints stand on their own; CvTdb dose linkage remains unresolved.'
        if 'protocol_reviewed' not in d:d['protocol_reviewed']=d['audit_status']!='not_reviewed'
        rows.append(d)
    out=pd.DataFrame(rows)
    if write_reports:out.to_csv(root/'reports/composite_audit.csv',index=False)
    summary=dict(total_candidates=len(out),assessed_candidates=int(out.audit_status.ne('not_reviewed').sum()),
                 reviewed_candidates=int(out.protocol_reviewed.sum()),
                 protocol_unresolved_candidates=int((~out.protocol_reviewed).sum()),
                 status_counts=out.audit_status.value_counts().to_dict(),
                 selection='All 112 candidates have evidence-based disposition or explicit unresolved protocol status; validation requires reviewed endpoint rows.',
                 expansion_decision='No individual-animal trajectories invented. Candidate groups and validated endpoint series have different granularity: formulations and strains may split a candidate.')
    if write_reports:(root/'reports/composite_audit.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary
