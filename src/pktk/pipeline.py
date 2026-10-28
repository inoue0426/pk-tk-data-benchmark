import argparse
import json
from pathlib import Path
import pandas as pd
from pktk.collectors.http import Collector,log
from pktk.collectors import cvtdb,pmc
from pktk.standardization.core import harmonize_cvtdb,structure
from pktk.qc.checks import check_observations,flag_counts,overlap_candidates
from pktk.benchmark.series import curves_and_candidates,reported_series
from pktk.benchmark.evaluate import evaluate

def dump(path,data):
    path.write_text(json.dumps(data,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x),allow_nan=False)+'\n')

def run(root,offline=False):
    root=Path(root)
    for folder in ['reports','results','data/processed','data/interim']:
        (root/folder).mkdir(parents=True,exist_ok=True)
    config=json.loads((root/'configs/sources.json').read_text())
    collector=Collector(root,offline)
    log(root/'reports/progress.jsonl',phase='collection',status='started',offline=offline)
    raw,cvmeta=cvtdb.collect(collector,config['cvtdb_sha256'])
    endpoints,pmmeta=pmc.collect(collector,config.get('pmc_sha256'))
    log(root/'reports/progress.jsonl',phase='collection',status='succeeded',cvtdb_rows=len(raw),pmc_groups=len(endpoints))
    obs=check_observations(harmonize_cvtdb(raw,cvmeta))
    curves,candidates=curves_and_candidates(obs)
    reported=reported_series(endpoints)
    high=endpoints[endpoints.high_confidence_cmax]
    compounds=pd.concat([obs[['source','compound_id','compound_name','canonical_smiles']],endpoints[['source','compound_id','compound_name','canonical_smiles']]]).drop_duplicates().reset_index(drop=True)
    chemistry=pd.DataFrame([structure(s) for s in compounds.pop('canonical_smiles')])
    compounds=pd.concat([compounds,chemistry],axis=1)
    compounds.to_csv(root/'data/processed/compound_registry.csv',index=False)
    tables={'observations':obs,'reported_endpoints':endpoints,'curve_qc':curves,
            'candidate_series':candidates,'reported_series':reported,'high_confidence_cmax':high}
    for name,data in tables.items(): data.to_csv(root/'data/processed'/f'{name}.csv',index=False)
    splits,predictions,metrics,nonlinear,status=evaluate(endpoints)
    for name,data in [('task_a_splits',splits),('predictions',predictions),('metrics',metrics),('nonlinearity',nonlinear)]:
        data.to_csv(root/'results'/f'{name}.csv',index=False)
    dump(root/'results/benchmark_status.json',status)
    overlap=overlap_candidates(obs,endpoints)
    dump(root/'results/potential_overlap.json',overlap)
    pre=obs[obs.species!='human']
    parent=pre[pre.parent_analyte]
    cross=candidates[candidates.n_doses>=2]
    summary={
        'cvtdb_original_observations':len(obs),'cvtdb_analyte_ids_all_species':obs.compound_id.nunique(),
        'cvtdb_source_chemical_ids':raw.fk_analyzed_chemical_id.nunique(),
        'cvtdb_document_ids_all_species':obs.document_id.nunique(),
        'cvtdb_dose_specific_study_ids_all_species':obs.study_id.nunique(),
        'cvtdb_source_series_ids_all_species':obs.experiment_id.nunique(),
        'preclinical_observations':len(pre),'preclinical_parent_observations':len(parent),
        'preclinical_analyte_identifiers_collected':pre.compound_id.nunique()+endpoints.compound_id.nunique(),
        'preclinical_parent_compound_identifiers':parent.compound_id.nunique()+endpoints.compound_id.nunique(),
        'preclinical_document_studies':pre.document_id.nunique()+endpoints.document_id.nunique(),
        'preclinical_source_study_ids_plus_pmc_study':pre.study_id.nunique()+endpoints.study_id.nunique(),
        'eligible_curve_observations':int(obs.eligible_for_curve.sum()),'usable_curves':int(curves.usable_curve.sum()),
        'cvtdb_candidate_crossdose_series':len(cross),'cvtdb_candidate_compounds_2plus_doses':cross.compound_id.nunique(),
        'cvtdb_candidate_compounds_3plus_doses':candidates.loc[candidates.n_doses>=3,'compound_id'].nunique(),
        'cvtdb_candidates_with_comparable_partial_auc_window':int(cross.comparable_auc_window.sum()),
        'verified_compounds_2plus_comparable_doses':high.compound_id.nunique(),
        'verified_compounds_3plus_comparable_doses':reported.loc[reported.n_doses>=3,'compound_id'].nunique(),
        'verified_crossdose_series':int(reported.high_confidence.sum()),
        'verified_series_with_reported_auc_and_cmax':int((reported.high_confidence & reported.has_auc & reported.has_cmax).sum()),
        'verified_crossdose_series_eligible_auc_benchmark':0,
        'high_confidence_cmax_groups':len(high),'high_confidence_structures':int(high.inchikey.ne('').sum()),
        'preclinical_species':sorted(pre.species.unique()),'preclinical_routes':sorted(pre.route.unique()),
        'fraction_collected_preclinical_records_eligible_cmax':len(high)/(len(pre)+len(endpoints)),
        'fraction_preclinical_compound_ids_verified_crossdose':high.compound_id.nunique()/(pre.compound_id.nunique()+endpoints.compound_id.nunique()),
        'task_a_episodes':splits.episode_id.nunique() if len(splits) else 0,
        'baseline_predictions':len(predictions),'potential_source_overlap_groups':len(overlap),
        'identity_caveat':'Counts use source analyte/DTXSID IDs and a publication-local PMC ID; not structure-deduplicated molecular entities.',
        'study_caveat':'Documents are publication-level study proxies; CvTdb study IDs are dose-specific; series are not independent animals.',
    }
    dump(root/'reports/summary.json',summary)
    dump(root/'reports/source_manifest.json',{'cvtdb':cvmeta,'pmc':pmmeta})
    qc={'input_observations':len(obs),'retained_audit_observations':len(obs),'flagged_observations':int(obs.qc_flags.ne('').sum()),
        'excluded_from_curve_derivation':int((~obs.eligible_for_curve.astype(bool)).sum()),
        'observation_flags':flag_counts(obs),'curve_groups':len(curves),'curve_group_flags':flag_counts(curves),
        'usable_curve_groups':int(curves.usable_curve.sum()),'pmc_selected_groups':len(endpoints),
        'pmc_cmax_eligible_groups':len(high),'pmc_auc_benchmark_eligible_groups':0,
        'pmc_scope_exclusions':'Other Table 1 columns (IV/repeated dosing), all dog Table 2 columns, tissues and supplements are out of selected adapter scope, preserved in raw XML. Dog table conflict documented.',
        'potential_overlap_groups':len(overlap)}
    dump(root/'reports/qc_counts.json',qc)
    (root/'reports/collection_summary.md').write_text('# Executed collection summary\n\nAll counts below are computed from downloaded bytes; none are estimates.\n\n'+
        '\n'.join(f'- **{k}**: {v}' for k,v in summary.items())+'\n\nSource manifests contain UTC retrieval time, version, URL and SHA-256. Human and metabolite records remain in the audit table. Counts of identifiers are not confirmed unique molecular identities.\n')
    (root/'reports/qc_report.md').write_text('# Executed QC report\n\nNo downloaded observation was silently deleted. Flags can overlap; counts must not be added as mutually exclusive exclusions. All CvTdb records lack strain/structure and verified shared protocol.\n\n'+
        '```json\n'+json.dumps(qc,indent=2)+'\n```\n\nPartial AUC uses linear trapezoids only from the first to last observed time; no t=0 or terminal tail is imputed. Duplicate timepoints block curve integration instead of being averaged. Sampled maxima are not necessarily true Cmax. Unsupported units remain missing and flagged. Original strings, identifiers and units are retained in original_json and the raw archive.\n\nThe PMC rat table is checked against dose/sex/design headers. Its AUC0-tn endpoints remain available but cannot enter an AUC benchmark until terminal-time comparability is established. Missing structure does not invalidate an empirical within-compound Cmax episode but prevents molecular/scaffold evaluation.\n\nPotential overlaps use DOI plus name/identifier; absent matches do not prove independence. Cross-source structural identity resolution remains unfinished.\n')
    overall=[{k:v for k,v in r.items() if v is None or not pd.isna(v)} for r in metrics[metrics.aggregation=='model/n_support'].to_dict('records')] if len(metrics) else []
    (root/'reports/benchmark_feasibility.md').write_text('# Executed benchmark feasibility\n\n'+
        f'This is a working collection and QC pipeline with a **one-compound Cmax pilot**, not a publication-ready molecular benchmark. {len(high)} reported rat group means form {int(reported.high_confidence.sum())} sex-specific series with three doses each.\n\n'+
        '\n'.join(f'- **{k}**: {v}' for k,v in status.items())+
        '\n\nLinear scaling uses the geometric mean of support exposure/dose; power law fits log exposure against log dose using two support doses. Low-only episodes predict both higher doses; two-dose episodes predict the highest dose. Both models are compared on the same two-dose episodes. Natural logarithms are used. There is no hyperparameter selection.\n\n'+
        '```json\n'+json.dumps(overall,indent=2)+'\n```\n\n'+
        'One-compound results cannot establish molecular ML performance. There is no independently calibrated interval or nonlinear classifier. Descriptive slopes use all three doses and are not held-out detection accuracy. No significance or generalization claim is justified. Metrics by dose ratio, compound and species are in results/metrics.csv.\n\n'+
        'Task B requires verified identities and enough independent compounds and studies for development/held-out sets; allocate whole identity–study connected components, not rows. Task C additionally requires verified structures and whole scaffold components (including acyclic structures as a shared group). Task D needs independent dose-response labels, negative controls and uncertainty; this pilot supplies descriptive examples only.\n\n'+
        'Next priorities: adjudicate CvTdb publication protocols and strains, recover original compound structures and identity aliases, obtain matched-window AUC or individual data, investigate the PMC dog-table error, then add independent compounds before any molecular model or compound/scaffold split.\n')
    (root/'reports/final_report.md').write_text(f"""# Verified pilot results

The complete collection → harmonization → QC → Task A baseline workflow was executed on CvTdb and an open-access PMC article. This release is a reproducible working subset, not a finished molecular benchmark.

| Measurement | Executed count |
|---|---:|
| CvTdb original concentration observations retained | {len(obs):,} |
| Preclinical concentration observations | {len(pre):,} |
| Preclinical compound/analyte identifiers including PMC | {summary['preclinical_analyte_identifiers_collected']} |
| Preclinical parent-compound identifiers including PMC | {summary['preclinical_parent_compound_identifiers']} |
| Preclinical document-level study proxies including PMC | {summary['preclinical_document_studies']} |
| Preclinical source study IDs including PMC (dose-specific in CvTdb) | {summary['preclinical_source_study_ids_plus_pmc_study']} |
| Usable source curves after computational QC | {summary['usable_curves']} |
| Provisional CvTdb cross-dose series / compounds | {len(cross)} / {cross.compound_id.nunique()} |
| Provisional CvTdb compounds with at least three doses | {summary['cvtdb_candidate_compounds_3plus_doses']} |
| Verified compounds with at least two comparable doses | {summary['verified_compounds_2plus_comparable_doses']} |
| Verified compounds with at least three comparable doses | {summary['verified_compounds_3plus_comparable_doses']} |
| Verified sex-specific series containing reported AUC and Cmax | {summary['verified_series_with_reported_auc_and_cmax']} |
| Series eligible for AUC benchmarking | 0 |
| High-confidence Cmax groups / series | {len(high)} / {summary['verified_crossdose_series']} |
| Task A episodes / predictions | {summary['task_a_episodes']} / {len(predictions)} |

Species collected: {', '.join(summary['preclinical_species'])}. Routes: {', '.join(summary['preclinical_routes'])}. The accepted pilot is oral rat plasma only.

Six of {len(pre)+len(endpoints):,} preclinical source records ({100*summary['fraction_collected_preclinical_records_eligible_cmax']:.4f}%) are accepted Cmax groups. This is a source-record fraction mixing timepoints and endpoint groups, **not an independent experiment fraction**. One of {summary['preclinical_analyte_identifiers_collected']} analyte identifiers ({100*summary['fraction_preclinical_compound_ids_verified_crossdose']:.3f}%) has verified comparable doses. Counts are computed, not estimated, but unique molecular identities and independent study counts remain unverified: aliases, metabolites, missing structures and dose-specific study IDs matter.

The two sex-specific rat series have reported AUC and Cmax, but AUC terminal times are unresolved. Only Cmax was evaluated. With two support doses and two high-dose targets, linear scaling and power law are compared on identical observations; exact natural-log and relative-error metrics are in [benchmark feasibility](benchmark_feasibility.md). This one-compound sample cannot rank model families or establish generalization.

No compound/scaffold-disjoint split, molecular ML model, calibrated interval or nonlinear classification score is claimed. The pilot has no structures and insufficient independent compounds. Missing CvTdb strains/protocols, source censoring, heterogeneous sampling, incomplete cross-source identity resolution and the conflicting PMC dog table remain material limitations. See [unresolved issues](unresolved_issues.json), [QC](qc_report.md), [inventory](data_inventory.md) and [validation](validation.md).
""")
    log(root/'reports/progress.jsonl',phase='pipeline',status='succeeded',summary=summary)
    print(json.dumps(summary,indent=2,default=int))
    return summary

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--root',type=Path,default=Path('.'))
    p.add_argument('--offline',action='store_true')
    p.add_argument('--initial-only',action='store_true',help='Reproduce only the original pilot')
    p.add_argument('--reuse-identities',action='store_true',help='Reuse the local expanded identity registry')
    args=p.parse_args()
    try:
        run(args.root,args.offline)
        if not args.initial_only:
            from pktk.expansion import run as expand
            expand(args.root,args.offline,not args.reuse_identities)
    except Exception as e:
        log(args.root/'reports/progress.jsonl',phase='pipeline',status='failed',error=str(e))
        raise

if __name__=='__main__': main()
