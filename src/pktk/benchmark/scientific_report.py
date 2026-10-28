"""Scientific reports generated from executed results; descriptive versus predictive use explicit."""
import json
import numpy as np
import pandas as pd

def md_table(frame,columns=None):
    x=frame[columns] if columns else frame
    def cell(v):
        if isinstance(v,(float,np.floating)):return '' if pd.isna(v) else f'{v:.3f}'
        return str(v).replace('|','/').replace('\n',' ')
    return '| '+' | '.join(x.columns)+' |\n| '+' | '.join(['---']*len(x.columns))+' |\n'+ '\n'.join('| '+' | '.join(cell(v) for v in row)+' |' for row in x.itertuples(index=False,name=None))

def write_reports(root,table,episodes,predictions,metrics,paired,nonlinear,feasibility):
    s=json.loads((root/'reports/expansion_summary.json').read_text())
    overall=metrics[metrics.aggregation.eq('overall')]
    main=overall[(overall.policy=='leave_one_dose_out') & (overall.n_support==2)&overall.model.isin(['linear','power'])]
    one=overall[(overall.policy=='one_extreme')&overall.model.eq('linear')]
    comparison=paired[paired.model.eq('power') & paired.n_support.eq(2)]
    coverage=pd.DataFrame(feasibility['task_coverage'])
    endpoint_counts=table.groupby('endpoint').agg(dose_groups=('source_record_id','count'),series=('endpoint_series_id','nunique'),compounds=('molecular_group_id','nunique'),studies=('study_id','nunique'),components=('component_id','nunique')).reset_index()
    model_coverage=predictions.groupby(['model','status']).size().reset_index(name='attempts')
    pool=paired[paired.model.isin(['pooled_power','shrinkage_power'])]
    near=nonlinear.groupby(['endpoint','descriptive_trend']).size().reset_index(name='series')
    examples=nonlinear.copy();examples['deviation']=abs(np.log(examples.normalized_high_low_ratio))
    strong=examples.sort_values('deviation',ascending=False).groupby('compound_name',sort=False).head(1).head(8)
    weak=examples[examples.n_doses>=3].sort_values('deviation').groupby('compound_name',sort=False).head(1).head(5)
    cols=['endpoint','task','model','n_components','n_series','n_predictions','mae_ln','rmse_ln','mean_absolute_relative_error','within_twofold']
    budget=pd.read_csv(root/'results/research/support_budget_comparisons.csv')
    ratio=metrics[(metrics.aggregation=='ratio_bin') & (metrics.policy=='one_extreme') & (metrics.task=='low_to_high') & (metrics.model=='linear')]
    examples_pred=metrics[(metrics.aggregation=='compound_name') & (metrics.policy=='leave_one_dose_out') & (metrics.n_support==2) & (metrics.task=='low_to_high') & (metrics.model=='linear')]
    findings=[]
    for endpoint in ['auc','cmax']:
        v=main[(main.endpoint==endpoint)&(main.task=='low_to_high')].set_index('model')
        c=comparison[(comparison.endpoint==endpoint)&(comparison.task=='low_to_high')].iloc[0]
        findings.append(f"- **{endpoint.upper()}, two-support low→high:** linear MAE(log) {v.loc['linear','mae_ln']:.3f}; power {v.loc['power','mae_ln']:.3f}. Paired power-minus-linear difference {c.delta_mae_ln:.3f}; exploratory component-bootstrap interval [{c.ci_low:.3f}, {c.ci_high:.3f}]. " + ('This interval includes zero; no reliable superiority is established.' if c.ci_low<=0<=c.ci_high else 'This conditional interval excludes zero, but the small selected cohort and multiple exploratory comparisons limit inference.'))
    nonlinear_ok=predictions[predictions.model.isin(['saturating','superlinear']) & predictions.status.eq('ok')]
    findings.append(f"- Nonlinear curves produce {len(nonlinear_ok)} predictions across {nonlinear_ok.component_id.nunique()} independent components. This cannot support a general nonlinear-model performance claim.")
    findings.append('- Interpolation and extrapolation are different tasks; small interpolation errors cannot justify high-dose extrapolation. Mean relative error is asymmetric and can favor a different model than log-error metrics.')
    text=f'''# Scientific analysis: limited-observation preclinical dose extrapolation

## Executive summary

Executed snapshot: **{s['primary_comparable_series']} validated multi-dose series, {s['primary_comparable_compound_identifiers']} compounds and {s['primary_independent_studies']} source studies**. Endpoint-specific evaluation uses **{len(table)} dose/endpoint observations across {len(nonlinear)} endpoint series**. These are aggregated group endpoints, not individual animals. The union contains **{feasibility['independence_components']} compound–study connected components**; that is the upper bound on independent resampling units in this collected cohort, not the number of predictions.

The benchmark executed **{len(episodes)} distinct evaluation episodes**, with **{int(predictions.status.eq('ok').sum())} successful model predictions**. All episodes pass source-record, dose-direction and reference-component leakage checks. Complex neural or molecular models were not trained. Results characterize these selected studies and do not establish performance across the chemical universe.

### Main measured findings

{chr(10).join(findings)}

## Dataset construction and provenance

{md_table(endpoint_counts)}

Collection, endpoint eligibility, compound identity and source exclusions are documented in [final dataset report](final_report.md), [inventory](data_inventory.md), [QC](qc_report.md) and [audit](composite_audit.json). {s['primary_both_endpoints_eligible_series']} series support both endpoints; AUC and Cmax tasks are evaluated separately. Species: {', '.join(s['primary_species'])}; routes: {', '.join(s['primary_routes'])}; matrices: {', '.join(s['primary_matrices'])}. Source data include different sexes, formulations, repeated-dose days and sampling schedules; each remains a separate comparable series. Shared studies/compounds remain linked for evaluation.

Only rows with endpoint-specific high-confidence flags enter the analysis. AUC definitions and exposure units must agree within each endpoint series. Values are reported group means/composite endpoints; missing data are not imputed. Whole blood is never relabeled plasma. Observed Cmax is not replaced by predicted peaks. NCA AUC0-inf includes a terminal extrapolation and does not imply a known small tail fraction. Source AUC and Cmax uncertainty/covariance is often missing or ambiguous. Lorlatinib retains CC BY-NC restrictions; generated endpoint data and figures remain local/untracked.

## Cross-dose tasks and evaluation protocols

- **One-extreme low→high:** only the lowest-dose exposure is supplied; predict every higher measured dose separately.
- **One-extreme high→low:** only the highest-dose exposure is supplied; predict every lower measured dose separately.
- **Leave-one-dose-out (LODO):** for series with at least three eligible doses, hold out each dose and fit on all remaining doses. Extreme targets are extrapolation; interior targets are interpolation. Two-dose series use the one-extreme tasks only to avoid identical duplicated episodes.
- **Ratio:** extrapolation distance is target/max-support for low→high or min-support/target for high→low, both ≥1. Interpolation has no extrapolation ratio and is stored as 1 with a distinct task label. Support span is retained separately.

{md_table(coverage)}

The same held-out dose may occur in separate support-budget tasks; those tasks are never pooled to create extra independent biological replication. All source-record memberships and reference-training memberships are saved in `results/research/`. Target exposure is used only after prediction to score error. Fits use supplied support observations only.

Pooled models are leave-compound–study-component-out: the entire connected component of the target compound/study is removed from reference training. This also excludes other compounds in the same publication. Local support for the held-out compound is explicitly permitted by the limited-observation task; this is **not zero-shot compound prediction**. Resolved molecular identities and original identifiers define grouping, with unresolved identity limitations retained.

## Baseline methods and equations

Let dose be D, exposure Y and support observations S. All fitting/scoring logarithms are natural logarithms.

1. **Linear scaling:** log Ŷ(D*) = log D* + mean_S[log Y − log D]. One fitted normalization; exponent fixed at 1. With one support this is simple proportional scaling.
2. **Power law:** log Y = a + b log D, fitted by least squares using ≥2 distinct support doses. Two support doses exactly determine a and b; there are no residual degrees of freedom for an individual prediction interval.
3. **Pooled power:** estimate each reference series' exponent, average within dependency component, then average components equally. Reference data must match endpoint, species, route, matrix, AUC definition, units and dosing regimen. At least three disjoint reference components are required. Calibrate the intercept using target-series support only. Reference slopes may differ in strain, sex, formulation and sampling method; sharing a slope across those factors is an empirical generalization hypothesis evaluated out of component, not experimental pooling.
4. **Shrinkage power:** fixed b = 0.5 b_local + 0.5 b_pooled, with ≥2 local support doses. This is a fixed empirical shrinkage rule, not a fitted hierarchical variance model; its weight is not tuned using held-out outcomes.
5. **Saturating curve:** Y = a D/(1+D/K); **superlinear curve:** Y = a D(1+D/K). Positive a,K are fitted by log-residual nonlinear least squares with ≥3 support doses. A prediction is withheld unless K lies strictly within the observed support-dose range, optimization converges and the Jacobian condition number is <1000. These empirical functions do not identify a biological mechanism. Non-identifiable attempts remain in the output, without fallback masquerading as a fitted nonlinear model.

## Quantitative prediction results

Primary comparisons below use the same two-support LODO episodes for linear and power law. First average losses over episodes within endpoint series, then over series within a dependency component, then across components equally. RMSE is the square root of the component-averaged squared log error. Relative error is |Ŷ/Y−1| (fraction, not percent); within-twofold is the corresponding weighted fraction. This prevents compounds with many sex/formulation strata dominating the overall result.

{md_table(main,cols)}

One-support linear scaling:

{md_table(one,cols)}

### Value of additional initial observations

The table matches the same extreme target and compares multiple-support models with one-extreme linear scaling. Negative delta favors the additional-support prediction. This is a matched-target descriptive comparison; it is separate from model comparisons at equal support budget.

{md_table(budget[budget.n_support.eq(2) & budget.model.isin(['linear','power'])])}

### Dose-distance dependence

One-support low→high linear performance by dose ratio (observational strata, confounded by compound/study composition):

{md_table(ratio,['endpoint','group','n_components','n_predictions','mae_ln','rmse_ln','within_twofold'])}

Performance by compound, study, species, route and dose-ratio bin is in `results/research/metrics.csv`. Scores from differently sized support budgets or differently available model cohorts are not ranked as if they were paired comparisons.

### Paired power-law minus linear error

Negative delta favors power law. Intervals are percentile bootstrap summaries from 2,000 resamples of **whole compound–study components**, after averaging paired episode differences within series/component. They are conditional exploratory summaries of these out-of-fold errors, not animal-level CIs, population guarantees, hypothesis-test p-values, or uncertainty from refitting all reference models on bootstrap samples. At least five components are required; otherwise intervals are withheld. Comparisons are multiple, exploratory and unadjusted.

{md_table(comparison,['endpoint','task','n_support','n_components','delta_mae_ln','ci_low','ci_high','leave_one_component_out_min','leave_one_component_out_max'])}

### Pooled/shrinkage comparisons on common episodes

{md_table(pool,['endpoint','policy','task','n_support','model','n_components','delta_mae_ln','ci_low','ci_high'])}

### Model availability and identifiability

{md_table(model_coverage)}

Attempt counts are not independent test cases. A sparse successful nonlinear subset cannot establish superiority over scaling. All model statuses and failure reasons are preserved in `predictions.csv`.

Sensitivity scores for known-structure-only and noncommercial-excluded evaluation cohorts are in `sensitivity_metrics.csv`. These rescore fixed predictions; they do not refit pooled references after exclusions. Leave-one-component-out paired mean-error ranges above assess whether a single biological unit changes the apparent direction of improvement.

### Where linear extrapolation succeeds or fails in this cohort

Two-support low→high errors, by compound; these examples are descriptive and are not independent confirmation of model selection.

{md_table(examples_pred.sort_values(['endpoint','mae_ln']),['endpoint','group','n_series','mae_ln','rmse_ln','within_twofold'])}

## Dose proportionality and nonlinear trends

These analyses use all observed doses **for description only**, never as inputs to held-out predictions. We report full-series power exponent, adjacent log slopes, high/low dose-normalized exposure ratio, and maximum/minimum normalized ratio. A normalized high/low ratio outside 0.8–1.25 is a descriptive flag, not an equivalence or significance test. An apparently proportional endpoint ratio can still conceal non-monotone intermediate behavior, hence adjacent slopes and max/min ratios are retained.

{md_table(near)}

Largest illustrative deviations, at most one series per compound (selected post hoc):

{md_table(strong,['compound_name','endpoint','species','sex','route','n_doses','power_exponent','normalized_high_low_ratio','normalized_max_min_ratio'])}

Examples with small endpoint-normalized deviation (also selected post hoc):

{md_table(weak,['compound_name','endpoint','species','route','power_exponent','normalized_high_low_ratio','normalized_max_min_ratio'])}

Observed sub/superproportional trends may reflect absorption, clearance, formulation, assay limits, sampling density, terminal extrapolation or study noise. They do **not** prove saturable clearance or another mechanism. For composite profiles, individual longitudinal covariance is unavailable. The published power-model framework supports explicit dose-range-dependent inference when appropriate uncertainty data exist; we do not infer statistical nonlinearity from regression fits to only two/three aggregate points. See [Smith et al., 2000](https://pubmed.ncbi.nlm.nih.gov/11145235/) and [Hummel et al., 2009](https://pubmed.ncbi.nlm.nih.gov/18386766/).

## Statistical robustness and benchmark feasibility

- **Independent units:** {feasibility['independence_components']} linked components overall; task-specific counts appear above. Study/compound counts are not summed. Sex, matrix, route and formulation strata within a component are correlated.
- **Compound-disjoint reference learning:** feasible for simple fixed pooled slopes in sufficiently populated exact strata. Evaluation retains initial measurements for the target compound. The current cohort is too small for a credible separately tuned molecular ML train/validation/test exercise.
- **Scaffolds:** {feasibility['resolved_compounds']} endpoint compounds have resolved structures. Full scaffold-group counts are in `feasibility.json`; acyclic structures share one empty-scaffold group, never artificial singleton scaffolds. Missing structures are not interpreted as acyclic. Scaffold-disjoint molecular-model comparisons were not conducted.
- **Prediction intervals/calibration:** not reported. Source aggregate SE/SD definitions, within-animal correlations and sampling uncertainty are insufficient; most power fits have zero residual degrees of freedom. Cross-component conformal calibration has too few exchangeable units in most strata. Bootstrap performance intervals are not exposure prediction intervals.
- **Nonlinearity classification:** unsupported without independent labels and uncertainty. Descriptive exponent flags are not ground truth for classifier scoring.
- **Selection sensitivity:** all primary results retain source quality flags and missing structures; noncommercial and unknown-identity subsets can be excluded explicitly. Dataset selection and variable measurement schedules limit generalization regardless of a narrow bootstrap interval.

## Important negative findings and limitations

Two-support power laws can be fitted but not meaningfully validated for curvature within their support data. Three-dose series generally cannot train and validate a nonlinear two-parameter curve with residual information simultaneously. Restricting pooled reference strata prevents incompatible exposure comparisons but often leaves too few independent reference compounds; missing predictions remain visible. High-to-low and low-to-high relative errors are asymmetric, even with the same fold error. Linear scaling can outperform a more adaptive fit when a noisy/intermediate dose distorts the local slope; interpolation performance does not establish extrapolation performance.

The available endpoints are selected public group summaries with publication/curation bias, frequently missing AUC tail fractions, source variances, full structures and individual animal data. The candidate audit still has unresolved original protocols. Compound-specific examples are exploratory and chosen after viewing results. No clinical or cross-species prediction claim is made.

## Targeted data collection after the initial analysis

The initial 57-series benchmark was executed and preserved in `research_initial_snapshot.json`. Parent-specific review of the public CEBS K04003B workbook then added two observed-Cmax series for 8:2 fluorotelomer alcohol, increasing diversity by one compound and one validated study. Male/female single-gavage groups share the specified formulation and dose volume within each stratum. Metabolites and uncertain AUC were excluded; the K04003B/C20619 alias remains unresolved, so a separate document ID is used. The complete benchmark was rerun on the updated cohort. No known blocked source was retried without a new approach. Original protocols, row selections and limitations are in `cebs_targeted.json` and `targeted_research_collection.json`.

## Figures and reproducibility

All figures are generated from executed tables as 300-dpi PNG, vector PDF and SVG in `results/figures/`:

1. `01_dataset_composition`: sources, species, compounds and studies.
2. `02_dose_exposure_auc/cmax`: dose–exposure strata.
3. `03_dose_normalized_exposure_auc/cmax`: normalized exposure patterns.
4. `04_predicted_observed`: two-support predictions versus source values.
5. `05_error_dose_ratio`: extrapolation error versus dose distance.
6. `06_performance_by_compound_name` and `07_performance_by_study_id`: correlated strata aggregated by biological unit.
7. `08_proportionality_examples`: post-hoc strong/weak observed proportionality.

Run `PYTHONPATH=src python -m pktk.benchmark.research` after collection; run the complete pipeline with `PYTHONPATH=src python -m pktk.pipeline --offline`. Config, input/output hashes and leakage-check counts are in [research_manifest.json](research_manifest.json). Fixed random seeds and plotting metadata make local replay deterministic. Source licenses apply to derived figures/data; no dataset publication is performed.

## Recommendations for future ML development

Prioritize independent compounds with ≥4 comparable doses, verified molecular structures, both observed Cmax and defined-window/NCA AUC, and animal-level uncertainty. Freeze evaluation protocols before adding flexible models. Expand endpoint-specific independent components before tuning descriptors or scaffold models. Use current linear scaling as the mandatory paired comparator; evaluate whether molecular information improves predictions beyond the same initial PK measurements. Mechanistic/hybrid models require richer concentration–time and dose-history evidence than these aggregate endpoints alone.
'''
    (root/'reports/scientific_analysis.md').write_text(text)
    summary=f'''# One-page research summary

**Question:** When can initial preclinical PK/TK measurements support dose extrapolation?

**Data:** {s['primary_comparable_series']} comparable series, {s['primary_comparable_compound_identifiers']} compounds, {s['primary_independent_studies']} studies; {s['primary_auc_eligible_series']} AUC series and {s['primary_cmax_eligible_series']} Cmax series. {len(episodes)} dose-prediction episodes represent only {feasibility['independence_components']} linked compound–study units overall. Species: {', '.join(s['primary_species'])}. Matrices and dosing protocols remain separate.

**Design:** One-support scaling and LODO interpolation/extrapolation; power laws require two support doses. Pooled exponents exclude the entire target compound–study component and match endpoint/species/route/matrix/AUC definition/regimen. No target-exposure leakage or neural models.

**Two-support low→high results** (equal component weighting; natural-log errors; relative error is a fraction):

{md_table(main[main.task=='low_to_high'],['endpoint','model','n_components','mae_ln','rmse_ln','mean_absolute_relative_error','within_twofold'])}

**Measured findings:**\n\n{chr(10).join(findings[:2])}\n\n**Interpretation:** Dose-normalized exposure varies across selected compounds. A fitted local exponent can help or amplify errors; paired comparisons across both directions and interpolation are in the full report. Descriptive nonlinear trends are not statistical proof of saturation. Curvature and prediction intervals are largely unsupported by sparse aggregate measurements.

**Independent evaluation:** Compound/study-aware bootstrap summaries are exploratory; thousands of correlated predictions do not create thousands of independent experiments. Scaffold-based molecular ML remains premature. Sampling schedules, uncertain AUC tails, missing individual variability and source-selection bias remain material limitations.

**Next priority:** Independent compounds with at least four comparable doses and well-defined AUC/Cmax, before flexible molecular models. Figures: `results/figures/`; complete methods/results: [scientific_analysis.md](scientific_analysis.md). Exact configuration and hashes: [research_manifest.json](research_manifest.json).
'''
    (root/'reports/one_page_summary.md').write_text(summary)
    (root/'reports/benchmark_feasibility.md').write_text('# Executed empirical benchmark feasibility\n\n'+md_table(coverage)+'\n\n'+
        f"{feasibility['independence_components']} linked compound-study units overall. Local scaling/power and disjoint-reference empirical models were executed. No complex molecular ML, scaffold-disjoint model score, calibrated exposure interval or nonlinear classifier is claimed. See [scientific analysis](scientific_analysis.md), [one-page summary](one_page_summary.md), and results/research/feasibility.json.\n")
