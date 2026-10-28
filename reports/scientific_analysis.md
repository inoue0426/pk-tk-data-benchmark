# Scientific analysis: limited-observation preclinical dose extrapolation

## Executive summary

Executed snapshot: **59 validated multi-dose series, 21 compounds and 19 source studies**. Endpoint-specific evaluation uses **223 dose/endpoint observations across 75 endpoint series**. These are aggregated group endpoints, not individual animals. The union contains **19 compound–study connected components**; that is the upper bound on independent resampling units in this collected cohort, not the number of predictions.

The benchmark executed **513 distinct evaluation episodes**, with **1369 successful model predictions**. All episodes pass source-record, dose-direction and reference-component leakage checks. Complex neural or molecular models were not trained. Results characterize these selected studies and do not establish performance across the chemical universe.

### Main measured findings

- **AUC, two-support low→high:** linear MAE(log) 0.525; power 0.583. Paired power-minus-linear difference 0.058; exploratory component-bootstrap interval [-0.176, 0.308]. This interval includes zero; no reliable superiority is established.
- **CMAX, two-support low→high:** linear MAE(log) 0.601; power 0.529. Paired power-minus-linear difference -0.072; exploratory component-bootstrap interval [-0.391, 0.231]. This interval includes zero; no reliable superiority is established.
- Nonlinear curves produce 3 predictions across 1 independent components. This cannot support a general nonlinear-model performance claim.
- Interpolation and extrapolation are different tasks; small interpolation errors cannot justify high-dose extrapolation. Mean relative error is asymmetric and can favor a different model than log-error metrics.

## Dataset construction and provenance

| endpoint | dose_groups | series | compounds | studies | components |
| --- | --- | --- | --- | --- | --- |
| auc | 124 | 42 | 15 | 13 | 13 |
| cmax | 99 | 33 | 12 | 12 | 12 |

Collection, endpoint eligibility, compound identity and source exclusions are documented in [final dataset report](final_report.md), [inventory](data_inventory.md), [QC](qc_report.md) and [audit](composite_audit.json). 16 series support both endpoints; AUC and Cmax tasks are evaluated separately. Species: hamster, mouse, rat; routes: intravenous, oral; matrices: blood, plasma. Source data include different sexes, formulations, repeated-dose days and sampling schedules; each remains a separate comparable series. Shared studies/compounds remain linked for evaluation.

Only rows with endpoint-specific high-confidence flags enter the analysis. AUC definitions and exposure units must agree within each endpoint series. Values are reported group means/composite endpoints; missing data are not imputed. Whole blood is never relabeled plasma. Observed Cmax is not replaced by predicted peaks. NCA AUC0-inf includes a terminal extrapolation and does not imply a known small tail fraction. Source AUC and Cmax uncertainty/covariance is often missing or ambiguous. Lorlatinib retains CC BY-NC restrictions; generated endpoint data and figures remain local/untracked.

## Cross-dose tasks and evaluation protocols

- **One-extreme low→high:** only the lowest-dose exposure is supplied; predict every higher measured dose separately.
- **One-extreme high→low:** only the highest-dose exposure is supplied; predict every lower measured dose separately.
- **Leave-one-dose-out (LODO):** for series with at least three eligible doses, hold out each dose and fit on all remaining doses. Extreme targets are extrapolation; interior targets are interpolation. Two-dose series use the one-extreme tasks only to avoid identical duplicated episodes.
- **Ratio:** extrapolation distance is target/max-support for low→high or min-support/target for high→low, both ≥1. Interpolation has no extrapolation ratio and is stored as 1 with a distinct task label. Support span is retained separately.

| endpoint | policy | task | n_support | episodes | series | compounds | studies | independence_components |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| auc | leave_one_dose_out | high_to_low | 2 | 38 | 38 | 12 | 11 | 11 |
| auc | leave_one_dose_out | high_to_low | 3 | 1 | 1 | 1 | 1 | 1 |
| auc | leave_one_dose_out | interpolation | 2 | 38 | 38 | 12 | 11 | 11 |
| auc | leave_one_dose_out | interpolation | 3 | 2 | 1 | 1 | 1 | 1 |
| auc | leave_one_dose_out | low_to_high | 2 | 38 | 38 | 12 | 11 | 11 |
| auc | leave_one_dose_out | low_to_high | 3 | 1 | 1 | 1 | 1 | 1 |
| auc | one_extreme | high_to_low | 1 | 82 | 42 | 15 | 13 | 13 |
| auc | one_extreme | low_to_high | 1 | 82 | 42 | 15 | 13 | 13 |
| cmax | leave_one_dose_out | high_to_low | 2 | 33 | 33 | 12 | 12 | 12 |
| cmax | leave_one_dose_out | interpolation | 2 | 33 | 33 | 12 | 12 | 12 |
| cmax | leave_one_dose_out | low_to_high | 2 | 33 | 33 | 12 | 12 | 12 |
| cmax | one_extreme | high_to_low | 1 | 66 | 33 | 12 | 12 | 12 |
| cmax | one_extreme | low_to_high | 1 | 66 | 33 | 12 | 12 | 12 |

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

| endpoint | task | model | n_components | n_series | n_predictions | mae_ln | rmse_ln | mean_absolute_relative_error | within_twofold |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| auc | high_to_low | linear | 11 | 38 | 38 | 0.483 | 0.648 | 0.721 | 0.674 |
| auc | high_to_low | power | 11 | 38 | 38 | 0.625 | 0.906 | 0.877 | 0.689 |
| auc | interpolation | linear | 11 | 38 | 38 | 0.292 | 0.407 | 0.342 | 0.917 |
| auc | interpolation | power | 11 | 38 | 38 | 0.291 | 0.413 | 0.348 | 0.917 |
| auc | low_to_high | linear | 11 | 38 | 38 | 0.525 | 0.700 | 0.381 | 0.701 |
| auc | low_to_high | power | 11 | 38 | 38 | 0.583 | 0.808 | 0.619 | 0.667 |
| cmax | high_to_low | linear | 12 | 33 | 33 | 0.669 | 0.916 | 1.269 | 0.646 |
| cmax | high_to_low | power | 12 | 33 | 33 | 0.519 | 0.713 | 0.559 | 0.708 |
| cmax | interpolation | linear | 12 | 33 | 33 | 0.296 | 0.380 | 0.306 | 0.990 |
| cmax | interpolation | power | 12 | 33 | 33 | 0.257 | 0.352 | 0.281 | 0.990 |
| cmax | low_to_high | linear | 12 | 33 | 33 | 0.601 | 0.793 | 0.509 | 0.615 |
| cmax | low_to_high | power | 12 | 33 | 33 | 0.529 | 0.716 | 0.601 | 0.646 |

One-support linear scaling:

| endpoint | task | model | n_components | n_series | n_predictions | mae_ln | rmse_ln | mean_absolute_relative_error | within_twofold |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| auc | high_to_low | linear | 13 | 42 | 82 | 0.521 | 0.722 | 0.922 | 0.689 |
| auc | low_to_high | linear | 13 | 42 | 82 | 0.486 | 0.684 | 0.420 | 0.716 |
| cmax | high_to_low | linear | 12 | 33 | 66 | 0.623 | 0.882 | 1.217 | 0.630 |
| cmax | low_to_high | linear | 12 | 33 | 66 | 0.703 | 0.958 | 0.606 | 0.604 |

### Value of additional initial observations

The table matches the same extreme target and compares multiple-support models with one-extreme linear scaling. Negative delta favors the additional-support prediction. This is a matched-target descriptive comparison; it is separate from model comparisons at equal support budget.

| endpoint | task | model | n_support | n_common_targets | n_components | delta_mae_ln |
| --- | --- | --- | --- | --- | --- | --- |
| auc | high_to_low | linear | 2 | 38 | 11 | -0.157 |
| auc | high_to_low | power | 2 | 38 | 11 | -0.015 |
| auc | low_to_high | linear | 2 | 38 | 11 | -0.114 |
| auc | low_to_high | power | 2 | 38 | 11 | -0.056 |
| cmax | high_to_low | linear | 2 | 33 | 12 | -0.175 |
| cmax | high_to_low | power | 2 | 33 | 12 | -0.325 |
| cmax | low_to_high | linear | 2 | 33 | 12 | -0.243 |
| cmax | low_to_high | power | 2 | 33 | 12 | -0.315 |

### Dose-distance dependence

One-support low→high linear performance by dose ratio (observational strata, confounded by compound/study composition):

| endpoint | group | n_components | n_predictions | mae_ln | rmse_ln | within_twofold |
| --- | --- | --- | --- | --- | --- | --- |
| auc | <=2 | 8 | 18 | 0.320 | 0.469 | 0.844 |
| auc | >10 | 2 | 5 | 0.673 | 0.727 | 0.250 |
| auc | >2-4 | 10 | 31 | 0.572 | 0.755 | 0.688 |
| auc | >4-10 | 9 | 28 | 0.529 | 0.732 | 0.681 |
| cmax | <=2 | 8 | 20 | 0.512 | 0.717 | 0.734 |
| cmax | >10 | 1 | 2 | 0.862 | 0.877 | 0.000 |
| cmax | >2-4 | 9 | 27 | 0.676 | 0.934 | 0.653 |
| cmax | >4-10 | 6 | 17 | 0.763 | 1.024 | 0.542 |

Performance by compound, study, species, route and dose-ratio bin is in `results/research/metrics.csv`. Scores from differently sized support budgets or differently available model cohorts are not ranked as if they were paired comparisons.

### Paired power-law minus linear error

Negative delta favors power law. Intervals are percentile bootstrap summaries from 2,000 resamples of **whole compound–study components**, after averaging paired episode differences within series/component. They are conditional exploratory summaries of these out-of-fold errors, not animal-level CIs, population guarantees, hypothesis-test p-values, or uncertainty from refitting all reference models on bootstrap samples. At least five components are required; otherwise intervals are withheld. Comparisons are multiple, exploratory and unadjusted.

| endpoint | task | n_support | n_components | delta_mae_ln | ci_low | ci_high | leave_one_component_out_min | leave_one_component_out_max |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| auc | high_to_low | 2 | 11 | 0.141 | -0.085 | 0.423 | 0.031 | 0.211 |
| auc | interpolation | 2 | 11 | -0.001 | -0.019 | 0.018 | -0.008 | 0.005 |
| auc | low_to_high | 2 | 11 | 0.058 | -0.176 | 0.308 | -0.021 | 0.126 |
| cmax | high_to_low | 2 | 12 | -0.151 | -0.558 | 0.170 | -0.219 | -0.009 |
| cmax | interpolation | 2 | 12 | -0.039 | -0.119 | 0.010 | -0.048 | -0.004 |
| cmax | low_to_high | 2 | 12 | -0.072 | -0.391 | 0.231 | -0.167 | 0.041 |

### Pooled/shrinkage comparisons on common episodes

| endpoint | policy | task | n_support | model | n_components | delta_mae_ln | ci_low | ci_high |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| auc | leave_one_dose_out | high_to_low | 2 | pooled_power | 9 | 0.014 | -0.203 | 0.238 |
| auc | leave_one_dose_out | high_to_low | 2 | shrinkage_power | 9 | -0.092 | -0.348 | 0.150 |
| auc | leave_one_dose_out | high_to_low | 3 | pooled_power | 1 | -0.384 |  |  |
| auc | leave_one_dose_out | high_to_low | 3 | shrinkage_power | 1 | -0.452 |  |  |
| auc | leave_one_dose_out | interpolation | 2 | pooled_power | 9 | 0.008 | -0.040 | 0.056 |
| auc | leave_one_dose_out | interpolation | 2 | shrinkage_power | 9 | 0.005 | -0.024 | 0.035 |
| auc | leave_one_dose_out | interpolation | 3 | pooled_power | 1 | -0.020 |  |  |
| auc | leave_one_dose_out | interpolation | 3 | shrinkage_power | 1 | -0.022 |  |  |
| auc | leave_one_dose_out | low_to_high | 2 | pooled_power | 9 | 0.007 | -0.147 | 0.180 |
| auc | leave_one_dose_out | low_to_high | 2 | shrinkage_power | 9 | -0.168 | -0.329 | 0.002 |
| auc | leave_one_dose_out | low_to_high | 3 | pooled_power | 1 | -0.423 |  |  |
| auc | leave_one_dose_out | low_to_high | 3 | shrinkage_power | 1 | -0.419 |  |  |
| auc | one_extreme | high_to_low | 1 | pooled_power | 10 | -0.054 | -0.205 | 0.107 |
| auc | one_extreme | low_to_high | 1 | pooled_power | 10 | -0.007 | -0.196 | 0.200 |
| cmax | leave_one_dose_out | high_to_low | 2 | pooled_power | 12 | 0.102 | -0.010 | 0.209 |
| cmax | leave_one_dose_out | high_to_low | 2 | shrinkage_power | 12 | -0.145 | -0.341 | 0.039 |
| cmax | leave_one_dose_out | interpolation | 2 | pooled_power | 12 | 0.018 | -0.004 | 0.042 |
| cmax | leave_one_dose_out | interpolation | 2 | shrinkage_power | 12 | -0.011 | -0.056 | 0.018 |
| cmax | leave_one_dose_out | low_to_high | 2 | pooled_power | 12 | 0.150 | 0.039 | 0.261 |
| cmax | leave_one_dose_out | low_to_high | 2 | shrinkage_power | 12 | -0.139 | -0.384 | 0.054 |
| cmax | one_extreme | high_to_low | 1 | pooled_power | 12 | 0.136 | 0.036 | 0.239 |
| cmax | one_extreme | low_to_high | 1 | pooled_power | 12 | 0.098 | -0.004 | 0.196 |

### Model availability and identifiability

| model | status | attempts |
| --- | --- | --- |
| linear | ok | 513 |
| pooled_power | not_identifiable | 67 |
| pooled_power | ok | 446 |
| power | not_identifiable | 296 |
| power | ok | 217 |
| saturating | not_identifiable | 513 |
| shrinkage_power | not_identifiable | 323 |
| shrinkage_power | ok | 190 |
| superlinear | not_identifiable | 510 |
| superlinear | ok | 3 |

Attempt counts are not independent test cases. A sparse successful nonlinear subset cannot establish superiority over scaling. All model statuses and failure reasons are preserved in `predictions.csv`.

Sensitivity scores for known-structure-only and noncommercial-excluded evaluation cohorts are in `sensitivity_metrics.csv`. These rescore fixed predictions; they do not refit pooled references after exclusions. Leave-one-component-out paired mean-error ranges above assess whether a single biological unit changes the apparent direction of improvement.

### Where linear extrapolation succeeds or fails in this cohort

Two-support low→high errors, by compound; these examples are descriptive and are not independent confirmation of model selection.

| endpoint | group | n_series | mae_ln | rmse_ln | within_twofold |
| --- | --- | --- | --- | --- | --- |
| auc | methyl eugenol | 1 | 0.043 | 0.043 | 1.000 |
| auc | Naphthalene | 4 | 0.200 | 0.214 | 1.000 |
| auc | paeonol | 1 | 0.221 | 0.221 | 1.000 |
| auc | Benzophenone | 4 | 0.266 | 0.288 | 1.000 |
| auc | lumefantrine | 1 | 0.351 | 0.351 | 1.000 |
| auc | Formamide | 8 | 0.372 | 0.423 | 0.875 |
| auc | lorlatinib | 1 | 0.389 | 0.389 | 1.000 |
| auc | DA-6886 | 2 | 0.604 | 0.732 | 0.500 |
| auc | imperatorin | 1 | 0.746 | 0.746 | 0.000 |
| auc | 2,4-Dichlorophenoxyacetic acid | 3 | 0.808 | 0.854 | 0.333 |
| auc | Bromodichloromethane | 8 | 1.131 | 1.207 | 0.250 |
| auc | 2-Hydroxy-4-methoxybenzophenone | 4 | 1.133 | 1.332 | 0.250 |
| cmax | methyl eugenol | 2 | 0.027 | 0.035 | 1.000 |
| cmax | lorlatinib | 1 | 0.095 | 0.095 | 1.000 |
| cmax | 2-Methyltetrahydrofuran | 4 | 0.278 | 0.333 | 1.000 |
| cmax | lumefantrine | 1 | 0.383 | 0.383 | 1.000 |
| cmax | Anthraquinone | 4 | 0.462 | 0.613 | 0.500 |
| cmax | Perfluorooctanoic Acid | 2 | 0.481 | 0.545 | 0.500 |
| cmax | Isoeugenol | 2 | 0.513 | 0.517 | 1.000 |
| cmax | Fluorotelomer alcohol 8+2 | 2 | 0.527 | 0.582 | 0.500 |
| cmax | Bromodichloromethane | 8 | 0.577 | 0.736 | 0.625 |
| cmax | 2-Hydroxy-4-methoxybenzophenone | 4 | 1.028 | 1.181 | 0.250 |
| cmax | DA-6886 | 1 | 1.371 | 1.371 | 0.000 |
| cmax | HY-071085 | 2 | 1.474 | 1.479 | 0.000 |

## Dose proportionality and nonlinear trends

These analyses use all observed doses **for description only**, never as inputs to held-out predictions. We report full-series power exponent, adjacent log slopes, high/low dose-normalized exposure ratio, and maximum/minimum normalized ratio. A normalized high/low ratio outside 0.8–1.25 is a descriptive flag, not an equivalence or significance test. An apparently proportional endpoint ratio can still conceal non-monotone intermediate behavior, hence adjacent slopes and max/min ratios are retained.

| endpoint | descriptive_trend | series |
| --- | --- | --- |
| auc | near_proportional | 13 |
| auc | subproportional | 4 |
| auc | superproportional | 25 |
| cmax | near_proportional | 8 |
| cmax | subproportional | 10 |
| cmax | superproportional | 15 |

Largest illustrative deviations, at most one series per compound (selected post hoc):

| compound_name | endpoint | species | sex | route | n_doses | power_exponent | normalized_high_low_ratio | normalized_max_min_ratio |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Bromodichloromethane | auc | mouse | female | oral | 3 | 3.063 | 17.468 | 17.468 |
| 2-Hydroxy-4-methoxybenzophenone | auc | mouse | female | oral | 3 | 2.261 | 9.833 | 9.833 |
| HY-071085 | cmax | rat | female | oral | 3 | 2.610 | 9.322 | 9.322 |
| DA-6886 | cmax | rat | male | oral | 3 | 1.932 | 8.471 | 8.471 |
| Formamide | auc | mouse | male | intravenous | 3 | 0.366 | 0.249 | 4.023 |
| 2-Methyltetrahydrofuran | cmax | mouse | male | oral | 3 | 1.350 | 2.784 | 2.804 |
| Anthraquinone | cmax | mouse | female | oral | 3 | 0.552 | 0.370 | 2.700 |
| lumefantrine | cmax | rat | male | oral | 3 | 0.308 | 0.383 | 3.170 |

Examples with small endpoint-normalized deviation (also selected post hoc):

| compound_name | endpoint | species | route | power_exponent | normalized_high_low_ratio | normalized_max_min_ratio |
| --- | --- | --- | --- | --- | --- | --- |
| methyl eugenol | cmax | mouse | oral | 0.992 | 0.994 | 1.020 |
| Anthraquinone | cmax | rat | oral | 1.001 | 1.012 | 1.077 |
| 2-Methyltetrahydrofuran | cmax | rat | oral | 0.980 | 0.962 | 1.275 |
| Naphthalene | auc | mouse | intravenous | 1.021 | 1.062 | 1.466 |
| Benzophenone | auc | mouse | oral | 1.054 | 1.074 | 1.376 |

Observed sub/superproportional trends may reflect absorption, clearance, formulation, assay limits, sampling density, terminal extrapolation or study noise. They do **not** prove saturable clearance or another mechanism. For composite profiles, individual longitudinal covariance is unavailable. The published power-model framework supports explicit dose-range-dependent inference when appropriate uncertainty data exist; we do not infer statistical nonlinearity from regression fits to only two/three aggregate points. See [Smith et al., 2000](https://pubmed.ncbi.nlm.nih.gov/11145235/) and [Hummel et al., 2009](https://pubmed.ncbi.nlm.nih.gov/18386766/).

## Statistical robustness and benchmark feasibility

- **Independent units:** 19 linked components overall; task-specific counts appear above. Study/compound counts are not summed. Sex, matrix, route and formulation strata within a component are correlated.
- **Compound-disjoint reference learning:** feasible for simple fixed pooled slopes in sufficiently populated exact strata. Evaluation retains initial measurements for the target compound. The current cohort is too small for a credible separately tuned molecular ML train/validation/test exercise.
- **Scaffolds:** 19 endpoint compounds have resolved structures. Full scaffold-group counts are in `feasibility.json`; acyclic structures share one empty-scaffold group, never artificial singleton scaffolds. Missing structures are not interpreted as acyclic. Scaffold-disjoint molecular-model comparisons were not conducted.
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
