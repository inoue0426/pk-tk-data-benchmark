# Final molecular PK research summary

These are molecular-plus-known-dose predictions of reported rat plasma PK endpoints, not dose extrapolation, not a validation of the full Dataset A, and not Dataset B transfer. All trained records come from ChEMBL. Original A/B tables and original empirical cross-dose results are unchanged.

**Measured findings.** Strict scaffold+study mean log10 MAE: clearance GINE 0.337, XGBoost 0.370, descriptor Ridge 0.474; Cmax GINE 0.858, XGBoost 0.888, descriptor Ridge 0.880. These averages describe overlapping splits, not independent replications. Model rankings are exploratory; paired component comparisons below determine the limits of inference.

Random-compound Cmax: Random Forest 0.425, GINE 0.499; scaffold-only Cmax: XGBoost 0.589, GINE 0.703. These diagnostic splits share publications across partitions and do not establish unseen-study generalization. Negative R² compares squared error with the held-out test mean, not with the separately trained mean predictor. Ensemble uncertainty and empirical interval coverage are reported separately; small group counts and distribution shift prevent coverage guarantees.

**Interpretation.** The training-mean control has clearance MAE 0.326, below all molecular models on average. All models have negative mean strict-split log-space R² for both endpoints. Molecular inputs therefore do not establish useful unseen-study clearance prediction here. Of 5 prespecified GINE–XGBoost strict comparisons, 2 exploratory fixed-prediction component-bootstrap intervals exclude zero in favor of GINE; the remaining intervals include zero. These unadjusted, overlapping-split results do not establish general graph superiority, and outperforming XGBoost alone does not establish improvement over simpler controls. Cmax becomes substantially harder when both studies and scaffolds are excluded; partition populations differ, so this is not a causal estimate of a scaffold effect.

gine_ensemble: 90% nominal intervals use 12 calibration groups and 14 test groups; observed group coverage 100.0%, radius 2.149 log10 (×/÷141). Such wide intervals limit usefulness despite high observed coverage; distribution-shift exchangeability is unverified.

xgb_ensemble: 90% nominal intervals use 12 calibration groups and 14 test groups; observed group coverage 100.0%, radius 2.184 log10 (×/÷153). Such wide intervals limit usefulness despite high observed coverage; distribution-shift exchangeability is unverified.

## Actual scope and research questions

The conservative ML layer uses **413 unique structure-resolved compounds**, 438 original eligible endpoint records and 437 assay-specific publication/identity/dose cells. It does not train on all 2,713 Dataset A compounds. 190 primary model runs and 9 focused follow-up runs completed. Actual first/last experiment events: 2026-10-30T14:41:42.898766+00:00 / 2026-10-30T14:43:38.472901+00:00; summed fit/tuning wall times 117.1 seconds (overlapping processes make this different from elapsed session time). This was not 72 hours of execution.

| task | endpoint | species | route | matrix | unit | records | original_records | compounds | families | publications |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rat_iv_plasma_clearance | clearance | rat | intravenous | plasma | mL/min/kg | 172 | 173 | 169 | 163 | 64 |
| rat_oral_plasma_cmax | cmax | rat | oral | plasma | ug/mL | 265 | 265 | 259 | 257 | 47 |

RQ1 is assessed only for conditional rat plasma clearance/Cmax. AUC, half-life and volume are withheld for incompatible/insufficient definitions. RQ2 compares small graph networks with prespecified XGBoost and other baselines. RQ3 includes compound, scaffold and study-aware designs but cannot assess genuine cross-source generalization. RQ4 has a limited independent-calibration assessment. RQ5 transfer is unsupported and was not forced.

## Strict scaffold + study results

| task | split | model | mae_log10_mean | mae_log10_std | rmse_log10_mean | r2_log10_mean | mae_log10_count |
| --- | --- | --- | --- | --- | --- | --- | --- |
| rat_iv_plasma_clearance | scaffold_study | dose_ridge | 0.331 | 0.043 | 0.411 | -0.085 | 3 |
| rat_iv_plasma_clearance | scaffold_study | gcn | 0.339 | 0.031 | 0.409 | -0.078 | 3 |
| rat_iv_plasma_clearance | scaffold_study | gin | 0.335 | 0.035 | 0.413 | -0.098 | 3 |
| rat_iv_plasma_clearance | scaffold_study | gine | 0.337 | 0.035 | 0.424 | -0.159 | 3 |
| rat_iv_plasma_clearance | scaffold_study | mean | 0.326 | 0.033 | 0.398 | -0.021 | 3 |
| rat_iv_plasma_clearance | scaffold_study | median | 0.349 | 0.060 | 0.439 | -0.237 | 3 |
| rat_iv_plasma_clearance | scaffold_study | mlp | 0.331 | 0.022 | 0.402 | -0.043 | 3 |
| rat_iv_plasma_clearance | scaffold_study | rf | 0.403 | 0.039 | 0.487 | -0.546 | 3 |
| rat_iv_plasma_clearance | scaffold_study | ridge | 0.474 | 0.204 | 0.565 | -1.606 | 3 |
| rat_iv_plasma_clearance | scaffold_study | xgb | 0.370 | 0.080 | 0.442 | -0.297 | 3 |
| rat_oral_plasma_cmax | scaffold_study | dose_ridge | 0.844 | 0.281 | 1.022 | -0.348 | 2 |
| rat_oral_plasma_cmax | scaffold_study | gcn | 0.813 | 0.202 | 0.995 | -0.253 | 2 |
| rat_oral_plasma_cmax | scaffold_study | gin | 0.833 | 0.197 | 1.012 | -0.294 | 2 |
| rat_oral_plasma_cmax | scaffold_study | gine | 0.858 | 0.263 | 1.040 | -0.387 | 2 |
| rat_oral_plasma_cmax | scaffold_study | mean | 0.951 | 0.281 | 1.155 | -0.697 | 2 |
| rat_oral_plasma_cmax | scaffold_study | median | 0.941 | 0.258 | 1.151 | -0.674 | 2 |
| rat_oral_plasma_cmax | scaffold_study | mlp | 0.931 | 0.308 | 1.143 | -0.670 | 2 |
| rat_oral_plasma_cmax | scaffold_study | rf | 0.831 | 0.179 | 1.021 | -0.305 | 2 |
| rat_oral_plasma_cmax | scaffold_study | ridge | 0.880 | 0.243 | 1.072 | -0.462 | 2 |
| rat_oral_plasma_cmax | scaffold_study | xgb | 0.888 | 0.221 | 1.091 | -0.494 | 2 |

Scores use log10 endpoints (unlike the earlier Dataset B natural-log errors). Equal chemical-family weighting prevents repeated records from dominating. Standard deviations across split seeds are descriptive; test sets overlap and are not independent replications. Values in original units and family/study counts are in results/ml/metrics.csv. No target outliers were removed.

## Evidence for or against graph improvement

| task | split | seed | model | n_test_families | n_dependency_groups | delta_family_mae_log10 | delta_component_mae_log10 | ci_low | ci_high | interval_status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rat_iv_plasma_clearance | scaffold_study | 42 | gine | 23 | 14 | 0.009 | 0.006 | -0.021 | 0.032 | exploratory_component_bootstrap_fixed_predictions |
| rat_iv_plasma_clearance | scaffold_study | 137 | gine | 23 | 16 | -0.120 | -0.147 | -0.267 | -0.032 | exploratory_component_bootstrap_fixed_predictions |
| rat_iv_plasma_clearance | scaffold_study | 2026 | gine | 28 | 9 | 0.010 | 0.005 | -0.012 | 0.021 | exploratory_component_bootstrap_fixed_predictions |
| rat_oral_plasma_cmax | scaffold_study | 42 | gine | 28 | 11 | -0.060 | -0.028 | -0.118 | 0.073 | exploratory_component_bootstrap_fixed_predictions |
| rat_oral_plasma_cmax | scaffold_study | 137 | gine | 52 | 9 | -0.000 | -0.083 | -0.134 | -0.031 | exploratory_component_bootstrap_fixed_predictions |

The prespecified GINE comparison, rather than selecting a winning architecture after seeing test outcomes, is the basis for graph-versus-fingerprint interpretation. A negative difference favors GINE. Two-seed Cmax strict evaluation is particularly underpowered. No positive result from a single split is treated as established improvement. Full model rankings are descriptive, not selection for deployment.

## Limitations and recommendations

The original sources contain far more endpoints than the scientifically compatible ML cohorts. Unknown original matrix/study metadata, unspecified AUC windows and half-life phases prevent training several requested targets. Retained source-reported assays still vary in sampling schedule, formulation, sex and strain; study exclusion reduces leakage but cannot remove all label heterogeneity. Unspecified stereocenters are preserved and grouped conservatively. Single-source rat cohorts do not establish broader species/source generalization. Small held-out study groups and overlapping split seeds limit precision. No calibrated-interval guarantee under scaffold shift, mechanistic saturation inference or Dataset B transfer claim is made.

Prefer additional independent, metadata-complete studies and prospective validation over larger networks. Keep the known dose input and a dose-only control when measuring the incremental value of molecular features. Investigate failures with profile-level and experimental variance information.

## Deliverables and reproducibility

See [dataset audit](ml_dataset_audit.md), [full benchmark](molecular_ml_benchmark.md), [GNN experiments](gnn_experiments.md), [generalization](generalization_analysis.md), [uncertainty](uncertainty_analysis.md) and [execution log](experiment_log.md). Exact split/task tables, local model checkpoints, prediction rows and metrics are under data/processed/ml and results/ml. Figures are vector PDF plus 300-dpi PNG under results/ml/figures. The reproducible launcher uses no paid hardware; caches permit offline experiment replay.

Source: [ChEMBL](https://www.ebi.ac.uk/chembl/) 37, CC BY-SA 3.0; source publication IDs and hashes remain in the local task tables. Source-derived aggregate tables/figures in these reports are supplied under [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/) with attribution to ChEMBL and the linked source publications; they are transformed analyses, not original database releases; numerical endpoint tables, predictions and model checkpoints are not committed. Code remains under its existing license.
