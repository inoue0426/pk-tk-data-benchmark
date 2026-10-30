# Generalization analysis

These are molecular-plus-known-dose predictions of reported rat plasma PK endpoints, not dose extrapolation, not a validation of the full Dataset A, and not Dataset B transfer. All trained records come from ChEMBL. Original A/B tables and original empirical cross-dose results are unchanged.

**Measured findings.** Strict scaffold+study mean log10 MAE: clearance GINE 0.337, XGBoost 0.370, descriptor Ridge 0.474; Cmax GINE 0.858, XGBoost 0.888, descriptor Ridge 0.880. These averages describe overlapping splits, not independent replications. Model rankings are exploratory; paired component comparisons below determine the limits of inference.

Random-compound Cmax: Random Forest 0.425, GINE 0.499; scaffold-only Cmax: XGBoost 0.589, GINE 0.703. These diagnostic splits share publications across partitions and do not establish unseen-study generalization. Negative R² compares squared error with the held-out test mean, not with the separately trained mean predictor. Ensemble uncertainty and empirical interval coverage are reported separately; small group counts and distribution shift prevent coverage guarantees.

**Interpretation.** The training-mean control has clearance MAE 0.326, below all molecular models on average. All models have negative mean strict-split log-space R² for both endpoints. Molecular inputs therefore do not establish useful unseen-study clearance prediction here. Of 5 prespecified GINE–XGBoost strict comparisons, 2 exploratory fixed-prediction component-bootstrap intervals exclude zero in favor of GINE; the remaining intervals include zero. These unadjusted, overlapping-split results do not establish general graph superiority, and outperforming XGBoost alone does not establish improvement over simpler controls. Cmax becomes substantially harder when both studies and scaffolds are excluded; partition populations differ, so this is not a causal estimate of a scaffold effect.

gine_ensemble: 90% nominal intervals use 12 calibration groups and 14 test groups; observed group coverage 100.0%, radius 2.149 log10 (×/÷141). Such wide intervals limit usefulness despite high observed coverage; distribution-shift exchangeability is unverified.

xgb_ensemble: 90% nominal intervals use 12 calibration groups and 14 test groups; observed group coverage 100.0%, radius 2.184 log10 (×/÷153). Such wide intervals limit usefulness despite high observed coverage; distribution-shift exchangeability is unverified.

## Split constraints

Compound families link full connectivity and canonical neutral tautomer identities irrespective of stereochemistry. Study components link original publications transitively, and scaffold components use stereo-free Bemis–Murcko scaffolds (all acyclic compounds grouped). No target values determine partitions. All scales, imputers and calibration functions are fit within permitted partitions.

| task | split | seed | eligible | partition_compounds | n_groups | largest_group_records | studies_crossing_partitions | scaffolds_crossing_partitions |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rat_iv_plasma_clearance | compound | 42 | True | {'train': 89, 'validation': 23, 'calibration': 26, 'test': 25} | 163 | 3 | 26 | 16 |
| rat_iv_plasma_clearance | compound | 137 | True | {'train': 89, 'validation': 24, 'calibration': 26, 'test': 24} | 163 | 3 | 26 | 13 |
| rat_iv_plasma_clearance | compound | 2026 | True | {'train': 87, 'validation': 25, 'calibration': 25, 'test': 26} | 163 | 3 | 28 | 15 |
| rat_iv_plasma_clearance | scaffold | 42 | True | {'train': 96, 'validation': 25, 'calibration': 21, 'test': 21} | 104 | 10 | 22 | 0 |
| rat_iv_plasma_clearance | scaffold | 137 | True | {'train': 92, 'validation': 27, 'calibration': 23, 'test': 21} | 104 | 10 | 19 | 0 |
| rat_iv_plasma_clearance | scaffold | 2026 | True | {'train': 88, 'validation': 26, 'calibration': 25, 'test': 24} | 104 | 10 | 20 | 0 |
| rat_iv_plasma_clearance | study | 42 | True | {'train': 92, 'validation': 24, 'calibration': 24, 'test': 23} | 63 | 14 | 0 | 0 |
| rat_iv_plasma_clearance | study | 137 | True | {'train': 90, 'validation': 24, 'calibration': 26, 'test': 23} | 63 | 14 | 0 | 0 |
| rat_iv_plasma_clearance | study | 2026 | True | {'train': 89, 'validation': 25, 'calibration': 21, 'test': 28} | 63 | 14 | 0 | 0 |
| rat_iv_plasma_clearance | scaffold_study | 42 | True | {'train': 92, 'validation': 24, 'calibration': 24, 'test': 23} | 63 | 14 | 0 | 0 |
| rat_iv_plasma_clearance | scaffold_study | 137 | True | {'train': 90, 'validation': 24, 'calibration': 26, 'test': 23} | 63 | 14 | 0 | 0 |
| rat_iv_plasma_clearance | scaffold_study | 2026 | True | {'train': 89, 'validation': 25, 'calibration': 21, 'test': 28} | 63 | 14 | 0 | 0 |
| rat_oral_plasma_cmax | compound | 42 | True | {'train': 141, 'validation': 40, 'calibration': 39, 'test': 37} | 257 | 2 | 25 | 32 |
| rat_oral_plasma_cmax | compound | 137 | True | {'train': 141, 'validation': 40, 'calibration': 39, 'test': 37} | 257 | 2 | 28 | 29 |
| rat_oral_plasma_cmax | compound | 2026 | True | {'train': 141, 'validation': 39, 'calibration': 40, 'test': 37} | 257 | 2 | 23 | 33 |
| rat_oral_plasma_cmax | scaffold | 42 | True | {'train': 142, 'validation': 36, 'calibration': 41, 'test': 38} | 119 | 31 | 17 | 0 |
| rat_oral_plasma_cmax | scaffold | 137 | True | {'train': 139, 'validation': 37, 'calibration': 35, 'test': 46} | 119 | 31 | 20 | 0 |
| rat_oral_plasma_cmax | scaffold | 2026 | True | {'train': 159, 'validation': 33, 'calibration': 33, 'test': 32} | 119 | 31 | 21 | 0 |
| rat_oral_plasma_cmax | study | 42 | True | {'train': 136, 'validation': 31, 'calibration': 62, 'test': 28} | 46 | 107 | 0 | 0 |
| rat_oral_plasma_cmax | study | 137 | True | {'train': 142, 'validation': 37, 'calibration': 36, 'test': 42} | 46 | 107 | 0 | 1 |
| rat_oral_plasma_cmax | study | 2026 | False | {'train': 113, 'validation': 14, 'calibration': 118, 'test': 12} | 46 | 107 | 0 | 0 |
| rat_oral_plasma_cmax | scaffold_study | 42 | True | {'train': 136, 'validation': 31, 'calibration': 62, 'test': 28} | 45 | 107 | 0 | 0 |
| rat_oral_plasma_cmax | scaffold_study | 137 | True | {'train': 139, 'validation': 34, 'calibration': 32, 'test': 52} | 45 | 107 | 0 | 0 |
| rat_oral_plasma_cmax | scaffold_study | 2026 | False | {'train': 217, 'validation': 14, 'calibration': 12, 'test': 14} | 45 | 107 | 0 | 0 |

2 requested task/split seeds fail the minimum partition-size rule and were not trained. They were not replaced with more favorable random seeds. Random-compound and scaffold-only scores are diagnostics with explicit publication overlap, not leakage-free study generalization claims.

## Generalization contrasts

| task | split | model | mae_log10_mean | mae_log10_std | rmse_log10_mean | r2_log10_mean | mae_log10_count |
| --- | --- | --- | --- | --- | --- | --- | --- |
| rat_iv_plasma_clearance | compound | gine | 0.428 | 0.032 | 0.542 | 0.112 | 3 |
| rat_iv_plasma_clearance | compound | rf | 0.339 | 0.036 | 0.492 | 0.292 | 3 |
| rat_iv_plasma_clearance | compound | xgb | 0.417 | 0.008 | 0.566 | 0.049 | 3 |
| rat_iv_plasma_clearance | scaffold | gine | 0.466 | 0.101 | 0.591 | -0.126 | 3 |
| rat_iv_plasma_clearance | scaffold | rf | 0.375 | 0.116 | 0.510 | 0.173 | 3 |
| rat_iv_plasma_clearance | scaffold | xgb | 0.403 | 0.067 | 0.513 | 0.122 | 3 |
| rat_iv_plasma_clearance | scaffold_study | gine | 0.337 | 0.035 | 0.424 | -0.159 | 3 |
| rat_iv_plasma_clearance | scaffold_study | rf | 0.403 | 0.039 | 0.487 | -0.546 | 3 |
| rat_iv_plasma_clearance | scaffold_study | xgb | 0.370 | 0.080 | 0.442 | -0.297 | 3 |
| rat_oral_plasma_cmax | compound | gine | 0.499 | 0.016 | 0.632 | 0.494 | 3 |
| rat_oral_plasma_cmax | compound | rf | 0.425 | 0.071 | 0.561 | 0.586 | 3 |
| rat_oral_plasma_cmax | compound | xgb | 0.387 | 0.065 | 0.520 | 0.646 | 3 |
| rat_oral_plasma_cmax | scaffold | gine | 0.703 | 0.255 | 0.886 | -0.112 | 3 |
| rat_oral_plasma_cmax | scaffold | rf | 0.500 | 0.127 | 0.668 | 0.321 | 3 |
| rat_oral_plasma_cmax | scaffold | xgb | 0.589 | 0.171 | 0.788 | 0.122 | 3 |
| rat_oral_plasma_cmax | scaffold_study | gine | 0.858 | 0.263 | 1.040 | -0.387 | 2 |
| rat_oral_plasma_cmax | scaffold_study | rf | 0.831 | 0.179 | 1.021 | -0.305 | 2 |
| rat_oral_plasma_cmax | scaffold_study | xgb | 0.888 | 0.221 | 1.091 | -0.494 | 2 |
| rat_oral_plasma_cmax | study | gine | 0.824 | 0.575 | 0.947 | -1.619 | 2 |
| rat_oral_plasma_cmax | study | rf | 0.873 | 0.716 | 1.012 | -2.048 | 2 |
| rat_oral_plasma_cmax | study | xgb | 0.829 | 0.651 | 0.971 | -1.777 | 2 |

There is one eligible source (ChEMBL) and one species (rat). Source- and species-stratified scores therefore equal each task score; genuine cross-source or cross-species tests are unsupported. Differences across partition schemes mix changes in chemical distance, study composition and test difficulty; they do not isolate a causal scaffold effect.

## Learning curves

| model | train_fraction | n_train_families | mae_log10 | rmse_log10 | validation_mae_log10 |
| --- | --- | --- | --- | --- | --- |
| gine | 1.000 | 92 | 0.304 | 0.350 | 0.334 |
| gine | 0.350 | 41 | 0.295 | 0.351 | 0.347 |
| gine | 0.650 | 76 | 0.298 | 0.348 | 0.339 |
| xgb | 1.000 | 92 | 0.294 | 0.347 | 0.346 |
| xgb | 0.350 | 41 | 0.329 | 0.412 | 0.301 |
| xgb | 0.650 | 76 | 0.294 | 0.350 | 0.344 |

Learning curves use nested subsets of whole training study/scaffold groups, with the same held-out validation/calibration/test sets. Each fraction is retuned on validation with the same declared grid. One split and two smaller fractions are exploratory; curves need not be monotone.

## Applicability domain and transfer

Error histograms and the five largest absolute errors per strict model/split are saved in results/ml/figures/08_error_distributions and results/ml/outliers.csv. These post-hoc diagnostic selections never change task eligibility or refit models. Nearest-training Tanimoto similarity is calculated with cached chirality-aware Morgan fingerprints. Family-averaged error correlations and risk curves are descriptive, not proof of a universal similarity threshold. Source metadata features are not tested because source is constant and publication IDs are unavailable as meaningful features for a new study.

Dataset A→B transfer was not executed. The current molecular tasks are conditional CL/Cmax, not AUC dose-response representations; Dataset B has only 19 dependency components overall and 11/12 units in the established two-support comparison. A valid transfer test would need identity-disjoint representation training, disjoint model selection and calibration, plus enough target components to compare against linear scaling. This study does not establish that feasibility, and no transfer benefit is claimed.

Source: [ChEMBL](https://www.ebi.ac.uk/chembl/) 37, CC BY-SA 3.0; source publication IDs and hashes remain in the local task tables. Source-derived aggregate tables/figures in these reports are supplied under [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/) with attribution to ChEMBL and the linked source publications; they are transformed analyses, not original database releases; numerical endpoint tables, predictions and model checkpoints are not committed. Code remains under its existing license.
