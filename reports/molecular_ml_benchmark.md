# Molecular PK benchmark

These are molecular-plus-known-dose predictions of reported rat plasma PK endpoints, not dose extrapolation, not a validation of the full Dataset A, and not Dataset B transfer. All trained records come from ChEMBL. Original A/B tables and original empirical cross-dose results are unchanged.

## Executed task coverage

| task | endpoint | species | route | matrix | unit | records | original_records | compounds | families | publications |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rat_iv_plasma_clearance | clearance | rat | intravenous | plasma | mL/min/kg | 172 | 173 | 169 | 163 | 64 |
| rat_oral_plasma_cmax | cmax | rat | oral | plasma | ug/mL | 265 | 265 | 259 | 257 | 47 |

Scores use log10 endpoints (unlike the earlier Dataset B natural-log errors). Equal chemical-family weighting prevents repeated records from dominating. Standard deviations across split seeds are descriptive; test sets overlap and are not independent replications. Values in original units and family/study counts are in results/ml/metrics.csv. No target outliers were removed.

**Measured findings.** Strict scaffold+study mean log10 MAE: clearance GINE 0.337, XGBoost 0.370, descriptor Ridge 0.474; Cmax GINE 0.858, XGBoost 0.888, descriptor Ridge 0.880. These averages describe overlapping splits, not independent replications. Model rankings are exploratory; paired component comparisons below determine the limits of inference.

Random-compound Cmax: Random Forest 0.425, GINE 0.499; scaffold-only Cmax: XGBoost 0.589, GINE 0.703. These diagnostic splits share publications across partitions and do not establish unseen-study generalization. Negative R² compares squared error with the held-out test mean, not with the separately trained mean predictor. Ensemble uncertainty and empirical interval coverage are reported separately; small group counts and distribution shift prevent coverage guarantees.

**Interpretation.** The training-mean control has clearance MAE 0.326, below all molecular models on average. All models have negative mean strict-split log-space R² for both endpoints. Molecular inputs therefore do not establish useful unseen-study clearance prediction here. Of 5 prespecified GINE–XGBoost strict comparisons, 2 exploratory fixed-prediction component-bootstrap intervals exclude zero in favor of GINE; the remaining intervals include zero. These unadjusted, overlapping-split results do not establish general graph superiority, and outperforming XGBoost alone does not establish improvement over simpler controls. Cmax becomes substantially harder when both studies and scaffolds are excluded; partition populations differ, so this is not a causal estimate of a scaffold effect.

gine_ensemble: 90% nominal intervals use 12 calibration groups and 14 test groups; observed group coverage 100.0%, radius 2.149 log10 (×/÷141). Such wide intervals limit usefulness despite high observed coverage; distribution-shift exchangeability is unverified.

xgb_ensemble: 90% nominal intervals use 12 calibration groups and 14 test groups; observed group coverage 100.0%, radius 2.184 log10 (×/÷153). Such wide intervals limit usefulness despite high observed coverage; distribution-shift exchangeability is unverified.

## Models and selection

Training-log mean/median, dose-only Ridge, 17 RDKit descriptors + Ridge, 2048-bit radius-2 chirality-aware Morgan + Random Forest / XGBoost, descriptor+fingerprint MLP, GCN, GIN and GINE were evaluated. Dose is an input for every fitted predictive model; constants are deliberately unconditioned controls. XGBoost is the prespecified strong comparator, GINE the primary graph comparator. Two RF configurations, two XGBoost configurations, four Ridge penalties and two small neural configurations are selected using validation family-macro MAE only. XGBoost uses validation early stopping, maximum 600 trees; neural models use at most 100 epochs and patience 15. No refitting on calibration or test data. Exact grids and hashes are in experiment code/manifests.

## Scaffold + study excluded

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

## All partition designs

| task | split | model | mae_log10_mean | mae_log10_std | rmse_log10_mean | r2_log10_mean | mae_log10_count |
| --- | --- | --- | --- | --- | --- | --- | --- |
| rat_iv_plasma_clearance | compound | dose_ridge | 0.457 | 0.019 | 0.592 | -0.034 | 3 |
| rat_iv_plasma_clearance | compound | gcn | 0.432 | 0.034 | 0.555 | 0.076 | 3 |
| rat_iv_plasma_clearance | compound | gin | 0.454 | 0.016 | 0.598 | -0.048 | 3 |
| rat_iv_plasma_clearance | compound | gine | 0.428 | 0.032 | 0.542 | 0.112 | 3 |
| rat_iv_plasma_clearance | compound | mean | 0.454 | 0.026 | 0.591 | -0.024 | 3 |
| rat_iv_plasma_clearance | compound | median | 0.464 | 0.038 | 0.608 | -0.078 | 3 |
| rat_iv_plasma_clearance | compound | mlp | 0.367 | 0.010 | 0.524 | 0.198 | 3 |
| rat_iv_plasma_clearance | compound | rf | 0.339 | 0.036 | 0.492 | 0.292 | 3 |
| rat_iv_plasma_clearance | compound | ridge | 0.433 | 0.038 | 0.583 | -0.003 | 3 |
| rat_iv_plasma_clearance | compound | xgb | 0.417 | 0.008 | 0.566 | 0.049 | 3 |
| rat_iv_plasma_clearance | scaffold | dose_ridge | 0.499 | 0.188 | 0.599 | -0.127 | 3 |
| rat_iv_plasma_clearance | scaffold | gcn | 0.459 | 0.096 | 0.579 | -0.114 | 3 |
| rat_iv_plasma_clearance | scaffold | gin | 0.456 | 0.082 | 0.583 | -0.110 | 3 |
| rat_iv_plasma_clearance | scaffold | gine | 0.466 | 0.101 | 0.591 | -0.126 | 3 |
| rat_iv_plasma_clearance | scaffold | mean | 0.486 | 0.179 | 0.588 | -0.066 | 3 |
| rat_iv_plasma_clearance | scaffold | median | 0.465 | 0.203 | 0.578 | -0.014 | 3 |
| rat_iv_plasma_clearance | scaffold | mlp | 0.420 | 0.042 | 0.539 | 0.008 | 3 |
| rat_iv_plasma_clearance | scaffold | rf | 0.375 | 0.116 | 0.510 | 0.173 | 3 |
| rat_iv_plasma_clearance | scaffold | ridge | 0.475 | 0.157 | 0.586 | -0.101 | 3 |
| rat_iv_plasma_clearance | scaffold | xgb | 0.403 | 0.067 | 0.513 | 0.122 | 3 |
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
| rat_oral_plasma_cmax | compound | dose_ridge | 0.591 | 0.022 | 0.741 | 0.313 | 3 |
| rat_oral_plasma_cmax | compound | gcn | 0.507 | 0.028 | 0.645 | 0.471 | 3 |
| rat_oral_plasma_cmax | compound | gin | 0.507 | 0.018 | 0.644 | 0.472 | 3 |
| rat_oral_plasma_cmax | compound | gine | 0.499 | 0.016 | 0.632 | 0.494 | 3 |
| rat_oral_plasma_cmax | compound | mean | 0.683 | 0.044 | 0.903 | -0.017 | 3 |
| rat_oral_plasma_cmax | compound | median | 0.679 | 0.042 | 0.903 | -0.018 | 3 |
| rat_oral_plasma_cmax | compound | mlp | 0.434 | 0.013 | 0.577 | 0.574 | 3 |
| rat_oral_plasma_cmax | compound | rf | 0.425 | 0.071 | 0.561 | 0.586 | 3 |
| rat_oral_plasma_cmax | compound | ridge | 0.495 | 0.050 | 0.645 | 0.469 | 3 |
| rat_oral_plasma_cmax | compound | xgb | 0.387 | 0.065 | 0.520 | 0.646 | 3 |
| rat_oral_plasma_cmax | scaffold | dose_ridge | 0.660 | 0.189 | 0.851 | -0.043 | 3 |
| rat_oral_plasma_cmax | scaffold | gcn | 0.701 | 0.274 | 0.879 | -0.082 | 3 |
| rat_oral_plasma_cmax | scaffold | gin | 0.602 | 0.104 | 0.777 | 0.118 | 3 |
| rat_oral_plasma_cmax | scaffold | gine | 0.703 | 0.255 | 0.886 | -0.112 | 3 |
| rat_oral_plasma_cmax | scaffold | mean | 0.802 | 0.333 | 1.019 | -0.442 | 3 |
| rat_oral_plasma_cmax | scaffold | median | 0.775 | 0.313 | 0.983 | -0.344 | 3 |
| rat_oral_plasma_cmax | scaffold | mlp | 0.525 | 0.115 | 0.696 | 0.249 | 3 |
| rat_oral_plasma_cmax | scaffold | rf | 0.500 | 0.127 | 0.668 | 0.321 | 3 |
| rat_oral_plasma_cmax | scaffold | ridge | 0.642 | 0.067 | 0.842 | -0.106 | 3 |
| rat_oral_plasma_cmax | scaffold | xgb | 0.589 | 0.171 | 0.788 | 0.122 | 3 |
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
| rat_oral_plasma_cmax | study | dose_ridge | 0.818 | 0.492 | 0.939 | -1.543 | 2 |
| rat_oral_plasma_cmax | study | gcn | 0.782 | 0.433 | 0.906 | -1.361 | 2 |
| rat_oral_plasma_cmax | study | gin | 0.815 | 0.591 | 0.948 | -1.629 | 2 |
| rat_oral_plasma_cmax | study | gine | 0.824 | 0.575 | 0.947 | -1.619 | 2 |
| rat_oral_plasma_cmax | study | mean | 0.838 | 0.673 | 0.991 | -1.886 | 2 |
| rat_oral_plasma_cmax | study | median | 0.805 | 0.617 | 0.963 | -1.702 | 2 |
| rat_oral_plasma_cmax | study | mlp | 0.840 | 0.685 | 0.995 | -1.913 | 2 |
| rat_oral_plasma_cmax | study | rf | 0.873 | 0.716 | 1.012 | -2.048 | 2 |
| rat_oral_plasma_cmax | study | ridge | 0.830 | 0.684 | 0.996 | -1.912 | 2 |
| rat_oral_plasma_cmax | study | xgb | 0.829 | 0.651 | 0.971 | -1.777 | 2 |

Study-only clearance partitions are identical to scaffold+study and are executed once, not double-counted. Compound/scaffold-only diagnostics retain cross-partition publication overlap and cannot establish unseen-study performance. The sparse Cmax strict result has two admissible split seeds; do not interpret it as a robust ranking.

Source: [ChEMBL](https://www.ebi.ac.uk/chembl/) 37, CC BY-SA 3.0; source publication IDs and hashes remain in the local task tables. Source-derived aggregate tables/figures in these reports are supplied under [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/) with attribution to ChEMBL and the linked source publications; they are transformed analyses, not original database releases; numerical endpoint tables, predictions and model checkpoints are not committed. Code remains under its existing license.
