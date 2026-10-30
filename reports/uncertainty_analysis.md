# Uncertainty and reliability

These are molecular-plus-known-dose predictions of reported rat plasma PK endpoints, not dose extrapolation, not a validation of the full Dataset A, and not Dataset B transfer. All trained records come from ChEMBL. Original A/B tables and original empirical cross-dose results are unchanged.

**Measured findings.** Strict scaffold+study mean log10 MAE: clearance GINE 0.337, XGBoost 0.370, descriptor Ridge 0.474; Cmax GINE 0.858, XGBoost 0.888, descriptor Ridge 0.880. These averages describe overlapping splits, not independent replications. Model rankings are exploratory; paired component comparisons below determine the limits of inference.

Random-compound Cmax: Random Forest 0.425, GINE 0.499; scaffold-only Cmax: XGBoost 0.589, GINE 0.703. These diagnostic splits share publications across partitions and do not establish unseen-study generalization. Negative R² compares squared error with the held-out test mean, not with the separately trained mean predictor. Ensemble uncertainty and empirical interval coverage are reported separately; small group counts and distribution shift prevent coverage guarantees.

**Interpretation.** The training-mean control has clearance MAE 0.326, below all molecular models on average. All models have negative mean strict-split log-space R² for both endpoints. Molecular inputs therefore do not establish useful unseen-study clearance prediction here. Of 5 prespecified GINE–XGBoost strict comparisons, 2 exploratory fixed-prediction component-bootstrap intervals exclude zero in favor of GINE; the remaining intervals include zero. These unadjusted, overlapping-split results do not establish general graph superiority, and outperforming XGBoost alone does not establish improvement over simpler controls. Cmax becomes substantially harder when both studies and scaffolds are excluded; partition populations differ, so this is not a causal estimate of a scaffold effect.

gine_ensemble: 90% nominal intervals use 12 calibration groups and 14 test groups; observed group coverage 100.0%, radius 2.149 log10 (×/÷141). Such wide intervals limit usefulness despite high observed coverage; distribution-shift exchangeability is unverified.

xgb_ensemble: 90% nominal intervals use 12 calibration groups and 14 test groups; observed group coverage 100.0%, radius 2.184 log10 (×/÷153). Such wide intervals limit usefulness despite high observed coverage; distribution-shift exchangeability is unverified.

## Calibration protocol

Finite-sample construction follows [Angelopoulos and Bates](https://arxiv.org/abs/2107.07511), applied here to group-maximum residuals; exchangeability assumptions remain essential.

A separate calibration partition was reserved before any fitting; it is never used for hyperparameters or early stopping. For study-excluded splits, nonconformity scores are the maximum absolute log10 residual within each calibration dependency group. At nominal coverage 1−α, radius q is the ceil((n+1)(1−α))-th ordered score. If that rank exceeds n, no finite interval is reported. Intervals are symmetric in log10 space and multiplicative in original units. Coverage is evaluated for **all observations in each held-out group**, alongside family coverage. This targets group coverage, not animal-level variability.

Under exchangeable calibration/test groups this is the usual split-conformal construction. Scaffold/source-like distribution shifts, unequal study designs and purposive public-data selection do not guarantee exchangeability. Accordingly, these results are empirical coverage assessments, not calibrated guarantees on unseen scaffolds. Compound/scaffold-only splits with shared publications are excluded from this group-calibration claim.

## Three-member ensembles (fixed clearance scaffold/study split, seed 42)

| task | split | seed | model | n_members | n_test_families | mae_log10 | uncertainty_error_spearman |
| --- | --- | --- | --- | --- | --- | --- | --- |
| rat_iv_plasma_clearance | scaffold_study | 42 | gine | 3 | 23 | 0.298 | 0.303 |
| rat_iv_plasma_clearance | scaffold_study | 42 | xgb | 3 | 23 | 0.295 | -0.212 |

GINE and XGBoost ensembles use the selected base-model hyperparameters with three random initializations/seeds on the same training partition. The calibration and test molecules remain excluded from every member. Variation captures optimization/subsampling variability, not all epistemic or measurement uncertainty. Ensemble size is small. Spearman associations and selective-risk curves are descriptive and may be unstable. No threshold was chosen using test errors.

## Empirical 90% interval coverage

| task | split | seed | model | nominal | n_calibration_groups | n_test_groups | radius_log10 | status | group_coverage | family_coverage | width_log10 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rat_iv_plasma_clearance | scaffold_study | 137 | dose_ridge | 0.900 | 8 | 16 |  | insufficient_calibration_groups |  |  |  |
| rat_iv_plasma_clearance | scaffold_study | 137 | gcn | 0.900 | 8 | 16 |  | insufficient_calibration_groups |  |  |  |
| rat_iv_plasma_clearance | scaffold_study | 137 | gin | 0.900 | 8 | 16 |  | insufficient_calibration_groups |  |  |  |
| rat_iv_plasma_clearance | scaffold_study | 137 | gine | 0.900 | 8 | 16 |  | insufficient_calibration_groups |  |  |  |
| rat_iv_plasma_clearance | scaffold_study | 137 | mean | 0.900 | 8 | 16 |  | insufficient_calibration_groups |  |  |  |
| rat_iv_plasma_clearance | scaffold_study | 137 | median | 0.900 | 8 | 16 |  | insufficient_calibration_groups |  |  |  |
| rat_iv_plasma_clearance | scaffold_study | 137 | mlp | 0.900 | 8 | 16 |  | insufficient_calibration_groups |  |  |  |
| rat_iv_plasma_clearance | scaffold_study | 137 | rf | 0.900 | 8 | 16 |  | insufficient_calibration_groups |  |  |  |
| rat_iv_plasma_clearance | scaffold_study | 137 | ridge | 0.900 | 8 | 16 |  | insufficient_calibration_groups |  |  |  |
| rat_iv_plasma_clearance | scaffold_study | 137 | xgb | 0.900 | 8 | 16 |  | insufficient_calibration_groups |  |  |  |
| rat_iv_plasma_clearance | scaffold_study | 2026 | dose_ridge | 0.900 | 15 | 9 | 2.528 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 5.056 |
| rat_iv_plasma_clearance | scaffold_study | 2026 | gcn | 0.900 | 15 | 9 | 2.599 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 5.199 |
| rat_iv_plasma_clearance | scaffold_study | 2026 | gin | 0.900 | 15 | 9 | 2.593 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 5.187 |
| rat_iv_plasma_clearance | scaffold_study | 2026 | gine | 0.900 | 15 | 9 | 2.588 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 5.176 |
| rat_iv_plasma_clearance | scaffold_study | 2026 | mean | 0.900 | 15 | 9 | 2.718 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 5.436 |
| rat_iv_plasma_clearance | scaffold_study | 2026 | median | 0.900 | 15 | 9 | 2.883 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 5.765 |
| rat_iv_plasma_clearance | scaffold_study | 2026 | mlp | 0.900 | 15 | 9 | 2.688 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 5.377 |
| rat_iv_plasma_clearance | scaffold_study | 2026 | rf | 0.900 | 15 | 9 | 2.981 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 5.962 |
| rat_iv_plasma_clearance | scaffold_study | 2026 | ridge | 0.900 | 15 | 9 | 2.820 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 5.640 |
| rat_iv_plasma_clearance | scaffold_study | 2026 | xgb | 0.900 | 15 | 9 | 2.739 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 5.477 |
| rat_iv_plasma_clearance | scaffold_study | 42 | dose_ridge | 0.900 | 12 | 14 | 2.103 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.205 |
| rat_iv_plasma_clearance | scaffold_study | 42 | gcn | 0.900 | 12 | 14 | 2.128 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.257 |
| rat_iv_plasma_clearance | scaffold_study | 42 | gin | 0.900 | 12 | 14 | 2.167 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.335 |
| rat_iv_plasma_clearance | scaffold_study | 42 | gine | 0.900 | 12 | 14 | 2.104 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.208 |
| rat_iv_plasma_clearance | scaffold_study | 42 | gine_ensemble | 0.900 | 12 | 14 | 2.149 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.297 |
| rat_iv_plasma_clearance | scaffold_study | 42 | mean | 0.900 | 12 | 14 | 2.179 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.357 |
| rat_iv_plasma_clearance | scaffold_study | 42 | median | 0.900 | 12 | 14 | 2.335 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.671 |
| rat_iv_plasma_clearance | scaffold_study | 42 | mlp | 0.900 | 12 | 14 | 2.089 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.178 |
| rat_iv_plasma_clearance | scaffold_study | 42 | rf | 0.900 | 12 | 14 | 2.173 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.345 |
| rat_iv_plasma_clearance | scaffold_study | 42 | ridge | 0.900 | 12 | 14 | 1.695 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 3.389 |
| rat_iv_plasma_clearance | scaffold_study | 42 | xgb | 0.900 | 12 | 14 | 2.182 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.364 |
| rat_iv_plasma_clearance | scaffold_study | 42 | xgb_ensemble | 0.900 | 12 | 14 | 2.184 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.368 |
| rat_oral_plasma_cmax | scaffold_study | 137 | dose_ridge | 0.900 | 6 | 9 |  | insufficient_calibration_groups |  |  |  |
| rat_oral_plasma_cmax | scaffold_study | 137 | gcn | 0.900 | 6 | 9 |  | insufficient_calibration_groups |  |  |  |
| rat_oral_plasma_cmax | scaffold_study | 137 | gin | 0.900 | 6 | 9 |  | insufficient_calibration_groups |  |  |  |
| rat_oral_plasma_cmax | scaffold_study | 137 | gine | 0.900 | 6 | 9 |  | insufficient_calibration_groups |  |  |  |
| rat_oral_plasma_cmax | scaffold_study | 137 | mean | 0.900 | 6 | 9 |  | insufficient_calibration_groups |  |  |  |
| rat_oral_plasma_cmax | scaffold_study | 137 | median | 0.900 | 6 | 9 |  | insufficient_calibration_groups |  |  |  |
| rat_oral_plasma_cmax | scaffold_study | 137 | mlp | 0.900 | 6 | 9 |  | insufficient_calibration_groups |  |  |  |
| rat_oral_plasma_cmax | scaffold_study | 137 | rf | 0.900 | 6 | 9 |  | insufficient_calibration_groups |  |  |  |
| rat_oral_plasma_cmax | scaffold_study | 137 | ridge | 0.900 | 6 | 9 |  | insufficient_calibration_groups |  |  |  |
| rat_oral_plasma_cmax | scaffold_study | 137 | xgb | 0.900 | 6 | 9 |  | insufficient_calibration_groups |  |  |  |
| rat_oral_plasma_cmax | scaffold_study | 42 | dose_ridge | 0.900 | 14 | 11 | 2.199 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.398 |
| rat_oral_plasma_cmax | scaffold_study | 42 | gcn | 0.900 | 14 | 11 | 2.212 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.423 |
| rat_oral_plasma_cmax | scaffold_study | 42 | gin | 0.900 | 14 | 11 | 2.109 | empirical_shift_coverage_no_exchangeability_guarantee | 0.909 | 0.964 | 4.218 |
| rat_oral_plasma_cmax | scaffold_study | 42 | gine | 0.900 | 14 | 11 | 2.205 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.410 |
| rat_oral_plasma_cmax | scaffold_study | 42 | mean | 0.900 | 14 | 11 | 2.469 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.937 |
| rat_oral_plasma_cmax | scaffold_study | 42 | median | 0.900 | 14 | 11 | 2.423 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.846 |
| rat_oral_plasma_cmax | scaffold_study | 42 | mlp | 0.900 | 14 | 11 | 2.350 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.700 |
| rat_oral_plasma_cmax | scaffold_study | 42 | rf | 0.900 | 14 | 11 | 2.281 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.562 |
| rat_oral_plasma_cmax | scaffold_study | 42 | ridge | 0.900 | 14 | 11 | 2.155 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.310 |
| rat_oral_plasma_cmax | scaffold_study | 42 | xgb | 0.900 | 14 | 11 | 2.228 | empirical_shift_coverage_no_exchangeability_guarantee | 0.909 | 0.964 | 4.456 |
| rat_oral_plasma_cmax | study | 137 | dose_ridge | 0.900 | 14 | 4 | 2.013 | empirical_shift_coverage_no_exchangeability_guarantee | 0.750 | 0.976 | 4.026 |
| rat_oral_plasma_cmax | study | 137 | gcn | 0.900 | 14 | 4 | 2.078 | empirical_shift_coverage_no_exchangeability_guarantee | 0.750 | 0.976 | 4.155 |
| rat_oral_plasma_cmax | study | 137 | gin | 0.900 | 14 | 4 | 2.046 | empirical_shift_coverage_no_exchangeability_guarantee | 0.750 | 0.952 | 4.092 |
| rat_oral_plasma_cmax | study | 137 | gine | 0.900 | 14 | 4 | 2.051 | empirical_shift_coverage_no_exchangeability_guarantee | 0.750 | 0.952 | 4.103 |
| rat_oral_plasma_cmax | study | 137 | mean | 0.900 | 14 | 4 | 2.113 | empirical_shift_coverage_no_exchangeability_guarantee | 0.750 | 0.952 | 4.226 |
| rat_oral_plasma_cmax | study | 137 | median | 0.900 | 14 | 4 | 2.030 | empirical_shift_coverage_no_exchangeability_guarantee | 0.750 | 0.952 | 4.061 |
| rat_oral_plasma_cmax | study | 137 | mlp | 0.900 | 14 | 4 | 2.136 | empirical_shift_coverage_no_exchangeability_guarantee | 0.750 | 0.952 | 4.271 |
| rat_oral_plasma_cmax | study | 137 | rf | 0.900 | 14 | 4 | 2.188 | empirical_shift_coverage_no_exchangeability_guarantee | 0.750 | 0.929 | 4.376 |
| rat_oral_plasma_cmax | study | 137 | ridge | 0.900 | 14 | 4 | 2.125 | empirical_shift_coverage_no_exchangeability_guarantee | 0.750 | 0.952 | 4.250 |
| rat_oral_plasma_cmax | study | 137 | xgb | 0.900 | 14 | 4 | 2.148 | empirical_shift_coverage_no_exchangeability_guarantee | 0.750 | 0.976 | 4.296 |
| rat_oral_plasma_cmax | study | 42 | dose_ridge | 0.900 | 14 | 9 | 2.001 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.001 |
| rat_oral_plasma_cmax | study | 42 | gcn | 0.900 | 14 | 9 | 2.088 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.176 |
| rat_oral_plasma_cmax | study | 42 | gin | 0.900 | 14 | 9 | 2.102 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.204 |
| rat_oral_plasma_cmax | study | 42 | gine | 0.900 | 14 | 9 | 2.139 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.277 |
| rat_oral_plasma_cmax | study | 42 | mean | 0.900 | 14 | 9 | 2.469 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.937 |
| rat_oral_plasma_cmax | study | 42 | median | 0.900 | 14 | 9 | 2.423 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.846 |
| rat_oral_plasma_cmax | study | 42 | mlp | 0.900 | 14 | 9 | 2.432 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.864 |
| rat_oral_plasma_cmax | study | 42 | rf | 0.900 | 14 | 9 | 2.314 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.627 |
| rat_oral_plasma_cmax | study | 42 | ridge | 0.900 | 14 | 9 | 2.203 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.405 |
| rat_oral_plasma_cmax | study | 42 | xgb | 0.900 | 14 | 9 | 2.352 | empirical_shift_coverage_no_exchangeability_guarantee | 1.000 | 1.000 | 4.704 |

Very small calibration/test group counts produce coarse estimates and wide intervals. The same test data appear at several nominal levels; these are not independent trials. No animal-level prediction intervals, conformal guarantees under covariate shift, or prospective dose-selection decisions are supported.

Risk–coverage curves rank families by ensemble SD or by one minus nearest-training similarity. They report error among retained test families retrospectively; they are not validated deployment policies.

Source: [ChEMBL](https://www.ebi.ac.uk/chembl/) 37, CC BY-SA 3.0; source publication IDs and hashes remain in the local task tables. Source-derived aggregate tables/figures in these reports are supplied under [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/) with attribution to ChEMBL and the linked source publications; they are transformed analyses, not original database releases; numerical endpoint tables, predictions and model checkpoints are not committed. Code remains under its existing license.
