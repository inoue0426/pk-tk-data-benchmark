# Graph neural-network experiments

These are molecular-plus-known-dose predictions of reported rat plasma PK endpoints, not dose extrapolation, not a validation of the full Dataset A, and not Dataset B transfer. All trained records come from ChEMBL. Original A/B tables and original empirical cross-dose results are unchanged.

**Measured findings.** Strict scaffold+study mean log10 MAE: clearance GINE 0.337, XGBoost 0.370, descriptor Ridge 0.474; Cmax GINE 0.858, XGBoost 0.888, descriptor Ridge 0.880. These averages describe overlapping splits, not independent replications. Model rankings are exploratory; paired component comparisons below determine the limits of inference.

Random-compound Cmax: Random Forest 0.425, GINE 0.499; scaffold-only Cmax: XGBoost 0.589, GINE 0.703. These diagnostic splits share publications across partitions and do not establish unseen-study generalization. Negative R² compares squared error with the held-out test mean, not with the separately trained mean predictor. Ensemble uncertainty and empirical interval coverage are reported separately; small group counts and distribution shift prevent coverage guarantees.

**Interpretation.** The training-mean control has clearance MAE 0.326, below all molecular models on average. All models have negative mean strict-split log-space R² for both endpoints. Molecular inputs therefore do not establish useful unseen-study clearance prediction here. Of 5 prespecified GINE–XGBoost strict comparisons, 2 exploratory fixed-prediction component-bootstrap intervals exclude zero in favor of GINE; the remaining intervals include zero. These unadjusted, overlapping-split results do not establish general graph superiority, and outperforming XGBoost alone does not establish improvement over simpler controls. Cmax becomes substantially harder when both studies and scaffolds are excluded; partition populations differ, so this is not a causal estimate of a scaffold effect.

gine_ensemble: 90% nominal intervals use 12 calibration groups and 14 test groups; observed group coverage 100.0%, radius 2.149 log10 (×/÷141). Such wide intervals limit usefulness despite high observed coverage; distribution-shift exchangeability is unverified.

xgb_ensemble: 90% nominal intervals use 12 calibration groups and 14 test groups; observed group coverage 100.0%, radius 2.184 log10 (×/÷153). Such wide intervals limit usefulness despite high observed coverage; distribution-shift exchangeability is unverified.

## Executed architectures

| model | runs | min_parameters | max_parameters | median_best_epoch | training_seconds |
| --- | --- | --- | --- | --- | --- |
| gcn | 19 | 5665.000 | 19521.000 | 33.000 | 29.473 |
| gin | 19 | 7779.000 | 27843.000 | 30.000 | 24.940 |
| gine | 19 | 8675.000 | 29635.000 | 31.000 | 32.233 |
| mlp | 19 | 66689.000 | 134401.000 | 17.000 | 5.020 |

PyTorch Geometric implements two-layer GCN, GIN and edge-aware GINE with hidden width 32 or 64, ReLU, concatenated mean/sum pooling, and a dose-conditioned prediction head. Sum pooling is divided by a fixed 20; this is not a fitted normalization. Atom features encode element, degree, formal charge, hybridization, CIP chirality, aromaticity, ring membership, hydrogen count and mass. Bond features encode bond type, stereo, conjugation and ring membership. GCN/GIN use connectivity but not bond attributes; GINE uses attributes. The explicit GINE ablation retains only atom type and connectivity with constant edge features.

All runs use CPU, two threads, AdamW, weight decay 0.001, dropout 0.1, gradient clipping and validation-only early stopping. Target centering/scaling and descriptor transforms use training observations only. Training losses weight repeated chemical families equally. Parameter counts, runtimes and learning histories are recorded. Tests verify atom-order invariance and that changing held-out labels cannot change fitted graph predictions.

## Prespecified GINE versus XGBoost

| task | split | seed | model | n_test_families | n_dependency_groups | delta_family_mae_log10 | delta_component_mae_log10 | ci_low | ci_high | interval_status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rat_iv_plasma_clearance | scaffold_study | 42 | gine | 23 | 14 | 0.009 | 0.006 | -0.021 | 0.032 | exploratory_component_bootstrap_fixed_predictions |
| rat_iv_plasma_clearance | scaffold_study | 137 | gine | 23 | 16 | -0.120 | -0.147 | -0.267 | -0.032 | exploratory_component_bootstrap_fixed_predictions |
| rat_iv_plasma_clearance | scaffold_study | 2026 | gine | 28 | 9 | 0.010 | 0.005 | -0.012 | 0.021 | exploratory_component_bootstrap_fixed_predictions |
| rat_oral_plasma_cmax | scaffold_study | 42 | gine | 28 | 11 | -0.060 | -0.028 | -0.118 | 0.073 | exploratory_component_bootstrap_fixed_predictions |
| rat_oral_plasma_cmax | scaffold_study | 137 | gine | 52 | 9 | -0.000 | -0.083 | -0.134 | -0.031 | exploratory_component_bootstrap_fixed_predictions |

Differences are GINE minus XGBoost; negative favors GINE. Component-bootstrap intervals resample fixed paired errors over held-out linked study/identity components, only with at least eight groups. They do not refit models, are exploratory and unadjusted for multiple comparisons. Too-small splits have no interval. No graph-superiority claim follows from a single favorable score.

## Atom/bond feature ablation

| task | split | seed | model | ablate | mae_log10 | rmse_log10 | n_train_families |
| --- | --- | --- | --- | --- | --- | --- | --- |
| rat_iv_plasma_clearance | scaffold_study | 42 | gine | False | 0.304 | 0.350 | 92 |
| rat_iv_plasma_clearance | scaffold_study | 42 | gine | True | 0.290 | 0.344 | 92 |

This one-seed comparison cannot establish a general feature effect. No pretrained encoder or multitask model was used. Small, heterogeneous endpoint cohorts do not justify increased network complexity merely to pursue positive results.

Implementation references: [PyG GINEConv](https://pytorch-geometric.readthedocs.io/en/latest/generated/torch_geometric.nn.conv.GINEConv.html), [PyG installation](https://pytorch-geometric.readthedocs.io/en/2.7.0/install/installation.html).

Source: [ChEMBL](https://www.ebi.ac.uk/chembl/) 37, CC BY-SA 3.0; source publication IDs and hashes remain in the local task tables. Source-derived aggregate tables/figures in these reports are supplied under [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/) with attribution to ChEMBL and the linked source publications; they are transformed analyses, not original database releases; numerical endpoint tables, predictions and model checkpoints are not committed. Code remains under its existing license.
