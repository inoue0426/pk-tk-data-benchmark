# One-page research summary

This report concerns **Dataset B**. The separate 2026-10-29 [Dataset A molecular PK expansion](molecular_pk_expansion.md) reports broader structure/endpoint coverage; it does not increase validated cross-dose counts.

**Question:** When can initial preclinical PK/TK measurements support dose extrapolation?

**Data:** 59 comparable series, 21 compounds, 19 studies; 42 AUC series and 33 Cmax series. 513 dose-prediction episodes represent only 19 linked compound–study units overall. Species: hamster, mouse, rat. Matrices and dosing protocols remain separate.

**Design:** One-support scaling and LODO interpolation/extrapolation; power laws require two support doses. Pooled exponents exclude the entire target compound–study component and match endpoint/species/route/matrix/AUC definition/regimen. No target-exposure leakage or neural models.

**Two-support low→high results** (equal component weighting; natural-log errors; relative error is a fraction):

| endpoint | model | n_components | mae_ln | rmse_ln | mean_absolute_relative_error | within_twofold |
| --- | --- | --- | --- | --- | --- | --- |
| auc | linear | 11 | 0.525 | 0.700 | 0.381 | 0.701 |
| auc | power | 11 | 0.583 | 0.808 | 0.619 | 0.667 |
| cmax | linear | 12 | 0.601 | 0.793 | 0.509 | 0.615 |
| cmax | power | 12 | 0.529 | 0.716 | 0.601 | 0.646 |

**Measured findings:**

- **AUC, two-support low→high:** linear MAE(log) 0.525; power 0.583. Paired power-minus-linear difference 0.058; exploratory component-bootstrap interval [-0.176, 0.308]. This interval includes zero; no reliable superiority is established.
- **CMAX, two-support low→high:** linear MAE(log) 0.601; power 0.529. Paired power-minus-linear difference -0.072; exploratory component-bootstrap interval [-0.391, 0.231]. This interval includes zero; no reliable superiority is established.

**Interpretation:** Dose-normalized exposure varies across selected compounds. A fitted local exponent can help or amplify errors; paired comparisons across both directions and interpolation are in the full report. Descriptive nonlinear trends are not statistical proof of saturation. Curvature and prediction intervals are largely unsupported by sparse aggregate measurements.

**Independent evaluation:** Compound/study-aware bootstrap summaries are exploratory; thousands of correlated predictions do not create thousands of independent experiments. Scaffold-based molecular ML remains premature. Sampling schedules, uncertain AUC tails, missing individual variability and source-selection bias remain material limitations.

**Next priority:** Independent compounds with at least four comparable doses and well-defined AUC/Cmax, before flexible molecular models. Figures: `results/figures/`; complete methods/results: [scientific_analysis.md](scientific_analysis.md). Exact configuration and hashes: [research_manifest.json](research_manifest.json).
