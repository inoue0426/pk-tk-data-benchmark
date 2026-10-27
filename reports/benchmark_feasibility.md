# Executed benchmark feasibility

This is a working collection and QC pipeline with a **one-compound Cmax pilot**, not a publication-ready molecular benchmark. 6 reported rat group means form 2 sex-specific series with three doses each.

- **A**: executed_Cmax_pilot
- **B**: not_generated: only one high-confidence compound; molecular identity unresolved
- **C**: not_generated: no verified molecular structures; no scaffold coverage
- **D**: descriptive slopes only; no independent nonlinear labels or negative controls
- **leakage_checks**: support/target IDs disjoint; target doses strictly higher; protocol and identity match
- **dependence**: episodes reuse groups and targets; never treat predictions as independent replicates
- **intervals**: not estimated: no independent calibration set
- **complex_models**: RF/XGBoost/ChemProp and saturable models withheld: one compound, at most two support doses

Linear scaling uses the geometric mean of support exposure/dose; power law fits log exposure against log dose using two support doses. Low-only episodes predict both higher doses; two-dose episodes predict the highest dose. Both models are compared on the same two-dose episodes. Natural logarithms are used. There is no hyperparameter selection.

```json
[
  {
    "model": "linear_scaling",
    "n_support": 1,
    "aggregation": "model/n_support",
    "n_predictions": 4,
    "mae_log": 1.79643454213236,
    "rmse_log": 1.8401710943061522,
    "mean_absolute_relative_error": 0.8204026983638295,
    "n_independent_compounds": 1,
    "interval_coverage": null
  },
  {
    "model": "linear_scaling",
    "n_support": 2,
    "aggregation": "model/n_support",
    "n_predictions": 2,
    "mae_log": 1.4742221105046842,
    "rmse_log": 1.4793851682351744,
    "mean_absolute_relative_error": 0.7692952648609124,
    "n_independent_compounds": 1,
    "interval_coverage": null
  },
  {
    "model": "power_law",
    "n_support": 2,
    "aggregation": "model/n_support",
    "n_predictions": 2,
    "mae_log": 0.6444248632553506,
    "rmse_log": 0.727789527851821,
    "mean_absolute_relative_error": 1.0148879961665442,
    "n_independent_compounds": 1,
    "interval_coverage": null
  }
]
```

One-compound results cannot establish molecular ML performance. There is no independently calibrated interval or nonlinear classifier. Descriptive slopes use all three doses and are not held-out detection accuracy. No significance or generalization claim is justified. Metrics by dose ratio, compound and species are in results/metrics.csv.

Task B requires verified identities and enough independent compounds and studies for development/held-out sets; allocate whole identity–study connected components, not rows. Task C additionally requires verified structures and whole scaffold components (including acyclic structures as a shared group). Task D needs independent dose-response labels, negative controls and uncertainty; this pilot supplies descriptive examples only.

Next priorities: adjudicate CvTdb publication protocols and strains, recover original compound structures and identity aliases, obtain matched-window AUC or individual data, investigate the PMC dog-table error, then add independent compounds before any molecular model or compound/scaffold split.
