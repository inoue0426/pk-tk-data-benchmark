# Executed QC report

No downloaded observation was silently deleted. Flags can overlap; counts must not be added as mutually exclusive exclusions. All CvTdb records lack strain/structure and verified shared protocol.

```json
{
  "input_observations": 35669,
  "retained_audit_observations": 35669,
  "flagged_observations": 35669,
  "excluded_from_curve_derivation": 23216,
  "observation_flags": {
    "missing_strain": 35669,
    "missing_structure": 35669,
    "unverified_protocol": 35669,
    "nonparent_analyte": 13788,
    "invalid_concentration": 11785,
    "missing_doi": 10755,
    "missing_external_identifier": 3931,
    "missing_sex": 2744,
    "not_single_dose": 2628,
    "radiolabeled_analyte": 889,
    "not_preclinical": 762,
    "radiolabel_status_unknown": 694,
    "invalid_dose": 396,
    "zero_concentration_censoring_unknown": 96,
    "extreme_dose_review": 30,
    "invalid_time": 5
  },
  "curve_groups": 22168,
  "curve_group_flags": {
    "insufficient_sampling_times": 20791,
    "blocked_observation": 13035,
    "nonpositive_exposure": 36,
    "duplicate_sampling_times": 16
  },
  "usable_curve_groups": 151,
  "pmc_selected_groups": 6,
  "pmc_cmax_eligible_groups": 6,
  "pmc_auc_benchmark_eligible_groups": 0,
  "pmc_scope_exclusions": "Other Table 1 columns (IV/repeated dosing), all dog Table 2 columns, tissues and supplements are out of selected adapter scope, preserved in raw XML. Dog table conflict documented.",
  "potential_overlap_groups": 0
}
```

Partial AUC uses linear trapezoids only from the first to last observed time; no t=0 or terminal tail is imputed. Duplicate timepoints block curve integration instead of being averaged. Sampled maxima are not necessarily true Cmax. Unsupported units remain missing and flagged. Original strings, identifiers and units are retained in original_json and the raw archive.

The PMC rat table is checked against dose/sex/design headers. Its AUC0-tn endpoints remain available but cannot enter an AUC benchmark until terminal-time comparability is established. Missing structure does not invalidate an empirical within-compound Cmax episode but prevents molecular/scaffold evaluation.

Potential overlaps use DOI plus name/identifier; absent matches do not prove independence. Cross-source structural identity resolution remains unfinished.
