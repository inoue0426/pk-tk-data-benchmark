# Current QC report

No source records are silently deleted. Eligibility applies to curated endpoints, not removal from source caches.

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
  "potential_overlap_groups": 0,
  "expansion": {
    "curated_endpoint_groups": 195,
    "primary_endpoint_groups": 169,
    "primary_comparable_series": 57,
    "flagged_protocol_series": 10,
    "extractions_failed": 0,
    "source_overlap_groups": 11,
    "identity_status_counts": {
      "resolved_exact": 228,
      "unresolved_or_identifier_conflict": 21,
      "ambiguous_multiple_structures": 3,
      "not_found": 2
    },
    "reviewed_composite_candidates": 97,
    "validated_composite_series": 31
  }
}
```

Exact source/table checksums and protocol phrases protect reviewed mappings. All 112 exploratory groups have an audit disposition; see composite_audit.json for original-protocol review and unresolved counts. Group endpoints are not individual measurements. Dose/matrix conflicts are withheld. Formamide Cmax_pred is not Cmax_obs. Benzophenone uses actual dose comments, not nominal dose. Isoeugenol rat dose conflicts prevent primary eligibility. AUC units and definitions are endpoint-specific; 0-t without a verified terminal time is not compared. Source exceptions, license terms, original cells and predicted values are retained.
