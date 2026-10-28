# Current QC report

All original observations and all selected curated endpoints are retained. Exclusions apply to primary eligibility, not source deletion.

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
    "curated_endpoint_groups": 47,
    "primary_endpoint_groups": 27,
    "primary_comparable_series": 9,
    "flagged_protocol_series": 8,
    "extractions_failed": 0,
    "source_overlap_groups": 0,
    "identity_status_counts": {
      "not_queried": 238,
      "resolved_exact": 13,
      "not_found": 2
    }
  }
}
```

Article/table hashes, exact cell addresses, dose columns, evidence phrases and units are checked against reviewed configurations. Censored CvTdb values remain censored; audit-only recovery does not imply validated AUC. Primary grouping requires publication, protocol, compound, species, strain, sex, route, matrix, regimen and formulation conditions. Source DOI overlaps are flagged; no records from independent publications are pooled. Chemical resolution accepts only unique full-InChIKey agreement; conflicts retain candidate lists. Exact CAS checksum is validated.

Flags remain for unavailable fasting, AUC tail fractions, aggregate-only uncertainty and missing structures. Reported AUC0-inf is distinct from observed-window AUC. Rat HY-071085 AUC0-tn is ineligible. Oral formulation changes and ambiguous Cmax units in PMC6160927, incomplete tofacitinib oral formulation evidence, and memantine composite/dose-volume conditions are outside the primary subset.
