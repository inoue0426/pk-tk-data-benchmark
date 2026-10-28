# Why 35,669 observations produced ten candidates

```json
{
  "input_observations": 35669,
  "source_series_ids": 22168,
  "singleton_series": 17020,
  "series_with_at_least_three_rows": 1379,
  "ND_NQ_observations": 6521,
  "missing_normalized_concentrations": 5261,
  "source_study_ids_with_more_than_one_dose": 0,
  "strict_profiles": 151,
  "strict_crossdose_series": 10,
  "censor_tolerant_profiles": 163,
  "censor_tolerant_crossdose_series": 11,
  "censor_tolerant_compounds": 24,
  "strict_profiles_less_restrictive_metadata_groups": 9,
  "censor_tolerant_less_restrictive_metadata_groups": 10,
  "dose_study_cells_with_three_quantified_times": 453,
  "potential_crossdose_groups_if_composite_sampling_confirmed": 112,
  "dtxsid_with_multiple_source_chemical_ids": 32,
  "dtxsid_with_multiple_cas_values": 0,
  "interpretation": "Sensitivity counts are not validated new series. Never pool across publication IDs; subject-level and composite profiles remain distinct."
}
```

The original parser retains every row. Source series IDs often identify individual animals with one sample; requiring three times per ID cannot reconstruct a destructive-sampling population curve. Many multi-time series contain ND/NQ concentrations, causing the original all-or-nothing curve gate to fail. ND/NQ must not be silently converted to zero; numeric-only maxima may miss the true peak and truncated AUC changes its definition. The audit measures recovery potential but does not approve these endpoints.

The metadata sensitivity removes vehicle/fasting/age/duration from grouping while retaining publication, compound, species, sex, route, matrix and regimen. Its increase is an upper bound needing protocol review, not a new verified count. Dose-specific source study IDs contain no multi-dose series; publication-level matching is necessary but insufficient. Multiple source chemical IDs can share a DTXSID, and some DTXSIDs carry multiple CAS values; identity conflicts need explicit resolution.

Proposed composite-sampling cells preserve source study and animal IDs and are saved separately. No independent publications or individual profiles are merged into verified experiments. Prioritize original NTP protocol/strain/formulation retrieval before approving these cells.
