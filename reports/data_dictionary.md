# Data dictionary and interpretation

CSV files use UTF-8 with headers. Empty normalized values mean unavailable/unresolvable, never zero. `original_json` retains the source fields as strings (source missing values become empty strings); exact original R types and missingness remain in the cached archive. Normalized concentrations are ug/mL, doses mg/kg, times h, AUC ug*h/mL. A column suffix such as `_h` specifies its unit.

| Field(s) | Meaning |
|---|---|
| source, source_record_id | Source namespace and original observation ID or exact table/column address. |
| compound_id, compound_name, casrn | Source DTXSID where available, otherwise local chemical ID; PMC compound has a publication-local ID. Names/aliases are not definitive molecular identities. |
| canonical_smiles, inchikey | Blank until a defensible structure mapping exists. RDKit preserves stereo and fragments; it does not silently strip salts or merge tautomers. |
| species, strain, sex | Source species and sex normalized to lowercase; mixed sexes explicitly mixed. CvTdb strain is unavailable. |
| route, matrix | Oral/intravenous; blood/plasma/serum distinguished. Unknowns are flagged. Oral includes gavage; source wording retained. |
| dose_mg_kg, dose_unit | Body-weight normalized dose only. An absolute mg dose is never assumed to be mg/kg. |
| time_h, time_unit | Source normalized observation time. Nonfinite/negative time blocks curve derivation. |
| concentration_ug_ml, concentration_unit | Source normalized concentration, dimensionally converted; no molar conversion without molecular weight. Zero is retained and flagged; textual censoring is not replaced with zero. |
| study_id | CvTdb dose-specific study ID, or DOI-based publication study for PMC; do not interpret these as comparable across sources. |
| document_id | CvTdb extraction document / publication proxy. Not proof of shared protocol. |
| experiment_id, subject_id | CvTdb series and subject IDs as supplied; may represent summary series or single points, not independent animals. PMC table column is a dose/sex group mean. |
| protocol_id | Explicit article-specific protocol for rat single oral dosing. Absent for CvTdb. |
| n_animals | Reported subject/group count when numeric; not inferred from point counts. |
| regimen, vehicle, fasting, age, age_unit, dose_duration, dose_duration_unit | Comparability covariates preserved where available. These participate in CvTdb candidate grouping. |
| parent_analyte, radiolabeled | Whether dosed and analyzed chemical source IDs match; source radiolabel status. Metabolites and confirmed radiolabel records do not enter derived curves. |
| publication_doi, publication_url | Original study citation where supplied. |
| source_url, retrieved_at, dataset_version, raw_sha256 | Download provenance; UTC time, distribution/snapshot version and raw-file digest. |
| original_json | All original CvTdb export fields, or exact selected PMC endpoint cells/design row and table/column mapping. |
| qc_flags, eligible_for_curve | Semicolon-delimited concerns and computational eligibility, not automatic scientific approval. |
| auc, auc_sd, auc_unit, auc_definition, auc_terminal_time_h | Reported group-mean AUC and SD, with source definition. AUC0-tn endpoint is unknown per group; not silently equated to AUC0-8h or AUC0-inf. |
| cmax, cmax_sd, cmax_unit | Reported mean and SD for PMC; sampled curve maximum in curve_qc (no SD). These origins remain separate. |
| endpoint_origin | reported_group_mean vs computed_from_source_curve. |
| sampling_schedule_h | Planned sample times from article Table 5; not concentration observations or proof that every group had quantifiable values at every time. |
| has_numeric_curve | True only for usable downloaded numerical concentration observations, not a published plot. |
| high_confidence_cmax, high_confidence_auc | Endpoint-specific eligibility. Cmax only in the initial PMC subset; no AUC benchmark accepted. |
| candidate_group_id / candidate_series_id | Stable hash of grouping metadata. Candidate match does not establish scientific comparability. |
| n_doses, dose_min, dose_max, dose_ratio, doses_json | Number and range of distinct normalized doses; replicate curves do not create additional dose levels. |
| comparable_auc_window | True only if computed partial AUC definitions (including endpoints) agree exactly. Does not resolve missing protocol. |
| metadata_completeness | Fraction of grouping fields present for CvTdb; PMC protocol grouping fields fully present. Does not score structure availability or uncertainty. |
| time_start_h, time_end_h | Actual curve integration interval; no artificial t=0. |
| source_record_ids, experiment_ids, study_ids | JSON arrays linking aggregate outputs back to auditable records. |
| episode_id, role, regime | Within-series Task A episode, support/target role and extrapolation label. |
| observed, predicted, log_error | Cmax in ug/mL and natural log(predicted/observed). |
| absolute_relative_error | abs(predicted/observed - 1); unitless. |
| n_support, dose_ratio (predictions) | Number of lower support doses and target/max-support ratio. |

QC flags cover invalid or unsupported values, source/original conversion disagreements, ambiguous routes, missing metadata, duplicate record IDs, inconsistent metadata inside a source series, repeated/unknown regimen, human observations, nonparent analytes, radiolabeling and inadequate sampling. Missing structure blocks molecular evaluation; missing strain/protocol blocks CvTdb high-confidence claims. A missing flag is not a claim that assay accuracy, formulation, disease state or full protocol has been independently verified.

## Expanded outputs

`reported_endpoints_expanded.csv` contains all selected PMC groups, including explicitly flagged groups; `high_confidence_endpoints.csv` contains groups belonging to primary comparable series. Eligibility is endpoint-specific (`high_confidence_auc`, `high_confidence_cmax`). `crossdose_series_expanded.csv` separates `protocol_verified` (explicit within-publication protocol) from `primary_eligible` (acceptable exposure comparability). A publication can contribute correlated route/sex/compound strata; series counts do not equal independent publications.

`molecular_group_id` is `InChIKey:<full key>` for an accepted mapping, otherwise the unchanged source compound ID. It never removes the study/protocol from grouping. `resolved_cas_json` contains ChEBI's CAS references; source CAS values remain separate. `pubchem_cids_json` retains every exact-structure UniChem mapping; `pubchem_cid` is populated only for a singleton. `identity_status` distinguishes `resolved_exact`, `not_found`, `not_queried`, invalid CAS, identifier conflict, ambiguous structures, incomplete search and lookup failure. `compound_identity_candidates.json` retains rejected and ambiguous structure alternatives and provenance.

`structure_input_kind=source_inchi_fallback_invalid_smiles` means the source SMILES failed RDKit but the same record's standard InChI produced a molecule matching its full expected InChIKey. The invalid original SMILES is retained in the candidate record. No guessed structure, salt stripping or stereoisomer collapse is used.

Expanded AUC uses explicit AUC0-inf where eligible; mass concentrations are ug/mL and AUC ug*h/mL. Memantine's retained Cmax remains nM and its ambiguous AUC is not imported. No mass/molar conversion is inferred. `original_json` records source cell, row label, column and conversion factor; dose and animal-count selections are in the hashed curation configuration. `protocol_evidence_json` records matching source text fragments.

CvTdb `cvtdb_profile_audit.csv` and `cvtdb_composite_candidates.csv` are audit-only recovery candidates, not fabricated animal timecourses. Original subject, source series and study IDs remain in the observation table.

## Second expansion additions

- `source=cebs`: original XLSX group endpoints. `source_record_id` includes experiment, sheet and physical Excel row. `original_json` retains every column, including nominal dose, predicted Cmax, model choice and source exceptions.
- `dose_basis=analytically measured actual`: benzophenone dose from the source's actual-dose comment, with nominal dose retained in `original_json`.
- `dose_volume_ml_kg`: explicit CEBS administration volume, checked within source strata; not inferred for PMC records.
- `n_animals_per_timepoint`: source composite-sampling count. `n_animals` is the reported total group count, never the number of longitudinal profiles.
- `endpoint_origin=reported_composite_profile`: reported population-profile endpoint, with no invented individual AUC/Cmax.
- `license_url` / `license_noncommercial`: source reuse terms. Lorlatinib's CC BY-NC restriction persists in derived rows.
- `regimen` and protocol day: methyl eugenol day 1 versus repeated-dose day 28 remain separate correlated series, sharing the same study ID.
- `reports/composite_audit.csv`: all 112 candidate group IDs, sampled-study decisions, raw observation metadata and matched endpoint IDs. `not_reviewed` is not an exclusion or validation.
- `primary_independent_studies`: distinct source experiment IDs/PMC publications in the primary subset. `primary_independent_publications` counts PMC publications only; CEBS chemical-collection DOIs are not independent publications.
- `resolved_unique_parent_structures` and `unresolved_parent_identifiers`: report molecular resolution separately; their sum is not proven unique molecular diversity.
