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

- `identity_provider`: EPA_DSSTox or ChEBI for the accepted structure; candidates retain both sources. EPA acceptance requires original DTXSID/CAS and source full-key agreement, never QSAR/MS-ready substitution.
- `identity_qc_flags`: includes nonstandard source InChIKeys whose generation options were not verified. These candidates are rejected rather than normalized by editing key characters. An independently validated ChEBI mapping may still resolve the identifier; all rejected EPA evidence remains visible.
- `identity_source_member` / `identity_source_row`: master-CSV member and physical row in the bulk DSSTox ZIP. Split copies are not separately parsed.
- `pubchem_lookup_status=not_queried_optional_crossref`: structure resolved, optional CID lookup not performed. This is distinct from unresolved compound identity.

## Expanded composite audit and CEBS rows

- `curation_round=composite_expansion`: explicitly reviewed rows from eight additional CEBS workbooks. Physical row selections and expected metadata are in `cebs_expansion.json`.
- `dose_basis=actual`: oxybenzone's analytically administered dose parsed from original exception text; nominal dose remains in `original_json`.
- `endpoint_origin=reported_individual_parameter_mean`: naphthalene rat group means of individually estimated parameters. Mouse composite profiles retain `reported_composite_profile`; these methods are not pooled.
- `n_animals` and `n_animals_per_timepoint` are intentionally missing for the new imports when source pool sizes, replacements and dose-specific counts need separate reconciliation. Original count fields and sampling text remain available; no fixed count is assumed.
- `protocol_reviewed` in the audit means original source protocol material was inspected, not that all required conditions were resolved. `assessed_candidates` includes explicit unavailable-protocol decisions. `protocol_unresolved_candidates` counts candidates without an obtained/reconciled full protocol; other statuses can also withhold validation.
- `validated_composite_candidate_groups` counts linked exploratory candidate groups. Historical alias `validated_composite_series` has the same meaning; it is **not** the total validated endpoint-series count. A candidate can split into formulation-specific series, or source endpoints can have unresolved CvTdb dose linkage.
- `cebs_row_selection.csv`: one entry per workbook row, selected flag, selection rationale and raw hash. Alternative fits, metabolites and other regimens stay in raw caches.
- AUC0-8h, AUC0-inf and unspecified-terminal AUC remain distinct. Predicted Cmax is never substituted for observed Cmax. Source-reported Cmax may be affected by dose-dependent sampling schedules.

## Empirical research outputs

`results/research/endpoint_observations.csv` contains endpoint-specific eligible rows (AUC and Cmax separately). `endpoint_series_id` combines source series and endpoint; `component_id` is the transitive closure of shared molecular identity, original identity and study. It defines reference exclusion and resampling, not individual animals.

`episodes.csv` records policy, direction, support count, source IDs, target dose, observed exposure, dose ratio and support span. `one_extreme` uses one extreme support; `leave_one_dose_out` uses all other doses. Interior targets are labeled interpolation. `episode_membership.csv` is the complete support/target audit. IDs include direction so opposite-direction tasks at the same target cannot collide.

`reference_training.csv` lists disjoint reference components, records and endpoint series for every episode. `predictions.csv` includes all attempted models, predictions, fitted parameters, failure reasons, log errors, absolute relative error and fold error. `not_identifiable` is not a substituted prediction. `metrics.csv` uses equal series then equal dependency-component weighting; repeated pairs are not independent replication. `paired_comparisons.csv` compares common episodes only and records exploratory component-bootstrap and leave-one-component-out sensitivity summaries. `nonlinearity.csv` uses all doses descriptively; its labels must never be used as independent evaluation truth.

`curation_round=targeted_research_expansion` marks the parent-only K04003B addition. Its separate document ID deliberately avoids asserting the unresolved CvTdb C20619 alias. Source AUC cells remain preserved but are not benchmarked.

## Dataset A: molecular PK endpoints

`data/processed/molecular/endpoint_records.csv` retains every collected endpoint, including exclusions. `eligible_endpoints.csv` is the structure-resolved, numerical preclinical subset. `compounds.csv` deduplicates full standardized InChIKeys; connectivity keys provide conservative grouping across stereochemical ambiguity. `source_compound_ids_json` preserves aliases without name-based merging.

- `record_kind`: literature-abstracted endpoint, source aggregate endpoint, or reviewed group endpoint. These are not individual animals or concentration–time observations. Original CvTdb observations remain in `data/processed/observations.csv`.
- `original_value/unit`, `source_standard_value/unit`, `value/unit/conversion`: source and harmonized quantities. For inherited B records, original_json retains original experimental fields while the inherited endpoint value is already harmonized.
- `dose_normalized`: yes, no (reviewed B), or not_reported. Unknown is not false. A normalized reference dose is never an administered dose. `dose_mg_kg` is populated only for a single unambiguous description-level administered dose or a reviewed B group.
- `endpoint_definition`, `matrix`, `route`, `sex`: unknown remains blank/source-defined. AUC0-inf and reported windows remain separate; unspecified Vd and half-life phase are not renamed Vss or terminal half-life. Dose/route extraction from descriptions does not verify full protocols.
- `original_study_id`: only populated for reviewed B protocols. `publication_id` is a normalized DOI/PMID/document proxy, not an independent experiment; absent for unmapped PKSmart aggregates. `acquisition_publication` identifies a compilation separately.
- `original_smiles`, `source_inchikey`, `canonical_smiles`, `inchikey`, `connectivity_key`: preserve source structure while standardizing reversible protonation and retaining stereochemistry. No largest-fragment or name guess. Multicomponent structures remain excluded/unresolved.
- `eligible_A`, `strict_B_eligible`, `exclusion_reasons`, `qc_flags`: separate curation tiers; A eligibility never implies B eligibility. Overlapping exclusion reasons must not be summed as row counts.
- `raw_file`, `raw_sha256`, `source_url`, `retrieved_at`, `dataset_version`, `license`, `original_json`, `assay_parameters_json`: trace endpoints to cached source rows, original descriptions and assay metadata.

`dose_candidates_unverified.csv` contains discovery leads only. `dependency_components.csv` links connectivity groups and known publications for future leakage-aware evaluation. `prediction_cohort_feasibility.csv` separates endpoint, definition, species, route, matrix, units and normalization; it does not contain model results.

## Final presentation verification

`reports/final_verification.json` independently reconciles Dataset A full-InChIKey coverage, Dataset B protocol-series counts, and primary natural-log prediction errors with underlying CSV rows. `independent_components` counts transitive compound–study dependency groups, not episodes. `baseline_results.n_components` is endpoint/task specific. `mean_absolute_relative_error` is a fraction; presentation tables multiply it by 100. The two-page brief contains source-derived summary figures; it is not an additional endpoint dataset.
