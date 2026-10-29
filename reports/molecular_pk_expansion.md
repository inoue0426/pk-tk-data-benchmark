# Dataset A molecular PK expansion — 2026-10-29

## Verified result and dataset boundary

**Dataset A contains 2,713 unique standardized full-InChIKey structures with numerical preclinical PK endpoints** (2,647 connectivity groups). The 500-structure target is **met**. These are source-curated numerical measurements with explicit in vivo evidence and RDKit-validated structures, not an assertion that every original animal protocol has been reviewed.

Dataset B remains **59 validated series, 21 compounds and 19 independently reviewed source studies**. No new B series were promoted and its numerical benchmark was not rerun. 19 B compounds have resolved structures and contribute to A. Records missing structures are retained in the exclusion ledger. A compound appearing in A never establishes cross-dose eligibility.

## Coverage

| metric | count |
| --- | --- |
| unique_structured_compounds | 2713 |
| unique_connectivity_groups | 2647 |
| compounds_with_auc | 995 |
| compounds_with_cmax | 994 |
| compounds_with_both_auc_cmax | 355 |
| compounds_with_clearance | 1182 |
| compounds_with_volume | 427 |
| compounds_with_half_life | 1250 |
| collected_endpoint_records | 11592 |
| eligible_endpoint_records | 6459 |
| excluded_endpoint_records | 5133 |
| known_source_publications | 627 |
| independent_experimental_studies | None |
| original_protocol_verified_studies | 17 |
| compounds_with_verified_comparable_doses | 19 |
| dose_normalized_records | 24 |
| explicit_dose_records | 3843 |
| target_500_obtained | True |
| compounds_with_three_verified_comparable_doses | 17 |

`independent_experimental_studies = null` means unknown, not zero. Dataset A has 627 deduplicated publication/source-study proxies with eligible endpoints; compilation rows lacking original citations are not counted as additional studies. Only the inherited 19 Dataset B studies have original-protocol validation. A numeric dose extracted from an abstracted assay description is not a reliable original protocol.

Dataset B: 21 compounds have ≥2 comparable doses; 19 have ≥3. 15 compounds support eligible AUC, 12 support Cmax and 6 support both somewhere in the validated cohort; 16 *individual series* support both. AUC/Cmax eligibility is endpoint-specific.

| source | collected_records | eligible_records | compounds | auc_compounds | cmax_compounds | source_publications |
| --- | --- | --- | --- | --- | --- | --- |
| chembl | 10000 | 4885 | 2337 | 980 | 984 | 610 |
| pksmart | 1369 | 1369 | 371 | 0 | 0 | 0 |
| strict_b | 223 | 205 | 19 | 15 | 10 | 17 |

Per-source compound counts overlap and must not be summed. Identity uses standardized full InChIKey, with source structures retained; stereochemistry is preserved. Connectivity groups provide a conservative sensitivity count and future split grouping. No names alone are merged. Mixtures and invalid structures are excluded from A and remain in the ledger. Reversible protonation is normalized; original structures and source keys remain available. No tautomer enumeration or largest-fragment guessing is used.

## Verified source findings

- **ChEMBL 37:** public CC BY-SA 3.0 numerical activity API plus original assay/document metadata. Fixed first 10,000 activity IDs from a rat-target ADME query for CL, Vd, T1/2, AUC and Cmax; the query contains 56,033 records. This bounded tranche is neither exhaustive nor random. API filters alone are insufficient: assay organism, source validity, endpoint descriptions, units, censoring and explicit in vivo evidence are independently checked. [API documentation](https://chembl.gitbook.io/chembl-interface-documentation/web-services) and [license](https://chembl.gitbook.io/chembl-interface-documentation/about).
- **PKSmart:** downloaded the publisher's CC BY 4.0 workbook and imported Table S2 only (371 source rows). This table reports raw-scale CL and Vss despite the paper's description of later log-transformed modeling. Units are interpreted from table headers; values are not exponentiated. Human tables, fu-only data and Table S6 predicted animal values are excluded. Missing original dose/publication metadata prevents cross-dose use. [Paper](https://doi.org/10.1186/s13321-025-01066-5), [licensed supplement](https://doi.org/10.6084/m9.figshare.30224986.v1).
- **2026 9,098-compound PBPK study:** internal Roche gRED data shared with pRED; dose-normalized AUC and median repeated endpoints. The listed supplement provides ECCS thresholds, not 9,098 structures/measurements. No compound-level dataset was obtained or counted. [Original article](https://pmc.ncbi.nlm.nih.gov/articles/PMC12958287/).
- **NCATS FRDB:** official ZIP successfully downloaded and inspected: 13,455 PK rows, all Homo sapiens. None count toward preclinical A. The archive readme does not establish a blanket redistribution license; no human rows are republished. [Official download](https://drugs.ncats.io/downloads-public).
- **REFRAMEdb:** the public catalog advertises PK annotations but is a clinical-stage library, and catalog PK counts do not prove preclinical numerical availability. The accessible client page did not establish an openly licensed bulk numerical PK export. No annotation count or screening dose-response is imported as PK. [Portal](https://reframedb.org/).
- **CEBS/PMC and PK-DB:** existing reviewed CEBS/PMC endpoints are incorporated once via B; no new CEBS series claimed. PK-DB's documented account/download and reuse restrictions remain unchanged; no retry or bypass was attempted. [Source inventory](data_inventory.md).

## Harmonization and exclusions

| reason | records |
| --- | --- |
| ambiguous_route_or_matrix | 836 |
| censored_or_interval_value | 695 |
| endpoint_description_not_confirmed | 628 |
| invalid_or_missing_structure | 25 |
| missing_numeric_endpoint | 396 |
| multicomponent_structure_unresolved | 601 |
| no_explicit_in_vivo_evidence_or_non_pk_context | 3091 |
| nonpositive_or_nonfinite_endpoint | 34 |
| source_data_validity_flag | 107 |
| species_mismatch_or_missing | 43 |
| unsupported_or_conflicting_unit | 617 |

Flags overlap; count excluded rows from the summary rather than summing reasons. Original values/units, standard values/units, descriptions, assay parameters, API IDs, URLs, hashes, retrieval dates and source licenses remain available. Unsupported dimensions, nonpositive/nonfinite values, interval/censored measurements, invalid structures, conflicting identifiers and absent/contradictory in vivo context remain excluded. Source-suspected duplicates and incomplete experimental fields remain flagged. No endpoint is computed from a concentration curve and no missing endpoint is imputed.

Dose-normalization has three explicit states: yes, no (reviewed B), and not_reported. Bodyweight-normalized CL/Vd is not labeled dose-normalized. A dose-normalized exposure never supplies a measured administered dose; no multiplication by a guessed dose is performed. Cmax nM conversions use the source structure molecular weight and are flagged. AUC definition, matrix and route remain distinct. Oral clearance retains an apparent-clearance warning; unspecified Vd/half-life phase is not relabeled Vss/terminal half-life.

## Original-dose search and strict eligibility

The metadata screen found **11 unverified endpoint-specific multi-dose candidates**, written to `dose_candidates_unverified.csv`. Matching publication/compound/species/sex/route/matrix is only a discovery lead: original protocol, formulation, regimen, sampling and AUC definitions still require review. These candidates contribute **zero new validated B series**. PKSmart's original 2013 compilation is cited by the paper, but per-record original study/dose mapping is absent from Table S2; it cannot establish dose series.

## Compound-disjoint modeling feasibility

There are 1467 scaffold groups (including one shared empty scaffold for 18 acyclic compounds). The known-publication subset has 577 connected compound–publication components; the largest contains 109 connectivity groups. These counts describe potential dependency structure, not independent animal experiments.

Large-scale molecular modeling is now structurally possible for sufficiently populated homogeneous endpoint/species/route cohorts. `prediction_cohort_feasibility.csv` separates definitions, units and dose-normalization status. Before freezing evaluation, resolve within-source and cross-source duplicates, review original study metadata and restrict exposure modeling to compatible definitions and known doses. Group all records of a molecular identity, connectivity variants and linked publications together. PKSmart rows with unknown original studies cannot support a study-disjoint evaluation claim. AUC/Cmax at unknown doses cannot be treated as unconditional molecular properties. No models or train/test splits were created in this expansion.

## Reproduce and continue

Run `PYTHONPATH=src python -m pktk.molecular` once online after the existing strict pipeline; use `--offline` for cached replay. Run `PYTHONPATH=src python scripts/verify_molecular_reproducibility.py` to compare A artifacts and verify B numerical artifacts are unchanged. Raw caches and generated data remain ignored by Git. Source-specific licenses apply to derived datasets; do not distribute the combined data under the code MIT license.

Next: adjudicate the dose-candidate queue against original papers; map PKSmart original records; review larger, versioned ChEMBL tranches across additional species. Source study validation takes precedence over increasing the structure count. Dataset B analysis remains in [scientific_analysis.md](scientific_analysis.md).
