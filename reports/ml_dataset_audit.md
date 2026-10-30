# Molecular ML eligibility audit

This experiment branch preserves the original 2,713-compound Dataset A and 59-series Dataset B. `scripts/verify_final_results.py` reconciled the original snapshot before work. All original numerical artifacts are pinned in `reports/ml/original_artifacts.json`; new task tables live under `data/processed/ml/`. Dataset A eligibility is **not** automatically ML eligibility or cross-dose validation.

## Endpoint decisions before model outcomes

The raw endpoint audit found insufficient matrix, dose or definition specificity for many records. Models predict **reported endpoints conditional on molecular structure and known administered dose**, not a universal dose-invariant molecular property. No exposure is divided by dose; no linear-PK assumption is used to construct labels.

| Endpoint | Decision and justification |
|---|---|
| AUC | No training. Available plasma records with explicit route, dose and compatible integration windows form very small strata. Unspecified AUC windows, brain/plasma mixtures and differently normalized endpoints are not pooled. |
| Cmax | Rat, oral, explicit plasma, absolute positive numerical Cmax, explicit mg/kg administered dose and original publication. Exclude fixed-time ambiguities, uncharged/unbound analytes and special regimens. Dose is a common input to molecular methods and a dose-only control. |
| Clearance | Rat, IV, explicit systemic plasma clearance, positive mL/min/kg values, known mg/kg dose and publication. Exclude infusion rates, steady-state, disease and organ-specific clearance. Known dose is retained as input; clearance is not assumed invariant across the dose range. |
| Volume | No training. ChEMBL Vd generally lacks a validated Vss definition and matrix. PKSmart Vss is source-aggregated with missing original-study/matrix fields; not merged with generic Vd. |
| Half-life | No training. Most records do not distinguish terminal, distribution or other phases. IV/oral and differently defined phases are not pooled to increase sample size. |

PKSmart [methods](https://doi.org/10.1186/s13321-025-01066-5) identify measured IV animal CL and Vss, but collected rows lack original study IDs and explicit matrix fields. Those missing fields prevent a strict combined study-aware/matrix-compatible task here; their exclusion is a benchmark decision, not a claim that the source has no scientific value. No predicted PKSmart values are used. Source-defined endpoints are derived PK parameters, not raw animal observations. Original study protocols have not been fully re-reviewed for every ChEMBL assay: retained strain, sex, formulation and sampling heterogeneity remain limitations even after the critical endpoint/dose/route/matrix gates.

## Identity, duplicates and target construction

RDKit parses each source-standardized structure. Multicomponent or invalid molecules are excluded. Full InChIKeys preserve reported stereochemistry. Unspecified stereocenters are flagged, never invented. Neutralized, canonical tautomer, stereo-free identities and original connectivity keys are joined for conservative split families; incomplete tautomer enumeration is excluded rather than treating a partial search as authoritative. Bemis–Murcko scaffolds ignore stereochemistry; all acyclic molecules share one explicit acyclic group. Original identity and canonical structures are retained.

Exact full identity / endpoint / original publication / assay / description / sex / dose cells are collapsed only when numerical values agree within a relative tolerance of 1e-6. All original IDs and record counts remain attached. Conflicting duplicate cells are excluded with explicit flags; they are not averaged into a synthetic target. Metabolizer phenotypes and explicit solution-versus-suspension comparisons are excluded as unresolved protocol covariates. Different assays or conditions are never pooled merely because the molecular identity matches. Training and scoring weight chemical families equally, so repeated records do not inflate the compound count. No global target-value outlier removal is performed.

The complete row-by-row disposition is `data/processed/ml/eligibility.csv`. Original datasets are unchanged. Selected rows and target definitions are in `tasks.csv`; no missing values are imputed into targets. Original endpoint units must be on an explicit compatible mass/molar whitelist; `uM equiv` is excluded because its analyte basis is unresolved. Every dose is independently checked against an explicit mg/kg expression in the source assay description. This catches the distinction between administered mass and infusion rate.

## Split design and interpretation

Four group designs are frozen before training: compound family, scaffold + compound family, publication + compound family, and scaffold + publication + compound family. Groups are connected transitively. Allocation uses a fixed random group order and partition row quotas, without target values. Train/validation/calibration/test quotas are 55/15/15/15%; group constraints can make actual sizes uneven. Any split with fewer than 50 training or 15 validation/calibration/test families is withheld, without searching seeds until a favorable score appears.

Compound-only and scaffold-only results are **diagnostics**: publication overlap is explicitly measured and they do not establish unseen-study generalization. Primary evidence requires publication exclusion, with scaffold exclusion where feasible. Study and scaffold+study partitions happen to coincide for the clearance task; these are executed once and not counted as independent replications. The Cmax cohort is dominated by one source publication (107 records in the original 10 mg/kg subset), making study-blocked partitions especially uneven. Independent study counts and held-out groups accompany every score.

Cross-source evaluation is unsupported: only ChEMBL passes these task gates. Do not weaken matrix or study requirements to add a second source. All models share exact split manifests. Descriptor imputation, scaling, target standardization, early stopping and hyperparameter selection use training/validation only; calibration is untouched until uncertainty analysis. No source identifiers are predictive features.

Executed counts and partition decisions: [machine-readable audit](ml/audit.json), `data/processed/ml/split_audit.json`, `data/processed/ml/splits.csv`. Counts in the machine-readable audit supersede any pre-audit size estimates.
