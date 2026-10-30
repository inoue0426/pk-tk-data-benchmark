# Molecular PK research execution log

This sprint began on 2026-10-30. The 72-hour authorization is an upper limit, not elapsed time claimed. Actual event timestamps and per-run training seconds are recorded in `results/ml/events.jsonl` and run metrics. No paid services, proprietary data or remote merges are used.

## Checkpoint 1 — audit and environment

- Verified original counts and baseline metrics; created `experiments/molecular-gnn-benchmark` from finalized main.
- Audited all five endpoints, source metadata and original assay descriptions. Rejected heterogeneous AUC/half-life/volume tasks and PKSmart rows with unresolved critical matrix/study metadata for this strict benchmark.
- Detected infusion-rate descriptions among nominal mg/kg clearance rows and fixed-time/uncharged Cmax ambiguities; excluded them in the new ML layer without changing Dataset A.
- Installed scikit-learn, XGBoost, PyTorch and PyTorch Geometric into the local virtual environment. Hardware: macOS arm64, ten logical CPUs; PyTorch reports no available CUDA/MPS accelerator. Experiments use CPU with two threads.
- XGBoost initially failed to load because OpenMP was not on its library search path. Preloading with ctypes did not resolve it. The launch script exposes the existing PyTorch wheel's local `libomp.dylib` through `DYLD_LIBRARY_PATH` before process startup; no system installation or private-directory modification is needed.
- Fixed split rules, grids and focused learning-curve/feature/ensemble follow-ups before inspecting model outcomes. `configs/ml/benchmark.json` and committed experiment code define the budget. No test-based tuning.

## Checkpoint and continuation protocol

Run `PYTHONPATH=src .venv/bin/python -m pktk.ml.data` to reproduce the task audit, then `.venv/bin/python scripts/run_molecular_ml.py --stage baselines` and `--stage graphs`. Follow with `--stage ablations` and `--stage ensembles`. Each completed run has a manifest, validation-only tuning history, predictions and metrics under `results/ml/runs/`. Checkpoints and numerical data remain ignored by Git.

A completed run is reused only if source table, split, configuration and training-code hashes match. Interrupted runs are rerun from their declared seed; completed runs are not refitted. Incompatible code/data changes require a separately versioned experiment rather than silently reusing stale results. Candidate best-epoch neural weights and preprocessing objects are cached locally. The stage queues and failure JSON files record pending/failed jobs. No background process is claimed to survive termination of the execution environment.

## Implementation diagnostics

- The first baseline smoke pass exposed pandas read-only NumPy views during in-place weight normalization; replaced in-place division with new arrays before any successful fits. All failed attempts remain in the event log. The smoke-test limit now counts attempts, so a failing method cannot consume the whole queue.
- RDKit tautomer enumeration reached its limit for a few structures. These rows are explicitly excluded; canonicalization was not silently accepted from an incomplete enumeration.

## Checkpoint 2 — duplicate/protocol correction

Before interpreting model scores, review of every multi-record cell found metabolizer phenotype strata and solution/suspension comparisons that could not be collapsed. The v1 pilot was archived under results/ml_superseded_pilot_v1 and data/processed/ml_superseded_pilot_v1, and its analysis summaries were removed from current reports. Version 2 excludes phenotype/formulation-confounded records, preserves assay/description/sex grouping, and excludes conflicting duplicate cells instead of taking a median. Every split and model is rerun under the corrected eligibility gate. This correction responds to source metadata, not test performance. Original Dataset A/B artifacts are unchanged.

## Checkpoint 3 — explicit unit basis and stable split IDs

The final source-unit audit found three Cmax labels with original units `uM equiv`. Their parent-drug versus equivalent-analyte basis was not resolved, so version 3 excludes them via an explicit original-unit whitelist. The audit also found component IDs were derived from global row indices: this did not cause within-run leakage, but changing an unrelated endpoint could alter split randomization. Component IDs now hash sorted chemical families and are index invariant. All version 2 runs/reports remain archived locally and are superseded. Version 3 reruns every model with unchanged hyperparameter grids and declared seeds; no outcomes are used to select new cohorts or favorable splits. Original Dataset A/B files remain unchanged.

## Checkpoint 4 — completed version 3 experiments and analysis

The corrected cohort contains 413 full-InChIKey compounds (405 conservative chemical families), 438 original endpoint records, 437 assay-specific task rows and 105 distinct source publication identifiers. Endpoint cohorts are 169 clearance compounds / 172 rows / 64 publications and 259 Cmax compounds / 265 rows / 47 publications. The original A/B numerical datasets are unchanged.

Executed 190 primary runs and nine focused follow-ups: five training-size/feature ablations and four additional ensemble members. All four stage failure lists are empty. Version 3 experiment events run from 2026-10-30 14:41:42 UTC to 14:43:38 UTC; summed fit/tuning wall time is 117.1 seconds. This excludes environment setup, scientific auditing, superseded pilots, analysis and validation; no 72-hour runtime is claimed. Ten model methods per admissible task/split case were evaluated, with validation-only hyperparameter selection. Two requested Cmax split cases failed minimum-size constraints and were not replaced; three identical clearance study/scaffold-study cases were executed once.

The strict clearance training-mean control outperforms molecular models on average; Cmax errors increase under combined publication/scaffold exclusion. Two of five exploratory GINE–XGBoost component intervals exclude zero in favor of GINE, but the remaining intervals include zero, split tests overlap, and no general GNN advantage is established. All models have negative mean strict-split log-space R². Broad conformal-style intervals have high observed coverage without demonstrating useful precision or exchangeability under shift. AUC/volume/half-life training, cross-source generalization and Dataset A→B transfer remain unsupported and were not forced.

Validation: 56 tests passed; all 199 runs independently reconciled for saved metrics, split memberships and validation selection (maximum numerical discrepancy below 3e-12). Original cached cross-dose pipeline replay reproduced 70 artifacts; Dataset A replay reproduced 14 artifacts and preserved 44 checked Dataset B artifacts. Fresh XGBoost/GINE refits on the declared clearance study/scaffold seed reproduce predictions byte-for-byte. Reproduction checks also rebuild ML eligibility/splits and render figures/reports; exact counts and hashes are in `reports/ml/reproducibility.json`. Library warnings concern openpyxl extension handling and Python 3.14/PyG deprecations; no tests fail.

Publication review: only code, configurations, aggregate scientific results, source-attributed CC BY-SA figures and validation documentation are included. Raw/processed datasets, per-record predictions, chemical structures and model checkpoints remain ignored. The dedicated experiment branch may be pushed as authorized; no merge into main.
