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
