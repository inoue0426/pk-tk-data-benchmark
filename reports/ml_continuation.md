# Molecular PK checkpoint and continuation

Work is isolated on `experiments/molecular-gnn-benchmark`; do not merge automatically into main. Existing public Dataset A/B artifacts are immutable inputs, checked by `reports/ml/original_artifacts.json`. The first metadata-pooling pilot is superseded; never combine its results with the corrected `molecular-pk-v3-explicit-units-stable-groups` experiment.

## Reproduce or resume the completed scope

Use the README commands in order. `scripts/run_molecular_ml.py` reuses a completed run only when its task, split, configuration and training-code hashes match. Incomplete runs restart from their fixed seed; completed runs remain cached. The best validation epoch and preprocessing are saved locally. `results/ml/queue_*.json`, `events.jsonl`, `failures_*.json` and run manifests provide a resumable audit trail. No training process is promised to survive session termination.

Do not change the target gates or tune model grids after inspecting these test outcomes and call it the same benchmark. A metadata correction requires a new version and full replay; a methodological hypothesis informed by current results needs new held-out evaluation. Source data and large checkpoints remain ignored. Dependencies are pinned in `requirements-ml-lock.txt`.

## Deferred scientific tasks

1. Resolve original assay/protocol metadata and obtain additional independent publications before widening endpoint eligibility. In particular, recover AUC windows, terminal half-life definitions, Vss versus generic Vd, matrix and original studies behind compiled data. Dataset A coverage alone is insufficient.
2. Improve independent study diversity before interpreting strict Cmax rankings; two requested Cmax study/scaffold-study seed combinations failed minimum-size rules.
3. Seek enough independent calibration and test groups to evaluate useful interval widths and coverage under a declared target distribution. Current GINE variance is not a validated reliability score.
4. Test pretrained encoders, multitask learning or source shifts only after compatible labels, provenance and evaluation groups support them. No paid hardware is required by the present scope.
5. Dataset A→B transfer needs identity-disjoint representation training and sufficiently many target study components to compare against matched-support linear scaling. Current results do not justify a transfer claim.

These are research directions, not jobs silently queued or results already obtained. All feasible experiments in the frozen current scope were executed; requested analyses without defensible data were explicitly withheld.
