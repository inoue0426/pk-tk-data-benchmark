# Public preclinical PK/TK cross-dose benchmark

An auditable Python pipeline for collecting real public PK/TK data, preserving source records, identifying candidate dose series, and evaluating a conservative within-compound pilot. This is an initial research dataset, **not yet a molecular ML benchmark**.

The executed release uses CvTdb's original concentration–time table packaged in `invivoPKfit 2.0.2` and the CC BY article [Dong et al. (2022)](https://doi.org/10.3390/ph15070801). See [source inventory](reports/data_inventory.md), [actual counts](reports/collection_summary.md), [QC](reports/qc_report.md), [feasibility and baseline results](reports/benchmark_feasibility.md), and [dictionary](reports/data_dictionary.md).

## Reproduce

Python 3.11+ with RDKit wheels and internet access is required. Tested environment: Python 3.14 on macOS arm64; exact installed versions are recorded in `requirements-lock.txt`. R is **not** required.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-lock.txt
python -m pip install -e . --no-deps
pktk --root .
python -m pytest -q
```

For other platforms, `python -m pip install -e '.[test]'` resolves compatible dependencies from `pyproject.toml`; that is not an exact replication of the recorded environment. The first run downloads about 2 MB of scientific source files plus robots policies. Parsing and QC may take several minutes. Files are cached with SHA-256, URL, UTC retrieval timestamp and source version. The CvTdb checksum is pinned in `configs/sources.json`. Changed cached files fail integrity checks; do not silently replace snapshots.

```bash
pktk --root . --offline
```

Offline reruns require the cache **and** its JSON manifests. Reports and generated CSVs are overwritten deterministically from source data; progress logs append events. Network errors, denied robots rules, authentication and unreviewed redirects fail explicitly. The collector retries transient HTTP failures and respects Retry-After. It does not bypass access restrictions. A generic cached `results`/`next` paginator is available, but selected sources use whole-file downloads.

## Outputs and scope

- `data/raw/`: immutable downloaded bytes, per-file manifests and collection success/failure log.
- `data/interim/`: decoded original CvTdb table and upstream package metadata.
- `data/processed/observations.csv`: every original CvTdb observation, provenance, original fields and QC flags, including out-of-scope human/metabolite rows.
- `data/processed/curve_qc.csv`: every source series ID and its derivation eligibility. No assumption that a series ID represents an independent animal.
- `data/processed/candidate_series.csv`: CvTdb groups sharing publication, compound and available conditions, including one-dose groups marked insufficient. These remain **unverified** without strain/shared-protocol adjudication.
- `data/processed/reported_endpoints.csv` and `reported_series.csv`: article-specific rat single oral dose groups, parsed with structural checks against study-design tables.
- `data/processed/high_confidence_cmax.csv`: six group means for one compound, two sex-specific three-dose series. “High confidence” describes extracted within-protocol Cmax comparability, not statistical precision or molecular identity validation.
- `results/`: Task A support/target manifests, predictions, log/relative error metrics, descriptive nonlinearity and feasibility status for Tasks B–D.
- `reports/summary.json`, `qc_counts.json`, `source_manifest.json`, `progress.jsonl`, `unresolved_issues.json`: machine-readable counts, provenance and outstanding work.

No record is silently removed. Derived partial AUC integrates only observed times. AUC definitions, windows, routes and matrices stay separate. The paper's rat AUC0-tn is retained but **not benchmarked** because per-group terminal times are unverified. Its dog table conflicts with its prose and is not imported into the endpoint dataset; the raw XML preserves it. No figures are digitized here. None of the selected sources provides a validated SMILES mapping in this workflow; RDKit standardization is implemented and tested, but structures are left missing rather than guessed.

Task A uses low doses as support and strictly higher doses as targets; same-compound support is intentional. Linear and power-law methods are compared with identical two-dose support where possible. Episodes reuse targets, so prediction counts are not independent sample sizes. Only one compound is sufficiently annotated: compound/scaffold splits, nonlinear classifiers, molecular ML and calibrated prediction intervals are withheld. No model tuning is performed.

## Licenses and attribution

Repository code uses the existing MIT license. Source licenses do not change: the EPA CvTdb release advertises CC0, the CRAN package is GPL-3, and Dong et al.'s article is CC BY 4.0. Cite both the database and underlying original publications when using measurements. Raw, interim, processed and result files are ignored by Git to prevent accidental redistribution; aggregate reports and source manifests are committed. See the inventory for restrictions on sources not downloaded. Nothing is pushed to a remote by this workflow.
