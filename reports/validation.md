# Reproducibility validation

Complete offline collection → harmonization → QC replay executed. 70/70 selected data/report artifacts are byte-identical to the preceding run. Differences: [].

Tests: 45 passed; pytest exit code 0. Compilation and git diff whitespace checks passed. See validation.json for hashes. Expanded empirical benchmarks and all figures were regenerated; no complex molecular ML was trained.

The replay uses locally cached source bytes and identity API snapshots. It does not prove future endpoint availability or bitwise equivalence of fresh external API responses. Run `PYTHONPATH=src python scripts/verify_reproducibility.py` after an online full/expansion run.
