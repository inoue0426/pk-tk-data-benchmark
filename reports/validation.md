# Executed validation

- Online collection → harmonization → QC → Task A evaluation: passed.
- Offline full replay from checksum-verified source cache: passed.
- 10 existing generated CSV artifacts were byte-identical between online and offline executions. The newly added compound registry contains missing structures, not guessed mappings.
- `python -m pytest -q`: **11 passed, 0 failed, 0 skipped** (unit tests and actual downloaded-data integration checks).
- `python -m compileall -q src` and `git diff --check`: passed.
- Editable package installation and CLI help: checked.

Tests cover unit dimensionality, missing/censored inputs, RDKit identities, support-dose identifiability, extrapolation and matrix leakage, partial AUC bounds, duplicate timepoints, record preservation, cache corruption, robots denial, pagination cycles, actual XML extraction and generated splits. Statistical external validity and manual verification of all CvTdb protocols remain outside these checks.
