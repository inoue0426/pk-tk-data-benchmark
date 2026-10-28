# Executed validation — 2026-10-28, second expansion

- Online collection → full identity assessment → harmonization → QC: passed.
- Complete offline replay (`python -m pktk.pipeline --offline`): passed, including the historical pilot and the current expanded cohort.
- **22 data/report CSV and JSON artifacts** were byte-identical between final online and offline execution. Checksums and the exact artifact list are in [validation.json](validation.json). Timestamped progress/download logs are excluded.
- `python -m pytest -q`: **23 passed, 0 failed, 0 skipped**.
- `python -m compileall -q src` and `git diff --check`: passed.

Tests include actual cached XML/XLSX integration, original-cell values and unit conversions, source license gates, CAS/DTXSID agreement, full structure-key conflicts, invalid-SMILES fallback, duplicate registry evidence, offline cache failure, study/sex/matrix separation, source dose versus nominal dose, predicted versus observed Cmax, repeated-day grouping, and the historical pilot's extrapolation/leakage checks. Three openpyxl warnings concern unsupported spreadsheet extensions; workbooks are read only and original bytes are preserved.

The full registry contains 254 identifiers. All 238 previously unqueried identifiers were assessed: 213 resolved, 21 unresolved/identifier conflicts, three multi-structure ambiguities and one not found. Across the full registry, 228 identifiers resolve to 227 distinct structures; 26 remain unresolved. No `not_queried` identity or lookup failure remains. Optional PubChem CIDs were not requested for every resolved structure.

The composite protocol audit covers 18/112 candidates, of which 14 support validated group endpoints. The other 94 are unreviewed, and the observed sample yield is not extrapolated. Original sources contain unresolved dose, matrix, terminal-time and nonstandard-InChI issues. These are scientific limitations, not test failures. Expanded endpoints were not used to claim new model performance; existing splits/results are pilot-only.
