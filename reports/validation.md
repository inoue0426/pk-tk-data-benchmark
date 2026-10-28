# Executed validation — 2026-10-28

- Online expansion collection → identity resolution → harmonization → QC: passed.
- Full offline replay (`python -m pktk.pipeline --offline`) from checksum-verified source cache: passed, including the existing pilot evaluation and expanded curation.
- All **15 checked CSV/JSON data and summary artifacts** were byte-identical between the online expansion and full offline replay. Timestamped collection/progress logs are intentionally excluded.
- `python -m pytest -q`: **18 passed, 0 failed, 0 skipped**, including actual downloaded-data integration checks.
- `python -m compileall -q src` and `git diff --check`: passed.

Tests cover unit dimensionality, censored inputs, structure and CAS conflicts, ambiguous candidates, InChI fallback for invalid source SMILES, cache integrity and POST request identity, robots denial, pagination cycles, source XML extraction, unit conversion, grouping across publications/sex/matrices, and the existing pilot extrapolation splits. Expanded Task A episodes have not been generated or evaluated; existing baseline results remain limited to HY-071085. Compound/scaffold splits are withheld because diversity is insufficient.

These checks validate reproducibility and implemented data constraints, not statistical external validity or all CvTdb sampling protocols. Identity queries cover 15 prioritized identifiers; 238 others are explicitly not queried. Eight curated series remain outside primary eligibility.
