# Public preclinical PK/TK cross-dose benchmark

A reproducible collection, identity-resolution and QC pipeline for dose–exposure research. The current expansion prioritizes independent publications and experimentally comparable dose series; it is **not yet a sufficiently diverse molecular ML benchmark**.

See the [verified final report](reports/final_report.md), [machine-readable counts](reports/summary.json), [CvTdb attrition investigation](reports/cvtdb_attrition.md), [source inventory](reports/data_inventory.md), [QC report](reports/qc_report.md), and [benchmark feasibility](reports/benchmark_feasibility.md).

## Reproduce the complete workflow

Python 3.11+ is required; tested with Python 3.14 on macOS arm64. No R installation is needed.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-lock.txt
python -m pip install -e . --no-deps
pktk --root .
python -m pytest -q
```

For other platforms, `python -m pip install -e '.[test]'` resolves compatible dependencies but is not an exact recreation of the recorded environment.

The full command runs the original CvTdb collection/harmonization/QC, the original simple HY-071085 baseline pilot, and then expanded PMC collection, chemical identity checks and current reports. No complex ML models are trained. The first online identity pass takes several minutes because the public APIs' crawl delays are respected. All scientific downloads, lookup results and negative identity searches are cached with provenance and checksums. Source/table hashes pin the curated exposure inputs; chemical API responses are preserved as the local snapshot.

```bash
# Complete offline replay, including identity checks from raw lookup caches:
pktk --root . --offline

# Update only expansion outputs after the original pipeline already ran:
python -m pktk.expansion --root . --offline

# Explicitly reproduce only the historical initial pilot:
pktk --root . --initial-only --offline
```

An expansion replay requires the cached original pipeline outputs as inputs. `--reuse-identities` additionally reuses the processed identity registry and is convenient during report development; use the default raw-cache replay for validation. Raw and generated datasets are ignored by Git. Download manifests and aggregate reports are committed; nothing is pushed automatically.

## Data and eligibility

- **CvTdb:** all 35,669 original observations from `invivoPKfit 2.0.2` are retained with original fields, source IDs, units and QC flags. The audit explains why only ten strict publication-matched candidate series emerged: most source series are single-timepoint animal records, and missing/censored concentrations block many others. Potential composite profiles are retained separately and never promoted without protocol evidence.
- **PMC studies:** the original HY-071085 rat tables plus reviewed tofacitinib, DA-6886, lumefantrine, DA-9805 pure-component and memantine tables. Exact article/table hashes, selected columns and evidence fragments are recorded in `configs/curation/pmc_expansion.json`. Exposures are parsed from original cells. Formulation changes, ambiguous units and incomplete protocol details are flagged.
- **Supplementary data:** the DA-6886 supplementary ZIP/PDF is downloaded and inspected. It contains calibration, stability and assay controls; these are not counted as additional animal PK experiments.
- **PK-DB:** anonymous legacy statistics and a small public metadata page were accessible. The inspected study is human and closed-licence; no measurements are imported. Newer endpoints failed TLS. Access/terms are documented without bypassing restrictions.
- **Identities:** ChEBI provides exact CAS/name evidence, SMILES and InChIKey; UniChem supplies exact-structure PubChem CID cross-references. RDKit validates structures, preserving salts and stereochemistry. Conflicts and alternatives remain explicit. A source InChI fallback is permitted only when it regenerates the full expected InChIKey.

The default identity pass prioritizes the eight curated PMC compounds and seven strict CvTdb candidate compounds. `configs/sources.json` sets `identity_priority_limit` to 15; increase it to cover additional collected identifiers. Remaining registry rows explicitly say `not_queried`. Never interpret a source-ID count as a confirmed unique-molecule count.

Study grouping retains publication, explicit protocol, species, strain, sex, route, matrix, regimen and formulation conditions. Molecular identity resolution cannot merge independent publications. `protocol_verified` and `primary_eligible` are separate: an identifiable dose study can still be unsuitable for extrapolation if its vehicle changes with dose. AUC eligibility is definition-specific; AUC0-inf, observed-window AUC and uncertain AUC0-tn are not interchangeable.

## Generated artifacts

| Location | Contents |
|---|---|
| `data/raw/` | Original scientific files, API JSON, manifests, successful/failed request log |
| `data/interim/` | Original CvTdb export and extracted supplement text |
| `data/processed/observations.csv` | Every CvTdb observation and source/QC provenance |
| `data/processed/cvtdb_profile_audit.csv` | Strict versus censored-profile recovery audit |
| `data/processed/cvtdb_composite_candidates.csv` | Dose-study cells requiring destructive-sampling/protocol adjudication |
| `data/processed/reported_endpoints_expanded.csv` | All curated PMC dose groups, including flagged comparisons |
| `data/processed/crossdose_series_expanded.csv` | Comparable and excluded series, dose ranges, endpoint eligibility, metadata flags |
| `data/processed/high_confidence_endpoints.csv` | Current primary comparable subset; eligibility remains endpoint-specific |
| `data/processed/compound_identity_registry.csv` | Resolved identities, unresolved/ambiguous mappings and unqueried source IDs |
| `data/processed/compound_identity_candidates.json` | Accepted/rejected alternatives and exact lookup provenance |
| `results/` | Historical HY-071085 Task A splits and simple baseline outputs only |
| `reports/` | Current counts, source/identity manifests, QC, feasibility, validation and unresolved issues |

The initial `reported_endpoints.csv`, `high_confidence_cmax.csv`, and `reported_series.csv` remain historical pilot artifacts. Current expanded outputs have the names shown above. Existing baseline metrics are not results for the expanded dataset. Compound/scaffold splits, nonlinear classifiers, prediction intervals and molecular ML remain withheld until independent study and compound diversity is adequate.

## Licenses and attribution

Code uses the existing MIT license. Source terms remain separate: EPA CvTdb original release CC0; packaged CRAN distributions GPL-3; selected PMC articles CC BY 2.0/4.0; ChEBI CC BY 4.0. UniChem references retain underlying source rights. Cite original studies and database publications. See [the inventory](reports/data_inventory.md) for source URLs, citations and restrictions. Do not apply the repository code license to downloaded datasets.
