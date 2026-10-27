# Verified pilot results

The complete collection → harmonization → QC → Task A baseline workflow was executed on CvTdb and an open-access PMC article. This release is a reproducible working subset, not a finished molecular benchmark.

| Measurement | Executed count |
|---|---:|
| CvTdb original concentration observations retained | 35,669 |
| Preclinical concentration observations | 34,907 |
| Preclinical compound/analyte identifiers including PMC | 324 |
| Preclinical parent-compound identifiers including PMC | 246 |
| Preclinical document-level study proxies including PMC | 132 |
| Preclinical source study IDs including PMC (dose-specific in CvTdb) | 807 |
| Usable source curves after computational QC | 151 |
| Provisional CvTdb cross-dose series / compounds | 10 / 7 |
| Provisional CvTdb compounds with at least three doses | 2 |
| Verified compounds with at least two comparable doses | 1 |
| Verified compounds with at least three comparable doses | 1 |
| Verified sex-specific series containing reported AUC and Cmax | 2 |
| Series eligible for AUC benchmarking | 0 |
| High-confidence Cmax groups / series | 6 / 2 |
| Task A episodes / predictions | 4 / 8 |

Species collected: dog, hamster, monkey, mouse, rabbit, rat. Routes: intravenous, oral. The accepted pilot is oral rat plasma only.

Six of 34,913 preclinical source records (0.0172%) are accepted Cmax groups. This is a source-record fraction mixing timepoints and endpoint groups, **not an independent experiment fraction**. One of 324 analyte identifiers (0.309%) has verified comparable doses. Counts are computed, not estimated, but unique molecular identities and independent study counts remain unverified: aliases, metabolites, missing structures and dose-specific study IDs matter.

The two sex-specific rat series have reported AUC and Cmax, but AUC terminal times are unresolved. Only Cmax was evaluated. With two support doses and two high-dose targets, linear scaling and power law are compared on identical observations; exact natural-log and relative-error metrics are in [benchmark feasibility](benchmark_feasibility.md). This one-compound sample cannot rank model families or establish generalization.

No compound/scaffold-disjoint split, molecular ML model, calibrated interval or nonlinear classification score is claimed. The pilot has no structures and insufficient independent compounds. Missing CvTdb strains/protocols, source censoring, heterogeneous sampling, incomplete cross-source identity resolution and the conflicting PMC dog table remain material limitations. See [unresolved issues](unresolved_issues.json), [QC](qc_report.md), [inventory](data_inventory.md) and [validation](validation.md).
