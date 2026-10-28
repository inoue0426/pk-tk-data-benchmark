# Verified expanded dataset — 2026-10-28

The curated primary subset grew from **2 to 9 comparable multi-dose series**, from **1 to 7 compound identifiers**, and from **1 to 5 independent publications**. Series within the same publication are correlated strata, not independent studies. This is a validated dataset expansion, not a claim of adequate ML diversity.

| Compound | Route | Doses | Eligible endpoints | Study |
|---|---|---:|---|---|
| DA-6886 (male) | intravenous | 3 | AUC | 10.3390/pharmaceutics14040702 |
| DA-6886 (male) | oral | 3 | AUC,Cmax | 10.3390/pharmaceutics14040702 |
| lumefantrine (male) | oral | 3 | AUC,Cmax | 10.1186/1475-2875-10-293 |
| saikosaponin A (male) | intravenous | 2 | AUC | 10.3390/pharmaceutics10030133 |
| imperatorin (male) | intravenous | 3 | AUC | 10.3390/pharmaceutics10030133 |
| paeonol (male) | intravenous | 3 | AUC | 10.3390/pharmaceutics10030133 |
| tofacitinib (male) | intravenous | 4 | AUC | 10.3390/pharmaceutics11070318 |
| HY-071085 (female) | oral | 3 | Cmax | 10.3390/ph15070801 |
| HY-071085 (male) | oral | 3 | Cmax | 10.3390/ph15070801 |

All counts are computed from downloaded bytes. The pipeline retains 35,669 CvTdb observations (34,907 preclinical) and 47 curated PMC dose groups from 6 publications. Across sources there are 253 preclinical parent-compound identifiers, 331 analyte identifiers including metabolites, and 133 document-level study proxies after exact DOI normalization. These are not fully resolved unique molecules or proven independent experiments.

- Compounds with ≥2 comparable doses: **7**; with ≥3: **6**.
- Primary series with reported AUC and Cmax: **4**. AUC eligible: **7**; Cmax eligible: **4**; both eligible: **2**. AUC definitions are kept separate; no assumed AUC0-t windows.
- Primary dose groups: **27/47 (57.4%)** of curated PMC groups. Verified comparable compound fraction: **2.77%** of collected parent identifiers. No mixed timepoint/group denominator is presented as an experiment fraction.
- Collected species: dog, hamster, monkey, mouse, rabbit, rat; routes: intravenous, oral, subcutaneous. Primary subset: rat, intravenous, oral, plasma.

## Why CvTdb yielded only ten candidates

17,020 of 22,168 source series IDs have one observation, often individual animals in NTP datasets. A three-timepoint-per-ID requirement misses potentially reconstructable composite designs. ND/NQ omissions recover only 11 candidate groups versus ten under the original gate. There are 112 possible groups if composite sampling is established; none is promoted without source protocol review. Relaxed metadata grouping yields 9 groups rather than ten (some groups merge), so merely loosening grouping does not solve the problem. The parser did not discard observations. See [quantitative audit](cvtdb_attrition.md).

## Identity validation and collection

All 253 parent-identifier registry entries retain their source names/CAS values. 13 distinct structures resolve through ChEBI exact CAS/name evidence and RDKit agreement with the full source InChIKey; 13 entries have a single PubChem CID cross-reference through UniChem. Multiple structures, CAS/name conflicts, salts, stereoisomers and failed lookups remain explicit alternatives. 6 primary compounds have resolved structures. No name-only fuzzy match or automatic parent/salt collapse is accepted.

The DA-6886 supplementary ZIP was downloaded and its PDF inspected: calibration/stability controls and chromatograms, not extra in vivo dose groups. The legacy PK-DB statistics and single-study metadata endpoints were accessible after reducing page size; the inspected record is human and closed-licence and was not imported. New PK-DB endpoints had TLS failures. PubChem's direct REST path and ChEBI's FTP listing were disallowed by robots; public ChEBI/UniChem APIs were used with their own documented terms. No restrictions were bypassed.

## Remaining limitations and benchmark feasibility

8 additional within-publication series remain outside the primary subset: oral vehicle composition changes, unreported oral vehicle, dose-volume changes/composite sampling, or suspect units. The original HY-071085 AUC terminal time and dog table remain unresolved. Some primary IV studies lack fasting details and AUC tail fractions; these are flagged, not filled in. Primary evidence includes group summaries and no new numerical individual curves.

No complex models were trained. The seven primary compounds and five publications are too few for reliable compound/scaffold-disjoint evaluation or tuning. Existing baseline results concern only the original HY-071085 pilot; they are not expanded-dataset results. Task A is structurally feasible for eligible endpoints, but no additional model development was prioritized. Tasks B–D remain unsupported for meaningful generalization/classification. Identity and protocol curation, especially NTP composite profiles and new independent publications, remain the next priorities.

See [source inventory](data_inventory.md), [QC](qc_report.md), [validation](validation.md), [machine-readable statistics](expansion_summary.json), [identity provenance](identity_source_manifest.json), and [unresolved issues](unresolved_issues.json).
