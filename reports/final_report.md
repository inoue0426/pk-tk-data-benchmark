# Verified expanded dataset — 2026-10-28

The primary subset now contains **26 comparable multi-dose series, 12 compound identities and 10 independent source studies**. The previous run had 9 series, 7 compounds and 5 studies. Species, sex, route and repeated-dose day strata within a study are correlated; they are not independent studies.

## Counts after each stage

| Scope after stage | Validated compounds | Independent studies | Series | AUC | Cmax | Both eligible |
|---|---:|---:|---:|---:|---:|---:|
| identity_resolution | 7 | 5 | 9 | 7 | 4 | 2 |
| composite_audit | 10 | 8 | 23 | 19 | 6 | 2 |
| additional_publications | 12 | 10 | 26 | 21 | 9 | 4 |

These are reproducible scope snapshots using the final identity registry, not model results. Source studies are identified by publication DOI for PMC and experiment number for CEBS; a CEBS source and its CvTdb copy are one study.

- Collected parent identifiers: **254**; resolved unique full-InChIKey structures: **227**; unresolved identifiers: **26**. A fully resolved unique-molecule total cannot be claimed. Analyte identifiers including metabolites: 332.
- Collected document-level study proxies: **135**; this is not the number of independent experiments. PMC publications collected: 8.
- Compounds with ≥2 comparable doses: **12**; with ≥3: **11**.
- AUC-eligible series: **21**; Cmax-eligible: **9**; both eligible: **4**. 9 primary series report both endpoints, including AUC definitions not accepted for comparison.
- Primary dose groups: **78/104 (75.0%)**. Verified comparable parent-identifier fraction: **4.72%**. These denominators are dose groups and parent identifiers, not raw observations.
- Collected species: dog, hamster, monkey, mouse, rabbit, rat; routes: intravenous, oral, subcutaneous. Primary subset: mouse, rat; intravenous, oral; plasma.

## Validated series

| Compound | Species / sex | Route | Doses | Eligible endpoints | Study |
|---|---|---|---:|---|---|
| DA-6886 | rat / male | intravenous | 3 | AUC | doi:10.3390/pharmaceutics14040702 |
| DA-6886 | rat / male | oral | 3 | AUC, Cmax | doi:10.3390/pharmaceutics14040702 |
| Isoeugenol | mouse / female | oral | 3 | Cmax | CEBS:K10262 |
| Isoeugenol | mouse / male | oral | 3 | Cmax | CEBS:K10262 |
| lumefantrine | rat / male | oral | 3 | AUC, Cmax | doi:10.1186/1475-2875-10-293 |
| lorlatinib | rat / male | oral | 3 | AUC, Cmax | doi:10.2147/dddt.s565228 |
| saikosaponin A | rat / male | intravenous | 2 | AUC | doi:10.3390/pharmaceutics10030133 |
| imperatorin | rat / male | intravenous | 3 | AUC | doi:10.3390/pharmaceutics10030133 |
| Benzophenone | mouse / female | oral | 3 | AUC | CEBS:S0592 |
| Benzophenone | mouse / male | oral | 3 | AUC | CEBS:S0592 |
| Benzophenone | rat / female | oral | 3 | AUC | CEBS:S0592 |
| Benzophenone | rat / male | oral | 3 | AUC | CEBS:S0592 |
| paeonol | rat / male | intravenous | 3 | AUC | doi:10.3390/pharmaceutics10030133 |
| tofacitinib | rat / male | intravenous | 4 | AUC | doi:10.3390/pharmaceutics11070318 |
| Formamide | mouse / female | intravenous | 3 | AUC | CEBS:S0613 |
| Formamide | mouse / male | intravenous | 3 | AUC | CEBS:S0613 |
| Formamide | mouse / female | oral | 3 | AUC | CEBS:S0613 |
| Formamide | mouse / male | oral | 3 | AUC | CEBS:S0613 |
| Formamide | rat / female | intravenous | 3 | AUC | CEBS:S0613 |
| Formamide | rat / male | intravenous | 3 | AUC | CEBS:S0613 |
| Formamide | rat / female | oral | 3 | AUC | CEBS:S0613 |
| Formamide | rat / male | oral | 3 | AUC | CEBS:S0613 |
| methyl eugenol | mouse / male | oral | 3 | AUC, Cmax | doi:10.3390/cimb46100673 |
| methyl eugenol | mouse / male | oral | 3 | Cmax | doi:10.3390/cimb46100673 |
| HY-071085 | rat / female | oral | 3 | Cmax | doi:10.3390/ph15070801 |
| HY-071085 | rat / male | oral | 3 | Cmax | doi:10.3390/ph15070801 |

## Source contributions and low-yield approaches

CEBS/NTP group spreadsheets add 14 validated series from three already collected compounds and three source studies: benzophenone (4 AUC series), formamide (8 AUC series) and isoeugenol (2 mouse Cmax series). These validate existing CvTdb candidates; they are not new independent copies of those observations. Formamide's predicted Cmax values remain in original source cells and are not presented as observed Cmax.

Two additional PMC publications add three series: lorlatinib (rat, oral; AUC0-24h and Cmax) and methyl eugenol (mouse, oral; day 1 AUC0-inf/Cmax and day 28 Cmax). The two methyl eugenol time-course strata use the same animals and count as one study. Day 28 AUC0-t is preserved but not eligible without an explicit endpoint window. Lorlatinib is CC BY-NC 4.0; its local derived data retain the noncommercial restriction and are not committed. It can be excluded using the license flag.

The composite audit reviewed 18/112 candidates across three protocols: 14 supported group-endpoint comparisons, two had plasma/blood conflicts and two had 35-versus-70 mg/kg dose conflicts. The remaining 94 were not reviewed. This purposive sample is not an estimate of yield for all 112 candidates. No unreviewed group is promoted and no individual-animal trajectory is invented. See [audit decisions](composite_audit.csv).

The identity pass does not itself create experiments. Direct PubChem REST remains outside the allowed robots paths; PubChem CIDs are exact-structure UniChem cross-references. EPA DSSTox bulk matches require exact DTXSID/CAS and full source-SMILES/InChIKey agreement. Cached ChEBI evidence is checked for cross-source conflicts; remaining identifiers use ChEBI queries. Unsuccessful and conflicting candidates remain explicit. PK-DB previously yielded metadata but no eligible preclinical measurements; the DA-6886 supplement added no exposure series. The simmitecan XML endpoint failed after retries in this expansion and contributed no data. Additional search leads are not counted as collected studies.

## QC and remaining limitations

The previously unqueried 238 identifiers now have statuses: `{"ambiguous_multiple_structures": 3, "not_found": 1, "resolved_exact": 213, "unresolved_or_identifier_conflict": 21}`. Full registry identity status counts: `{"ambiguous_multiple_structures": 3, "not_found": 2, "resolved_exact": 228, "unresolved_or_identifier_conflict": 21}`. 58 registry entries have one PubChem CID and 11 primary compounds have resolved structures. Full stereochemical/isotopic keys are retained; salts and ambiguous mappings are not collapsed. CID lookups cover primary/candidate compounds and existing caches; other resolved compounds may have no optional CID query. See [EPA bulk manifest](dsstox_identity_manifest.json).

The pipeline retains all 35,669 original CvTdb observations, including 34,907 preclinical observations. The original audit found 17,020 singleton series out of 22,168; sparse sampling and censored observations explain much of the curve attrition. Endpoint imports preserve all original spreadsheet fields/cells. CEBS actual versus nominal dose, repeated model fits, predicted versus observed peaks, and source metadata errors are explicitly distinguished.

10 curated series remain outside the primary subset. Primary records still flag missing fasting information, unavailable AUC extrapolated fractions, composite/aggregate uncertainty and source exceptions. Different AUC definitions, study days and matrices are never pooled. Missing structures and identifier conflicts remain limitations. The 112-candidate universe was exploratory, not 112 established experimental series.

No complex models were trained. The cohort has only 12 compounds in 10 studies; several compounds share a publication. This is insufficient for meaningful compound-disjoint tuning or scaffold evaluation. Existing baseline results and splits concern only the original HY-071085 pilot. Expanded endpoints have not been used to claim new model performance.

See [inventory](data_inventory.md), [QC](qc_report.md), [validation](validation.md), [statistics](expansion_summary.json), [identity provenance](identity_source_manifest.json) and [unresolved issues](unresolved_issues.json).
