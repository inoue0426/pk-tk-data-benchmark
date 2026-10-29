# Preclinical dose–exposure prediction: final scientific report

Verified finalization snapshot: 2026-10-29. This is a preliminary investigation using public data only. No proprietary pharmaceutical data, new data collection, or new model development was used during finalization. Counts were independently recalculated from local endpoint tables and losses from saved predictions; see [verification](final_verification.json).

## 1. Research motivation

Dose selection requires anticipating exposure at doses that have not yet been measured. The central question is: **When can preclinical drug exposure be reliably predicted at unobserved doses, and what additional information is needed beyond simple empirical scaling?** A large collection of molecular PK endpoints and a valid cross-dose benchmark serve different purposes. The former provides chemical coverage; the latter requires comparable experiments that permit an honest held-out-dose test.

## 2. Scientific questions

1. How much publicly accessible preclinical PK information has usable structures and numerical endpoints?
2. How much of that information establishes within-study dose comparability?
3. With one or two initial dose measurements, does adapting a dose–exposure exponent improve extrapolation over proportional scaling?
4. Which apparent nonlinearities are identifiable, and which require additional doses, concentration–time profiles or experimental uncertainty?

The first three questions have descriptive answers in this cohort. Mechanism identification and the value of molecular ML remain future research questions.

## 3. Public data sources and provenance

Dataset A combines a fixed ChEMBL 37 rat ADME tranche, the **measured animal Table S2** of the [PKSmart supplement](https://doi.org/10.1186/s13321-025-01066-5), and structure-resolved endpoints inherited from Dataset B. PKSmart model predictions and human measurements are excluded. ChEMBL contributed 4,885 eligible records / 2,337 compounds; PKSmart contributed 1,369 / 371; inherited strict endpoints contributed 205 / 19. Source compound counts overlap and must not be summed.

Dataset B contains 47 series from 12 CEBS source experiments and 12 series from seven reviewed publications. Original citations, download addresses, retrieval metadata, pinned raw hashes and cell selections are retained in the [source inventory](data_inventory.md), collection manifests and `configs/curation/`. CvTdb provided concentration–time observations and discovery candidates; protocol-reviewed source endpoints are not counted again as independent replicas of CvTdb copies.

ChEMBL data retain CC BY-SA 3.0 terms; PKSmart retains CC BY 4.0. Selected PMC studies have article-specific terms, including lorlatinib CC BY-NC 4.0. CEBS attribution and source exceptions remain explicit. The repository MIT license applies to code, not all underlying data. Raw and harmonized numerical datasets remain excluded from Git. This report summarizes existing findings; public availability is not a blanket redistribution license.

## 4. Data collection and harmonization

Cached downloads, source record IDs, original values/units, dataset versions and hashes provide replayable provenance. RDKit validates molecular structures; standardized full InChIKeys deduplicate Dataset A while preserving stereochemistry. Connectivity groups remain available to prevent related representations crossing future evaluation partitions. Ambiguous identities are retained rather than forcibly resolved.

Endpoint harmonization preserves AUC integration windows, dose normalization, matrix and administration route. AUC0–t and AUC0–inf are not interchangeable. Dose-normalized exposure does not establish an administered dose. Group endpoints remain separate from original concentration–time observations and individual animals. No missing exposure is imputed.

Dataset A eligibility screens numerical values, in vivo evidence, species, endpoint descriptions, units and structures. It does not constitute full experimental protocol review. Dataset B additionally requires source-supported within-study comparability, retaining species, strain, sex, route, matrix, regimen and vehicle strata. Matching molecular identity alone cannot merge studies. The composite audit assigned dispositions to 112 candidates; 15 still lack obtained/reconciled full protocols. A disposition is not validation.

## 5. Dataset A: molecular PK coverage

| Verified quantity | Count |
|---|---:|
| Structure-resolved compounds (full standardized InChIKey) | 2,713 |
| Connectivity groups | 2,647 |
| Collected endpoint records | 11,592 |
| Eligible numerical preclinical endpoint records | 6,459 |
| Excluded endpoint records, retained for audit | 5,133 |
| Compounds with AUC | 995 |
| Compounds with Cmax | 994 |
| Compounds with both AUC and Cmax somewhere in the dataset | 355 |
| Compounds with clearance | 1,182 |
| Compounds with volume of distribution | 427 |
| Compounds with half-life | 1,250 |
| Known source-publication/study proxies | 627 |

Endpoint coverage overlaps. Both-endpoint coverage does not establish paired endpoints under the same conditions. Species are rat, dog, monkey, mouse and hamster; reported routes include IV, oral, subcutaneous, intraperitoneal and intraduodenal, with missing route metadata retained. There are 24 explicitly dose-normalized records. PKSmart animal values are source aggregates, often lacking original dose and study identifiers.

The 627 publication/study proxies are **not 627 verified independent experiments**; the total number of independent experimental studies in A is unknown. A inherits 19 structure-resolved cross-dose compounds (17 with at least three doses) and 17 protocol-verified studies from B. It does not add validated cross-dose series. Eleven endpoint-specific candidate groups across seven publications remain unverified. The 500-compound structure/PK coverage target was exceeded; the cross-dose evaluation problem remains small.

ChEMBL retrieval used the first 10,000 activity IDs from a specified rat ADME query, not a random or exhaustive sample. Exclusion flags overlap, so their counts cannot be summed to obtain 5,133. See [molecular QC](molecular_qc.md) and [coverage report](molecular_pk_expansion.md).

## 6. Dataset B: strict cross-dose benchmark

| Verified quantity | Count |
|---|---:|
| Comparable multi-dose series | 59 |
| Compounds with at least two comparable doses | 21 |
| Compounds with at least three comparable doses | 19 |
| Compounds with resolved structures | 19 |
| Independent source studies | 19 |
| AUC-eligible series | 42 |
| Cmax-eligible series | 33 |
| Series eligible for both endpoints | 16 |
| Endpoint-specific series | 75 |
| Dose/endpoint observations (group endpoints) | 223 |
| Prediction episodes | 513 |
| Successful predictions across available empirical models | 1,369 |
| Linked compound–study components | 19 |

B covers rat, mouse and hamster; oral and IV administration; plasma and whole blood. AUC is available for 15 compounds and Cmax for 12; six compounds have both somewhere in B. The 59 protocol series include correlated sex, formulation or regimen strata. The 19 dependency components are the operational upper bound on independent evaluation units, not proof that every source study is biologically exchangeable. The 513 episodes are not 513 independent experiments.

## 7. Evaluation protocols and leakage prevention

**One-support tasks:** predict higher doses from the lowest measurement or lower doses from the highest. **Leave-one-dose-out:** for at least three eligible doses, hold one out and fit only the remaining observations. Extreme held-out doses define low-to-high or high-to-low extrapolation; interior targets define interpolation. Two-dose series do not duplicate the same tests as leave-one-dose-out episodes.

For low-to-high extrapolation the distance ratio is target / highest support dose; for high-to-low it is lowest support / target. Both are at least one. Interpolation is separately labeled and does not represent extrapolation distance.

Saved membership records and tests enforce held-out record separation, dose direction and compatible endpoint strata. Changing a held-out exposure cannot change its prediction. Pooled empirical reference slopes exclude the entire linked compound–study component, including co-studied compounds and resolved molecular aliases. Local measurements for the target compound remain permitted: this is limited-observation extrapolation, not zero-shot molecular prediction.

Losses are averaged first within endpoint series, then across series within a dependency component, then equally across components. RMSE is the square root of the resulting mean squared natural-log error. Mean absolute relative error is mean |prediction / observation − 1| under the same weighting. Repeated dose pairs do not increase biological replication.

## 8. Baseline extrapolation results

Linear scaling fits `log Ŷ(D*) = log D* + mean_support(log Y − log D)`. A power law fits `log Y = a + b log D` using at least two support doses. Both use exactly the permitted measurements; target exposure enters only scoring. Pooled and fixed 50:50 local/pooled exponent baselines were also evaluated when disjoint compatible reference components existed. These are empirical baselines, not molecular ML.

**Primary comparison: two-support low-to-high extrapolation.** All logarithms are natural.

| Endpoint | Model | Log MAE | Log RMSE | Mean absolute relative error | Independent units |
|---|---|---:|---:|---:|---:|
| AUC | Linear | 0.525 | 0.700 | 38.1% | 11 |
| AUC | Power law | 0.583 | 0.808 | 61.9% | 11 |
| Cmax | Linear | 0.601 | 0.793 | 50.9% | 12 |
| Cmax | Power law | 0.529 | 0.716 | 60.1% | 12 |

AUC uses 38 endpoint series across 12 compounds / 11 linked units; Cmax uses 33 series across 12 compounds / 12 units. These are paired comparisons on the same episodes. Independent recalculation from predictions agrees with the saved metrics and previously reported rounded values; no numerical correction was required.

The power-minus-linear log-MAE difference is **+0.058 for AUC** (exploratory 95% component-bootstrap interval **−0.176 to +0.308**) and **−0.072 for Cmax** (**−0.391 to +0.231**). Both intervals include zero. Power law did not consistently outperform linear scaling; the evidence establishes neither superiority nor equivalence. Cmax illustrates metric sensitivity: its power-law log errors are smaller, but its mean absolute relative error is larger.

Intervals use 2,000 whole-component resamples of fixed out-of-fold errors. They do not refit all models, quantify individual-animal uncertainty or guarantee population coverage. Comparisons are exploratory and unadjusted for multiplicity. Leave-one-component-out differences also cross zero in these primary comparisons.

Other directions reinforce the need to specify the task: AUC high-to-low log MAE is 0.483 / 0.625 (linear / power), versus 0.292 / 0.291 for interpolation; Cmax is 0.669 / 0.519 high-to-low and 0.296 / 0.257 for interpolation. Easier interpolation does not establish reliable extrapolation. Full compound-, study-, species-, route- and ratio-stratified results remain in `results/research/metrics.csv` and [the detailed analysis](scientific_analysis.md).

## 9. Nonlinearity and identifiability analysis

Dose-normalized exposure and all-dose log–log exponents describe observed trends. An exponent near one is compatible with proportionality over the measured range; it is not proof beyond that range. The descriptive normalized-ratio range 0.8–1.25 is a screening convention, not a significance test.

Methyl eugenol is a near-proportional example: selected male-mouse oral exponents are approximately 0.96 for AUC and 0.99 for Cmax. Its two-support linear low-to-high errors are small (compound-level log MAE approximately 0.043 / 0.027). Naphthalene AUC is another relatively favorable example (0.200). In contrast, bromodichloromethane AUC has compound-level linear log MAE about 1.131; one female-mouse oral series has a 17.47-fold high/low dose-normalized AUC ratio. HY-071085 Cmax has linear log MAE about 1.474 and a selected female-rat normalized ratio of 9.32. These are post-hoc illustrations, not prospectively validated classes of predictability.

With two support points, power-law intercept and exponent are exactly determined; there are no residual degrees of freedom for within-series uncertainty. Parameter estimation is possible while robust generalization remains unresolved. Nonlinear empirical curves were restricted to at least three support doses, positive parameters, an in-range curvature scale and a well-conditioned Jacobian. No saturating fit passed; only three superlinear predictions in one component passed. General nonlinear-model comparison is unsupported.

Missing animal-level variances and composite-sampling covariance prevent convincing statistical nonlinearity claims. Nonlinear exposure alone cannot distinguish changing bioavailability, clearance, absorption, formulation, sampling artifacts or other causes. It does not establish mechanistic saturation.

## 10. Data limitations and potential biases

Public data are selected, chemically uneven and dominated by available study summaries. ChEMBL A measurements vary in dose, route, matrix, endpoint definition and experimental completeness. Database eligibility is not full protocol validation. PKSmart aggregation loses original dose-level context. Dataset B has only 19 components overall, with fewer in endpoint-specific tests and two unresolved compound structures. Its source mix may not represent a contemporary pharmaceutical discovery portfolio.

AUC tail fractions, assay uncertainty, animal-level measurements and detailed concentration–time profiles are often absent. Whole blood and plasma, repeated-dose and single-dose regimens, and incompatible AUC windows cannot be pooled to increase sample size. Source aliases and unavailable protocols leave residual overlap uncertainty. Error versus dose ratio is confounded by compound/study composition; there is no established universal safe extrapolation ratio.

Final analysis rules were fixed after an audited pilot, not preregistered. Bootstrap precision is limited by small, purposively selected components. Independent prediction-interval calibration, robust scaffold-disjoint model selection and statistically validated nonlinear classification are unsupported. Complex molecular ML was neither trained nor evaluated; no molecular PK prediction performance is claimed.

## 11. Opportunities for future ML research

Dataset A could support feasibility work on narrowly defined endpoint/species/route cohorts after resolving publication dependencies and metadata gaps. Structure coverage alone does not demonstrate that a compound-disjoint prediction task is scientifically coherent. Future splits should exclude molecular identity/connectivity relatives and linked studies, with scaffold evaluation only when endpoint-specific diversity is adequate.

For cross-dose research, the highest-value additions would be independent compounds with at least four comparable doses, explicit AUC windows, observed Cmax, concentration–time profiles and animal-level uncertainty. This would separate fitting from validation and improve identifiability. A prospective compound-disjoint evaluation should compare molecular or hybrid models against linear scaling with the **same** initial PK measurements and support budget.

Hypotheses for future work include molecular representations informing a prior dose exponent, mechanistic constraints reducing implausible extrapolation, and uncertainty-guided selection of the next informative dose. None is demonstrated here. Mechanistic hybrids require information beyond aggregate nonlinear trends; reliable interval evaluation requires sufficient independent calibration and test compounds. Experimental confirmation would be needed before dose-selection use.

## 12. Conclusions

Public molecular PK coverage can be scaled to thousands of structures, but protocol-comparable dose series remain scarce. In the current small cohort, proportional scaling is a useful mandatory comparator and flexible exponents do not show consistent superiority. Some observed series are nearly proportional; others pose substantial extrapolation challenges. The data support these preliminary empirical findings, not a general prediction guarantee or a mechanistic explanation.

Reproduce counts with `python scripts/verify_final_results.py`; regenerate the two-page brief with `python scripts/render_research_brief.py`. Full cached replay and artifact checks are documented in the [README](../README.md). See [executive summary](executive_summary.md), [interview brief](merck_interview_brief.md), [figure audit](figure_audit.md), and [validation](validation.md). No new models or measurements were introduced during finalization.
