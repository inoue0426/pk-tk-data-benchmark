# Preclinical dose–exposure prediction — executive summary

**Objective.** Determine when initial in vivo PK measurements support prediction of exposure at unobserved doses, and what evidence would justify more flexible empirical or molecular models. This independent exploration uses public data only.

**Two distinct datasets.** Dataset A contains **2,713 structure-resolved compounds and 6,459 eligible numerical preclinical PK records**; 5,133 endpoint records were excluded but retained for audit. Endpoint coverage is 995 compounds with AUC, 994 with Cmax, 355 with both, 1,182 with clearance, 427 with volume and 1,250 with half-life. Sources are a fixed ChEMBL 37 tranche, measured-animal PKSmart data and structure-resolved strict endpoints. Its 627 publication/study proxies do not establish 627 independent experiments. **Dataset A is not a validated cross-dose benchmark.**

Dataset B contains **59 protocol-reviewed multi-dose series, 21 compounds and 19 studies**: 42 series support AUC, 33 Cmax and 16 both. Nineteen compounds have at least three comparable doses. Rat, mouse and hamster studies cover oral/IV dosing and blood/plasma, kept as separate strata. There are 513 evaluation episodes but only **19 linked compound–study units** overall.

**Main benchmark findings.** Two-support low-to-high extrapolation, with natural-log errors and equal series/component weighting:

| Endpoint | Model | Log MAE | Log RMSE | Mean absolute relative error |
|---|---|---:|---:|---:|
| AUC | Linear | 0.525 | 0.700 | 38.1% |
| AUC | Power law | 0.583 | 0.808 | 61.9% |
| Cmax | Linear | 0.601 | 0.793 | 50.9% |
| Cmax | Power law | 0.529 | 0.716 | 60.1% |

Only **11 independent AUC units and 12 Cmax units** contribute. Power-minus-linear log-MAE differences are +0.058 (exploratory 95% component-bootstrap interval −0.176 to +0.308) and −0.072 (−0.391 to +0.231), respectively. Both intervals include zero; neither superiority nor equivalence is established. No complex molecular ML was trained or evaluated.

**Three scientific insights.**

1. **Comparability is the bottleneck.** Thousands of structured endpoints do not substitute for independent, within-study dose series with compatible routes, matrices and exposure definitions.
2. **Flexibility is not reliably beneficial with sparse support.** Two points determine a power law but leave no residual information for within-series uncertainty. Linear scaling performs well for methyl eugenol and poorly in selected bromodichloromethane AUC and HY-071085 Cmax examples; these examples are post-hoc.
3. **Identifiability and generalization require different evidence.** Nonlinear trends do not identify a mechanism. Molecular priors, mechanistic hybrids and uncertainty-guided experiments are promising research hypotheses, not demonstrated improvements.

**Limitations and next steps.** Public-data selection bias, heterogeneous endpoint metadata, missing animal-level uncertainty and few independent compounds limit inference. Interpolation is easier than extrapolation and must remain separate. Robust scaffold-disjoint evaluation and calibrated prediction intervals are unsupported. Prioritize independent compounds with at least four comparable doses, defined AUC windows, observed Cmax and concentration–time profiles; then freeze compound/study-disjoint evaluations against matched-support linear scaling.

Counts and metrics were recalculated from saved data and predictions. Details: [final scientific report](final_scientific_report.md), [verification](final_verification.json), [two-page PDF](pk_tk_research_brief.pdf).
