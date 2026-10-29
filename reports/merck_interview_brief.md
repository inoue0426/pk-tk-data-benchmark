# Interview brief: preclinical cross-dose exposure prediction

Personal preparation for discussing an independent public-data investigation. No Merck internal data were accessed. Adapt these talking points to your own voice; distinguish measured results from proposed research.

## A. 30-second introduction

“I wanted to understand what makes cross-dose exposure prediction difficult before discussing model choices. I independently collected public preclinical PK data and built a small, protocol-reviewed benchmark. The main lesson was that having many molecular PK records is different from having comparable dose experiments. I tested simple scaling methods first to see where they work and where additional molecular or mechanistic information might help.”

## B. 60-second research overview

“I built two separate datasets. The broader collection has 2,713 structure-resolved compounds with numerical preclinical PK endpoints, mainly from ChEMBL and measured-animal PKSmart data. The strict benchmark is much smaller: 59 comparable dose series across 21 compounds and 19 studies.

I held out doses and predicted them using only allowed initial measurements, separating interpolation from extrapolation. I also grouped related compounds and studies so repeated dose pairs would not inflate the effective sample size.

For two-support low-to-high prediction, linear scaling had a log MAE of about 0.53 for AUC and 0.60 for Cmax. Power law was worse for AUC and better on Cmax log error, but not on relative error. The uncertainty intervals included zero, with only 11 or 12 independent units. I see this as a preliminary baseline study, not evidence that one model generally wins. I did not train molecular ML.”

## C. Three research insights

1. **Dose-level evidence and comparability come first.** A compound can have AUC and Cmax in a database without any valid cross-dose comparison. I kept study, species, route, matrix, regimen and formulation differences explicit. The main limitation was independent experimental diversity, not the number of rows.
2. **Flexible models need enough information.** A two-point power law fits exactly, so a good support fit is not evidence of reliable extrapolation. It may adapt to a real trend or amplify measurement error. Linear scaling is a useful comparator; the current results do not establish a consistent advantage for power law.
3. **There is room for methodological innovation, with stronger validation.** I would investigate whether molecular representations provide useful priors beyond initial PK measurements, whether mechanistic constraints improve extrapolation, and whether uncertainty can guide the next dose experiment. Those are hypotheses to test on independent compounds, not outcomes of this project.

## D. Three questions for the interviewers

1. “How many independent compounds and studies have several comparable doses, and what proportion include concentration–time profiles and animal-level variability? How does that coverage vary by species and route?”
2. “Where is there the most room for methodological innovation: molecular representations, mechanistic hybrids, uncertainty estimation, or selecting the next informative experiment? What empirical baseline would a new method need to improve on?”
3. “Could predictions be tested prospectively in dose-escalation experiments, and how would success be judged? What opportunities are there to publish methodological findings while protecting confidential compound information?”

## E. Five likely technical questions

**1. What are AUC and Cmax, and why evaluate both?**

AUC integrates concentration over a specified time interval and summarizes exposure; Cmax is the peak concentration, often the highest sampled value in these sources. AUC0–t and AUC0–inf differ. Cmax is particularly sensitive to absorption and sampling schedule, so a model can perform differently on the two endpoints. I evaluate them separately and preserve matrix and endpoint definitions.

**2. What does dose proportionality mean?**

Under comparable conditions, exposure is proportional to dose: `Y = kD`. Dose-normalized exposure is then constant, and the log–log slope is one. In the benchmark, linear scaling estimates the normalization from permitted support measurements. Proportionality over observed doses does not prove it beyond that range.

**3. What might cause nonlinear exposure?**

Potential explanations include dose-dependent absorption or bioavailability, clearance, transport, binding, or formulation effects. Sampling differences can also create apparent trends. A nonlinear AUC or Cmax relationship alone does not identify which explanation applies, and does not prove saturation. I would examine profiles, protocols and variability before interpreting a mechanism.

**4. How does identifiability affect model choice?**

It asks whether the available measurements can determine the model parameters. Two distinct positive observations determine a power-law intercept and exponent, but leave no residual degrees of freedom to estimate within-series fit uncertainty. More complex curves can have poorly constrained parameters even when optimization converges. I withheld nonlinear predictions when support or numerical identifiability checks failed.

**5. What is compound-disjoint generalization?**

Training and test compounds must not share molecular identities, and linked studies should not leak across the split. Closely related representations and scaffolds need additional care. In a limited-observation task, a new test compound can have explicitly permitted initial PK measurements while remaining excluded from global training. That is different from zero-shot prediction from structure alone. The current pooled empirical evaluation excludes linked compound–study components, but no molecular ML or robust scaffold-disjoint performance was evaluated.

## Presentation notes

Lead with the two-page PDF and explain the A/B boundary before showing model results. Use `06_performance_by_compound_name` to show heterogeneous errors, and `08_proportionality_examples` to discuss observed trends. The latter uses independent y-axis ranges; the near-proportional panels are visually zoomed. Never call 513 episodes independent experiments, 2,713 compounds cross-dose validated, or exploratory numerical differences statistically established improvements. See [figure audit](figure_audit.md) for captions and source restrictions.
