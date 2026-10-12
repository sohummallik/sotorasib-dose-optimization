# Sotorasib exposure and acid-reducing medicines

Proton pump inhibitor (PPI) use is one of the covariates in the published model. Comparing that result with the prescribing information shows why a PK estimate needs its study context before it can inform medication use.

## What the model and interaction studies show

The [published PPI coefficient](../data/poppk_parameters_final.csv) is −0.178. The implementation applies it to relative bioavailability, giving a multiplier of 1 − 0.178 = 0.822. It holds the other parameters constant, including apparent clearance. The [saved results](../results/pk_covariate_checks.csv) therefore show about 18% lower AUC and Cmax in the typical subject at 960 mg once daily. This is a reproduction of a published covariate association, not a new interaction study. [Nagase et al., Table II and Figure 3](https://doi.org/10.1208/s12248-024-01013-6).

All values below are exposure ratios with versus without the acid reducer. The label ratios are calculated as 1 minus the reported fractional decrease.

| Evidence and setting | AUC ratio | Cmax ratio |
|---|---:|---:|
| Published popPK PPI contrast at 960 mg once daily | 0.820 | 0.828 |
| This repository's reproduction, day 30 | 0.822 | 0.822 |
| Label: repeated omeprazole with single-dose sotorasib, fasted | 0.580 | 0.430 |
| Label: repeated omeprazole with single-dose sotorasib, fed | 0.430 | 0.350 |

The interaction results come from the [U.S. prescribing information, Section 12.3](https://www.pi.amgen.com/united_states/lumakras/lumakras_pi_hcp_english.pdf), revised January 2025 and checked October 11, 2026. Sotorasib has pH-dependent solubility. The label also states that its exposure-response relationships are unknown. [Sections 7.1 and 12.2](https://www.pi.amgen.com/united_states/lumakras/lumakras_pi_hcp_english.pdf).

These estimates answer different questions. A covariate effect from a pooled model is not interchangeable with an omeprazole comparison under specified meal conditions. Single-dose and repeated-dose sotorasib also differ because its PK changes over time. The label results are therefore not additional validation targets for the existing model checks. AUC after an oral dose alone cannot separate bioavailability from clearance; assigning the PPI effect to bioavailability is part of this model's structure.

## Implication for medication review

The label advises avoiding PPIs and H2 receptor antagonists. When an acid reducer cannot be avoided, the spacing instruction is to give sotorasib at least four hours before or ten hours after a **local antacid**; it does not establish a spacing workaround for PPIs. [Sections 2.4 and 7.1](https://www.pi.amgen.com/united_states/lumakras/lumakras_pi_hcp_english.pdf).

The smaller modeled PPI effect does not override that advice. Neither evidence source supports calculating a compensating sotorasib dose from the exposure ratio or assuming that food cancels the interaction.

## Proposed follow-up analysis

**Question:** In an existing NSCLC monotherapy dataset, which treatment and patient factors explain measured sotorasib exposure, and how does early exposure relate to later tolerability?

This is a proposed research question. No patient-level dataset was available for this project, and the analysis below has not been performed. It would use existing records rather than deliberately introduce a label-discouraged drug combination.

| Records needed | Why they matter |
|---|---|
| Actual dose and administration times; missed doses, interruptions, and reductions | Assigned starting dose does not describe every dose received. |
| Concentrations, sample times, treatment day, and assay quality | A charted dose or one untimed concentration cannot provide interval AUC. |
| Meals and acid-reducer type, dose, dates, and indication | Exposure conditions and the reason for treatment may differ between patients. |
| Baseline albumin, ECOG status, tumor burden, liver tests; dated adverse events and dose changes | Disease severity may affect both exposure and outcomes; event timing determines which exposure could precede toxicity. |

First, assess whether the concentration sampling and dosing histories support an exposure analysis. AUC would require adequate timed sampling or a separately evaluated patient-level PK model. For a model-based approach, compare predictions with observed concentrations across treatment days and acid-reducer groups, and examine residual patterns before using exposure estimates in an outcome analysis. The current typical-subject checks do not establish that suitability.

Prespecify an early exposure window and a subsequent tolerability endpoint before examining their association. Report patients missing that window and any earlier toxicity. Otherwise, toxicity that led to a dose reduction could misleadingly appear associated with the lower exposure measured afterward.

With enough data, consider baseline prognosis and acid-reducer indication alongside exposure. If observations or events are sparse, keep the analysis descriptive. Even an adjusted association would leave residual confounding and would not establish that changing exposure improves outcomes or define a therapeutic target.
