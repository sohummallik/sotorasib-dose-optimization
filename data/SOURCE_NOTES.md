# Sources and inputs

Evidence cutoff: September 26, 2026. The full randomized dose-comparison article was reviewed on September 27, 2026. Full-text journal PDFs are not distributed in this repository.

| Input | Source and location | Use |
|---|---|---|
| [poppk_parameters_final.csv](poppk_parameters_final.csv) | [Nagase et al. AAPS J. 2025;27:26](https://doi.org/10.1208/s12248-024-01013-6), Table II, PDF p7 | Published parameter estimates for the deterministic PK checks. |
| [covariate_targets.json](covariate_targets.json) | Same article, Figure 3, PDF p9 | Twenty published point ratios. The uncertainty bands are not reproduced. |
| [fda_figure26_quartiles.csv](fda_figure26_quartiles.csv) | [FDA NDA 214665 review](https://www.accessdata.fda.gov/drugsatfda_docs/nda/2021/214665Orig1s000MultidisciplineR.pdf), Figure 26 left panel, PDF p248 | Group boundaries, displayed coordinates, response counts, and denominators. The regression uses the displayed coordinates or specified alternatives. |
| Response probabilities of 35% and 25% | Hypothetical design assumptions | Not estimated from the trial or documented sponsor planning assumptions. |

The FDA figure's four groups total 228 patients, while nearby methods describe 248. The figure does not state whether the exposure coordinates are means, medians, or another summary. Both limitations affect interpretation of the grouped analysis.

## Randomized dose-comparison article

[Hochmair et al. EJC. 2024;208:114204](https://doi.org/10.1016/j.ejca.2024.114204), Figure 3, page 7, reports rounded 960/240 ratios: day-1 Cmax 1.4 and AUC0–24h 1.5; day-8 Cmax and AUC0–24h both 1.3. The concentration profiles are labeled mean (SD). The caption gives a PK cutoff of September 9, 2022 and 208 patients overall.

The article does not report exact AUC and Cmax summaries or specify whether the ratios compare arithmetic or geometric means. It also does not give PK denominators for each arm and day, or uncertainty intervals for the ratios. The SD bars describe concentration variability; they are not confidence intervals for a ratio. The supplement and detailed protocol have not been reviewed. These gaps limit the interpretation of the low-dose projections excluded from this project's conclusions.

The efficacy and safety cutoffs were June 23, 2023 and January 18, 2023, respectively. Table 2 reports a stratified ORR difference of 6.9 percentage points (90% CI −3.4 to 17.1), compared with an unstratified difference of 7.9 (−2.3 to 18.2). The methods describe stratification by prior treatment lines, CNS history, race, and ECOG status, and permit dose reduction only in the 960-mg arm. These findings provide clinical context; they are not inputs to the active calculations.

## Regulatory context

[FDA's August 2024 dosage guidance](https://www.fda.gov/media/164555/download), Sections III.A–C, especially III.B, printed p7, informs the interpretation of dose selection and trial design. It is a later, nonbinding framework, not a retrospective requirement or a numerical input.

[Singh, Vellanki, and Pazdur. JCO. 2025;43:248–250](https://doi.org/10.1200/JCO.24.00310), published online October 7, 2024, describes the regulatory dose decision summarized in the README. This source supplies no numerical input, and colorectal combination-therapy evidence is not pooled with the NSCLC monotherapy data.
