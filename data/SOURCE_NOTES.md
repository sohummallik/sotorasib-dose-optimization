# Source and input lineage

Evidence-publication cutoff: September 26, 2026. Source-access update: September 27, 2026. No copyrighted full-text PDFs are included.

| Input | Source and exact location | Verification and use |
|---|---|---|
| `poppk_parameters_final.csv` | Nagase et al. AAPS J. 2025;27:26, Table II, PDF p7; [DOI](https://doi.org/10.1208/s12248-024-01013-6) | Numeric final-estimate values preserved from the baseline; table and model figure inspected in the source audit. Annotation overclaims corrected. Point estimates drive deterministic checks. |
| `covariate_targets.json` | Same article, Figure 3, PDF p9 | All twenty point ratios visually checked. This is not a transcription of the uncertainty bands. |
| `fda_figure26_quartiles.csv` | [FDA NDA 214665 review](https://www.accessdata.fda.gov/drugsatfda_docs/nda/2021/214665Orig1s000MultidisciplineR.pdf), Figure 26 left panel, PDF p248 | Bounds, labels, responder counts and denominators visually checked separately. The active regression uses labels or explicitly chosen alternatives. |
| Hypothetical 35% and 25% response probabilities | Analyst-specified illustration | Not estimated from the trial and not documented sponsor planning assumptions. |

The FDA review's surrounding methods mention 248 patients, whereas the displayed bins total 228. The source does not explicitly identify what summary statistic its central labels represent. Those limitations remain visible in outputs.

## Full trial article inspected after the initial correction

A lawfully supplied copy of [Hochmair et al., EJC 2024;208:114204](https://doi.org/10.1016/j.ejca.2024.114204) was inspected on September 27, 2026. Figure 3, page 7, reports rounded 960/240 ratios: day-1 Cmax 1.4 and AUC0–24h 1.5; day-8 Cmax and AUC0–24h both 1.3. The plotted concentration profiles are labeled mean (SD). The caption specifies the PK cutoff of September 9, 2022 and 208 patients overall. This resolves the existence, day assignment and rounded magnitude of the historical calibration ratios.

It does not supply exact AUC/Cmax parameter values, explicitly identify the arithmetic/geometric convention of the ratios, give per-arm/per-day PK denominators or quantify ratio uncertainty. Concentration SD bars are not standard errors or confidence intervals for the ratio. The supplement and detailed protocol remain uninspected. No article PDF is placed in the public repository.

The full text also verifies the June 23, 2023 efficacy and January 18, 2023 safety cutoffs; Table 2 distinguishes the stratified ORR difference 6.9 percentage points (90% CI −3.4 to 17.1) from the unstratified 7.9 (−2.3 to 18.2). The methods describe stratification by prior lines, CNS history, race and ECOG and permit dose reduction only in the 960-mg arm. These clinical clarifications do not change active computational inputs.

The source-access limitation at the initial audit is retained in the correction history. The calibration and dependent projections remain archived because source access does not resolve the exactly identified fit, population-versus-typical estimand, assumed target uncertainty or algebraically imposed exposure ratio.

`reports/execution_report.json` records input and code hashes. `archive/baseline-043a805/MANIFEST_SHA256.json` preserves the reviewed baseline. Current source corrections are factual transcriptions and operational interpretation statements, not a claim to possess the sponsor control stream.
