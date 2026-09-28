# Source and input lineage

Evidence-publication cutoff: September 26, 2026. Latest documented source inspection: September 27, 2026. No copyrighted full-text PDFs are included.

| Input | Source and exact location | Verification and use |
|---|---|---|
| `poppk_parameters_final.csv` | Nagase et al. AAPS J. 2025;27:26, Table II, PDF p7; [DOI](https://doi.org/10.1208/s12248-024-01013-6) | Numeric final-estimate values preserved from the baseline; table and model figure inspected. Point estimates drive deterministic checks. Operational choices and reporting ambiguities are documented in [parameter notes](PARAMETER_RESOLUTION.md). |
| `covariate_targets.json` | Same article, Figure 3, PDF p9 | All twenty point ratios visually checked. This is not a transcription of the uncertainty bands. |
| `fda_figure26_quartiles.csv` | [FDA NDA 214665 review](https://www.accessdata.fda.gov/drugsatfda_docs/nda/2021/214665Orig1s000MultidisciplineR.pdf), Figure 26 left panel, PDF p248 | Bounds, labels, responder counts and denominators visually checked separately. The active regression uses labels or explicitly chosen alternatives. |
| Hypothetical 35% and 25% response probabilities | Analyst-specified illustration | Not estimated from the trial and not documented sponsor planning assumptions. |

The FDA review's surrounding methods mention 248 patients, whereas the displayed bins total 228. The source does not explicitly identify what summary statistic its central labels represent. Those limitations remain visible in outputs.

## Trial dose-comparison evidence

A lawfully supplied copy of [Hochmair et al., EJC 2024;208:114204](https://doi.org/10.1016/j.ejca.2024.114204) was inspected on September 27, 2026. Figure 3, page 7, reports rounded 960/240 ratios: day-1 Cmax 1.4 and AUC0–24h 1.5; day-8 Cmax and AUC0–24h both 1.3. The plotted concentration profiles are labeled mean (SD). The caption specifies the PK cutoff of September 9, 2022 and 208 patients overall.

It does not supply exact AUC/Cmax parameter values, explicitly identify the arithmetic/geometric convention of the ratios, give per-arm/per-day PK denominators or quantify ratio uncertainty. Concentration SD bars are not standard errors or confidence intervals for the ratio. The supplement and detailed protocol remain uninspected.

The full text also verifies the June 23, 2023 efficacy and January 18, 2023 safety cutoffs; Table 2 distinguishes the stratified ORR difference 6.9 percentage points (90% CI −3.4 to 17.1) from the unstratified 7.9 (−2.3 to 18.2). The methods describe stratification by prior lines, CNS history, race and ECOG and permit dose reduction only in the 960-mg arm. These clinical clarifications do not change active computational inputs.

The calibration and dependent projections remain archived because the full article does not resolve the exactly identified fit, population-versus-typical estimand, assumed target uncertainty or algebraically imposed exposure ratio. The [correction record](../docs/CORRECTION_LOG.md) preserves the source-access history.

The [baseline manifest](../archive/baseline-043a805/MANIFEST_SHA256.json) preserves the reviewed original files. Current source notes distinguish factual transcriptions from operational interpretations; the sponsor control stream was not obtained.

## Interpretive sources

[FDA final August 2024 dosage guidance](https://www.fda.gov/media/164555/download): Sections III.A–C, especially III.B, printed p7. Used for interpretation, not numerical inputs; no retrospective binding requirement is asserted.

[Singh, Vellanki and Pazdur, JCO 2025;43:248–250](https://doi.org/10.1200/JCO.24.00310), online October 7, 2024: publisher-indexed text inspected, not the full PDF. Used to establish the completed regulatory dose-decision context described in README. It supplies no active numerical input. The colorectal evidence is not pooled with the NSCLC aggregate data.

The reported computation is pinned to [46c83b90f396e942452241fdbe09ae7f1546a7c7](https://github.com/sohummallik/sotorasib-dose-optimization/tree/46c83b90f396e942452241fdbe09ae7f1546a7c7) and the [manifest at that snapshot](https://github.com/sohummallik/sotorasib-dose-optimization/blob/46c83b90f396e942452241fdbe09ae7f1546a7c7/reports/execution_report.json). The [execution report](../docs/EXECUTION_REPORT.md) explains the run and subsequent documentation checks.
