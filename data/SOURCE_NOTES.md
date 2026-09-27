# Source and input lineage

Evidence cutoff: September 26, 2026. No copyrighted full-text PDFs are included.

| Input | Source and exact location | Verification and use |
|---|---|---|
| `poppk_parameters_final.csv` | Nagase et al. AAPS J. 2025;27:26, Table II, PDF p7; [DOI](https://doi.org/10.1208/s12248-024-01013-6) | Numeric final-estimate values preserved from the baseline; table and model figure inspected in the source audit. Annotation overclaims corrected. Point estimates drive deterministic checks. |
| `covariate_targets.json` | Same article, Figure 3, PDF p9 | All twenty point ratios visually checked. This is not a transcription of the uncertainty bands. |
| `fda_figure26_quartiles.csv` | [FDA NDA 214665 review](https://www.accessdata.fda.gov/drugsatfda_docs/nda/2021/214665Orig1s000MultidisciplineR.pdf), Figure 26 left panel, PDF p248 | Bounds, labels, responder counts and denominators visually checked separately. The active regression uses labels or explicitly chosen alternatives. |
| Hypothetical 35% and 25% response probabilities | Analyst-specified illustration | Not estimated from the trial and not documented sponsor planning assumptions. |

The FDA review's surrounding methods mention 248 patients, whereas the displayed bins total 228. The source does not explicitly identify what summary statistic its central labels represent. Those limitations remain visible in outputs.

The inspected EJC abstract reports modest exposure separation; it does not verify the exact historical day-specific calibration inputs, averaging method, denominators or uncertainty. Those inputs and their dependent projections are archived and excluded from active claims.

`reports/execution_report.json` records input and code hashes. `archive/baseline-043a805/MANIFEST_SHA256.json` preserves the reviewed baseline. Current source corrections are factual transcriptions and operational interpretation statements, not a claim to possess the sponsor control stream.
