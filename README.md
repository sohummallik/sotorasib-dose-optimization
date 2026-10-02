# Sotorasib dose evidence
### Reproducing published pharmacokinetic results and examining the 240 mg versus 960 mg comparison

## Project overview

I used the sotorasib 240 mg versus 960 mg comparison to examine what published exposure data can tell us about dose selection. The sources were public FDA reviews, clinical literature, and a published population pharmacokinetic (popPK) model.

I checked selected published model results, examined how the grouped exposure-response estimates change under different assumptions, and calculated sample size for a hypothetical response-rate comparison. The PK and grouped exposure-response analyses use published summaries; the sample-size example uses stated hypothetical response rates. Individual patient records were not available.

FDA considered the dose-related postmarketing requirement fulfilled in December 2023 and retained 960 mg using the randomized comparison and other supporting evidence. This project examines selected calculations within that evidence; it does not determine a preferred starting dose. [Regulatory perspective](https://doi.org/10.1200/JCO.24.00310).

## Analyses and main findings

| Analysis | Finding | Interpretation |
|---|---|---|
| **Published PK covariate comparisons** | All 10 comparisons, comprising 20 exposure ratios, were reproduced within **1.901%** of the published point values. | These results check selected parts of the implementation. They do not validate the model for clinical decisions. |
| **Grouped exposure-response analysis** | The odds ratio per doubling of the assigned exposure was **0.553** using the displayed group coordinates, **0.677** using geometric midpoints, and **0.653** using arithmetic midpoints. | The estimate changes with the choice of representative exposure. Four aggregate groups cannot account for patient-level differences or establish causation. |
| **Hypothetical response-rate comparison** | **329 participants per arm, 658 total**, for response rates of 35% versus 25%, equal allocation, two-sided alpha 0.05, and 80% power under a normal approximation. | This is one hypothetical superiority design. It does not reconstruct the sponsor's plan or set a minimum sample size for dose optimization. |

The conditional 95% confidence intervals were **0.364–0.842** for the displayed coordinates, **0.505–0.907** for geometric midpoints, and **0.480–0.887** for arithmetic midpoints. They do not include uncertainty about the representative exposures or confounding.

## Figures

![Published and reproduced PK covariate comparisons](figures/pk_covariate_checks.png)

**Figure 1.** Published and reproduced covariate effects on exposure at 960 mg once daily. Only the published point values are compared; their uncertainty bands are not reproduced.

![Grouped exposure-response analysis](figures/grouped_coordinate_sensitivity.png)

**Figure 2.** Left: response counts in the four predicted-AUC groups from FDA Figure 26, with exact binomial 95% intervals calculated from those counts. Right: exposure-response estimates using three choices of representative exposure, with conditional Wald 95% intervals. These are unadjusted associations.

## What I take from the analysis

The model checks helped me understand the equations and test my implementation. The grouped analysis showed how the exposure assigned to each group affects the estimate. These calculations could not recover patient-level information or settle the clinical dose decision.

Dose selection brings together activity, safety, tolerability, and exposure. A descriptive randomized comparison can inform that decision without being designed for a formal superiority test. The sample-size calculation here addresses a separate, explicitly stated testing scenario. [FDA dosage guidance, Section III.B](https://www.fda.gov/media/164555/download).

## Files and reproduction

```text
data/       Published parameters, extracted results, and source notes
src/        Python analyses and figures; R model and statistical checks
results/    Saved analysis tables and summaries
figures/    PK comparisons and grouped exposure-response results
tests/      Numerical and statistical checks
docs/       Methods and contributions
reports/    Numerical checks; local logs are generated when run
```

The Python analysis was tested with Python 3.13.5 and the package versions in `requirements-tested.txt`. From the repository root, install those versions in an isolated environment and run:

```bash
python -m pip install -r requirements-tested.txt
python src/run_all.py
```

The runner executes the tests, model checks, exposure-response analysis, sample-size calculation, and figure generation. The calculations are deterministic. The saved [numerical checks](reports/numerical_checks.json) report grid-refinement and integration results. Each run also creates local command logs and an execution record with package versions and input hashes in `reports/`; these generated records are not tracked in Git.

### R checks

After running the Python pipeline, run the two R checks below. Both were tested with R 4.6.1 and use the same published inputs as Python. The PK check requires `rxode2` and `jsonlite`; tested package versions are in [requirements-R.txt](requirements-R.txt). The statistical check uses only base R.

```bash
Rscript src/check_pk.R
Rscript src/check_statistics.R
```

- [check_pk.R](src/check_pk.R) uses `rxode2` for the ten covariate comparisons at 960 mg once daily. It checks a finer time grid and compares AUC and Cmax ratios with Python.
- [check_statistics.R](src/check_statistics.R) fits the three grouped binomial models, calculates exact response intervals, and uses `power.prop.test` for the hypothetical sample-size example.

The R tables are saved as `results/r_*.csv`.

Agreement between implementations checks the calculations under shared assumptions. It does not provide independent clinical evidence or resolve the source limitations below.

## Limitations

- The model uses published parameters and aggregate results. It was not fitted to individual trial data, and agreement with selected published values is not clinical validation.
- The exposure-response analysis cannot adjust for patient characteristics. FDA Figure 26 totals 228 patients, while nearby methods describe 248; that difference remains unresolved.
- Representative exposures and model-parameter interpretations require assumptions. Details are in the [methods](docs/ANALYSIS_NOTES.md) and [parameter notes](data/PARAMETER_RESOLUTION.md).
- Earlier exposure-overlap projections, toxicity/utility calculations, and an incomplete R parameter-estimation example are excluded from the findings above. The [methods](docs/ANALYSIS_NOTES.md#analyses-not-used-in-the-conclusions) explain their limitations.

## Primary sources

1. [FDA. LUMAKRAS NDA 214665 multidisciplinary review (2021)](https://www.accessdata.fda.gov/drugsatfda_docs/nda/2021/214665Orig1s000MultidisciplineR.pdf).
2. [Nagase et al. AAPS Journal. 2025;27:26](https://doi.org/10.1208/s12248-024-01013-6).
3. [Hochmair et al. European Journal of Cancer. 2024;208:114204](https://doi.org/10.1016/j.ejca.2024.114204).
4. [Singh, Vellanki, and Pazdur. Journal of Clinical Oncology. 2025;43:248–250](https://doi.org/10.1200/JCO.24.00310).
5. [FDA. Optimizing the Dosage of Human Prescription Drugs and Biological Products for the Treatment of Oncologic Diseases (August 2024)](https://www.fda.gov/media/164555/download).

Source locations and access limitations are recorded in the [source notes](data/SOURCE_NOTES.md). Evidence cutoff: September 26, 2026; source-access update: September 27, 2026.

Project led by Sohum Mallik. [Contributions](docs/CONTRIBUTIONS.md).
