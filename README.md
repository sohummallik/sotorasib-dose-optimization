# Sotorasib dose selection: evidence and reproducible analyses

What can published pharmacokinetic and clinical evidence establish about sotorasib dose selection in KRAS G12C-mutated non-small-cell lung cancer?

This project examines that question through implementation checks of a published population pharmacokinetic (PK) model, sensitivity analysis of grouped exposure-response data, and a hypothetical response-rate testing scenario. The [randomized 960-mg versus 240-mg trial](https://doi.org/10.1016/j.ejca.2024.114204) provides clinical context. The analyses use published aggregate evidence; they do not fit individual trial data or develop a new population PK model.

## Supported findings

| Analysis and purpose | Retained result | Interpretation |
|---|---|---|
| **Published-model implementation checks:** compare AUC and peak-concentration ratios for ten covariate contrasts with Nagase Figure 3. | All 20 point ratios fall within the retained 5% numerical comparison threshold; the largest difference is **1.901%**. | Agreement supports this limited implementation check. It does not establish clinical prediction accuracy or reproduce the original model-development process. |
| **Grouped exposure-coordinate sensitivity:** fit the same four FDA response-count groups using three representations of exposure. | Odds ratios per coordinate doubling are **0.553**, **0.677**, and **0.653** for displayed labels, geometric midpoints, and arithmetic midpoints, respectively. | The fitted magnitude depends on the assigned bin coordinates. These unadjusted associations do not establish a causal exposure effect. |
| **Hypothetical formal response-rate test:** calculate sample size under explicit assumptions. | **329 participants per arm, 658 total**, for 35% versus 25% response, equal allocation, two-sided alpha 0.05, and 80% power. | This is an illustration of formal hypothesis testing, not the sponsor's design or a minimum sample size for dosage selection. |

The conditional 95% intervals for the grouped odds ratios are 0.364–0.842, 0.505–0.907, and 0.480–0.887. They omit uncertainty in coordinate choice, within-bin exposure, and confounding. Figure 26's four groups total 228 patients, while the nearby FDA methods refer to 248; this difference remains unresolved.

![Selected PK covariate point checks](figures/pk_covariate_checks.png)

**Figure 1.** Published and implemented point contrasts at 960 mg once daily, day 30. AUC is the area under the concentration-time curve. Comparisons use a 0.03125-hour sampling grid; solver and grid checks are recorded in the [numerical report](reports/numerical_checks.json). Source covariate labels are retained. Published uncertainty bands are not reproduced.

![Grouped exposure-coordinate sensitivity](figures/grouped_coordinate_sensitivity.png)

**Figure 2.** Left: response counts from FDA Figure 26, with Clopper–Pearson exact binomial 95% intervals calculated from those counts. Right: conditional Wald 95% intervals from grouped-binomial logistic fits. Midpoints are assumed coordinates, not recovered patient exposures. Source AUC units are h·ng/mL. This is a project-generated display of aggregate data.

## Limits of interpretation

Dosage selection integrates activity, safety, tolerability, and exposure. It differs from establishing superiority or noninferiority, as described in [FDA's nonbinding August 2024 guidance, Section III.B](https://www.fda.gov/media/164555/download). FDA authors subsequently reported that the dose-related postmarketing requirement was fulfilled in December 2023, retaining 960 mg on the broader evidence package ([Singh et al.](https://doi.org/10.1200/JCO.24.00310)). These limited calculations do not reconstruct that assessment or identify an optimal dose.

- PK checks use deterministic typical-subject scenarios with random effects set to zero. Published model-notation ambiguities, parameter uncertainty, residual variability, and the original control stream are not resolved by selected point agreement.
- Grouped response data cannot recover individual exposures or adjust for disease burden and other confounders.
- The trial report verifies rounded 960/240 AUC ratios of 1.5 on day 1 and 1.3 on day 8, but leaves the ratio averaging convention, per-arm/per-day PK denominators, and ratio uncertainty unresolved. Its supplement and detailed protocol remain uninspected.
- Calibrated low-dose exposure-overlap projections, toxicity/utility calculations, efficacy thresholds, observed-effect post hoc power, and incomplete R estimation results are excluded from current findings. Their original files remain identifiable in the [superseded baseline archive](archive/README.md).

The [methods](docs/ANALYSIS_NOTES.md) explain these boundaries and the [parameter notes](data/PARAMETER_RESOLUTION.md) document source ambiguities.

## Reproduce the results

The reported calculations are fixed at [commit `46c83b90f396e942452241fdbe09ae7f1546a7c7`](https://github.com/sohummallik/sotorasib-dose-optimization/tree/46c83b90f396e942452241fdbe09ae7f1546a7c7). Its [execution manifest](https://github.com/sohummallik/sotorasib-dose-optimization/blob/46c83b90f396e942452241fdbe09ae7f1546a7c7/reports/execution_report.json) records Python 3.13.5, dependency versions, six successful commands, and 14 passing tests. Later documentation changes do not change these calculations.

From a fresh clone on macOS or Linux, using Python 3.13:

```bash
git clone https://github.com/sohummallik/sotorasib-dose-optimization.git
cd sotorasib-dose-optimization
git checkout 46c83b90f396e942452241fdbe09ae7f1546a7c7
python3.13 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-tested.txt
python src/run_all.py
```

On Windows, create the environment with `py -3.13 -m venv .venv` and activate with `.venv\Scripts\activate` instead. The isolated-environment instructions are a reproduction route; the recorded run used the existing Anaconda environment identified in the manifest.

The deterministic runner executes tests, PK checks, numerical checks, grouped sensitivity, the hypothetical design calculation, and figure generation. It writes results, figures, logs, and a new execution manifest. It never runs archived analyses. Individual commands, tolerances, and checks not performed are listed in the [execution documentation](docs/EXECUTION_REPORT.md).

## Repository map and sources

| Location | Contents |
|---|---|
| [`src/`](src/) and [`tests/`](tests/) | Retained analyses and focused numerical/statistical tests. |
| [`data/`](data/) | Published point estimates, grouped counts, and [source locations and verification notes](data/SOURCE_NOTES.md). |
| [`results/`](results/) and [`figures/`](figures/) | Retained numerical outputs and figures in PNG/PDF. |
| [`docs/ANALYSIS_NOTES.md`](docs/ANALYSIS_NOTES.md) | Methods, assumptions, and exclusions. |
| [`docs/EXECUTION_REPORT.md`](docs/EXECUTION_REPORT.md) and [`reports/`](reports/) | Execution details, original logs, and verification records. |
| [`docs/CORRECTION_LOG.md`](docs/CORRECTION_LOG.md) and [`archive/`](archive/README.md) | Correction history and the unchanged, superseded baseline. |

Primary sources are [Nagase et al., AAPS Journal 2025](https://doi.org/10.1208/s12248-024-01013-6), [FDA's 2021 multidisciplinary review](https://www.accessdata.fda.gov/drugsatfda_docs/nda/2021/214665Orig1s000MultidisciplineR.pdf), [Hochmair et al., EJC 2024](https://doi.org/10.1016/j.ejca.2024.114204), [FDA's 2024 dosage guidance](https://www.fda.gov/media/164555/download), and [Singh et al., JCO 2025](https://doi.org/10.1200/JCO.24.00310). [Aung et al., JCO Oncology Practice 2026](https://doi.org/10.1200/OP-25-01315) provides subsequent context. Exact source locations and access limitations are in the source notes. Evidence-publication cutoff: September 26, 2026; source-access update: September 27, 2026. Journal PDFs are not distributed here.

## Provenance and assistance

Sohum Mallik is the sole human author and conducted the original project work with Claude assistance in writing Python and R code. Codex later assisted with source checking, analytical and code corrections, testing, figure generation, and substantial drafting and editing of repository documentation. This statement describes repository work; separate manuscripts and writing samples have their own assistance histories. Tool-assisted checks do not constitute independent human scientific review or confirm the author's final approval.
