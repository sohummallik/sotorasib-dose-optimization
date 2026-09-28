# Scientific correction record, v0.2

Baseline: `043a80508e35306f2ddaac6edc9a621c12a25b21`. All 45 tracked baseline files were preserved byte-for-byte under `archive/baseline-043a805`; a SHA-256 manifest protects that record. The original four Python analysis scripts were rerun successfully before correction. `reports/baseline_defects.json` and baseline command logs record reproduction of the defects against that source. This correction branch does not modify default-branch history.

## Confirmed corrections

| Issue | Baseline behavior | Current action and consequence |
|---|---|---|
| FDA Figure 26 bounds | Incorrect lower/upper arrays, with correct plotting coordinates and counts | Correct bounds to 5000–16000, 16000–22000, 22000–35000 and 35000–85000. Separate verified bounds from displayed coordinates. A source-linked regression test protects all fields. The old fit used coordinates, so correcting bounds alone does not change its coefficient. |
| Interval endpoints | Exact endpoint assigned to the next dose by floor(time/24), excluding the final segment | Closed-interval selection includes both endpoints and rejects extrapolation. At unchanged 0.25-h grid, typical day-8 AUC changes from 29,749.027 to 29,801.897 h·ng/mL and terminal concentration from 213.294 to 209.663 ng/mL. |
| Nonbinary sample grid | dt=0.05 h can place a sample beyond the solver interval | Build relative offsets and append the exact endpoint; test 30 doses at 0.05 h and a nondividing 0.7-h grid. |
| Terminal nomenclature | Cmin actually meant last sampled concentration | Use C_tau; document that it is not necessarily the interval minimum. |
| Peak and trapezoid resolution | Active grid 0.25 h | Use 0.03125 h with documented grid refinement. The maximum retained point-contrast discrepancy changes from the baseline's listed 2.4% to 1.9007%; this does not eliminate source/model ambiguities or establish clinical validity. |
| Grouped model description | Implied continuous patient-level footing | Describe a four-bin binomial model at chosen coordinates; add geometric/arithmetic range-midpoint sensitivity and conditional intervals. No patient-level exposures or covariate adjustment exist. |
| Source/parameter assertions | Categorical author-error, nonidentifiability and selection-bias statements | Replace with operational choices, reporting ambiguities and evidence boundaries. Numeric parameter estimates are unchanged. |

The active grouped fits give OR per coordinate doubling of 0.5534, 0.6768 and 0.6525, with conditional intervals 0.3637–0.8421, 0.5048–0.9073 and 0.4799–0.8872. The alternative representations do not establish true within-bin mean exposures or a causal effect.

## Components excluded rather than retuned

**Toxicity/utility.** Reproduction confirmed population predictions of 58.4214% and 49.0000% despite text claiming calibration to 61.5% and 49.0%. The first target had been anchored at a reference exposure, not its population mean. Matching both means with two fitted parameters would correct that numerical description but would not justify the causal exposure-toxicity shape, treatment-emergent-event endpoint or utility weights. Clinical efficacy and safety summaries also have different cutoffs. The component, its outputs and dependent utility panel are therefore archived and excluded. There is no replacement utility recommendation.

**Efficacy threshold.** The old OR 5.1 calculation matched a 7.9-point difference while predicting marginal response probabilities of 36.90% and 29.00%, rather than the observed rates. Its reference anchor and unverified exposure projection determine the result. It is excluded instead of presented as a biological requirement.

**240-mg calibration and projected overlap.** Two F parameters were calibrated to two ratios with the remainder fixed. At the initial correction, the final article had not been inspected and the day-specific ratios and their source definitions were unverified. The source-access update below resolves the rounded ratios, timepoints and overall PK analysis count, but not exact parameter values, per-arm/per-day counts, the explicit ratio averaging convention or empirical target uncertainty. The shared-subject steady-state ratio is algebraically determined by the assumed F ratio; simulation adds no independent corroboration of that ratio. The historical interval used assumed independent target perturbations and omitted most parameter/model uncertainty. The claimed induction-pattern consistency was false (19.18% versus the cited adjacent-dose range 7.5%–16.6%). These analyses and figures remain visible only as excluded historical work.

**Post hoc power and R.** Observed-effect post hoc power is removed from active evidence. The explicit prospective-style 35% versus 25% illustration remains with all assumptions. R outputs are archived: their incomplete covariance/SE step, diagonal random effects and reduced-scale synthetic run are not represented as successful clinical estimation or validation. No R rerun was required by the retained analysis.

## Dependency reconciliation

| Current input/method | Current results | Figure/document dependencies |
|---|---|---|
| `poppk_parameters_final.csv`, `covariate_targets.json`, corrected model/interval logic | `pk_covariate_checks.csv`, `pk_check_summary.json`, `reports/numerical_checks.json` | `pk_covariate_checks.png/.pdf`; README; methods and numerical report |
| `fda_figure26_quartiles.csv`, grouped likelihood and coordinate schemes | `grouped_coordinate_sensitivity.csv`, `grouped_exposure_response.json` | `grouped_coordinate_sensitivity.png/.pdf`; README; methods |
| Explicit hypothetical response probabilities and alpha/power | `design_scenario.json` | README; methods; no clinical recommendation figure |
| Supported runner and focused tests | Command logs and `execution_report.json` | Execution report and reproducibility instructions |

All prior results and figures have been removed from the active output directories and retained in the baseline archive. Current figure generation uses only current retained outputs. No preprint, formal release or application has been submitted by these corrections. No institutional, funding or conflict declaration is inferred from missing information.


## Source-access addendum, September 27, 2026

After the initial correction, the author supplied the full Hochmair EJC article. Figure 3 was visually inspected and supports the rounded day-1 AUC ratio of 1.5 and day-8 ratio of 1.3 used historically; day-1 Cmax is 1.4-fold and day-8 Cmax 1.3-fold. Its profiles are labeled mean (SD), and its caption gives a September 9, 2022 cutoff and 208 patients overall. The full text also confirms stratified/unstratified response differences, randomization strata, distinct efficacy/safety cutoffs and the restriction of dose reductions to the 960-mg group.

This corrects the documentation's previous present-tense source-access gap. It does not convert rounded means into individual or typical-patient observations, derive ratio uncertainty from profile SD bars, validate the exactly identified calibration, or independently confirm the simulated GMR. The active analysis remains unchanged. Exact parameter summaries, the explicit averaging convention of the ratios, per-arm/per-day PK sample counts and ratio uncertainty remain unresolved; the supplement and detailed protocol remain uninspected.

Six documentation files and an explicitly labeled execution-metadata addendum changed in this update. No source code, numerical parameter, active result, figure or archived baseline file changed, and a numerical rerun was not warranted. No copyrighted full-text PDF was copied into the repository.


## v0.3 interpretation and reproducibility update

The added JCO regulatory perspective places this audit after the completed dose decision; it does not provide new computational inputs. The FDA framework now distinguishes dosage selection from a formal efficacy claim. The 658-person illustration is explicitly not a dosage-optimization minimum. Numerical results and exclusions are unchanged.

The full active implementation, tests, numerical input files, dependency declarations, output tables and figure code were inspected again. No material uncovered defect in the retained fixed-input analysis justified another analysis, test or numerical change. This decision does not claim comprehensive software correctness or clinical validation.

Reproducibility links now identify the exact computational snapshot `46c83b90f396e942452241fdbe09ae7f1546a7c7` and its execution manifest; `1e46c8cbe37b6e9bbdf9e483908f8e1c591aaa46` remains the separate source-access documentation update. The original execution record is preserved. A fresh content-hash check, recorded in `reports/documentation_review_v03.json`, verifies unchanged code, tests, numerical inputs, dependencies, results, figures and baseline archive. The pipeline was not rerun for these documentation-only changes.


### Repository presentation

The README now uses a single concise contribution statement. Repeated assistance wording and author-directed review instructions were removed from the analysis notes. Scientific limitations and the assistance record remain explicit. This presentation change does not alter code, data or results.


## Documentation cleanup, September 27, 2026

The README now leads with the research question, retained findings and their limitations, followed by figures, reproduction instructions, navigation and one contribution statement. Active methods, source notes and execution documentation consolidate repeated update narratives; the archive index clearly separates historical files from current results. The contribution statement distinguishes the original project, original coding assistance, later computational/source/documentation assistance and separate writing artifacts.

This update changes documentation only. Code, tests, numerical inputs, dependencies, saved results, figures, execution records and all 45 baseline files remain unchanged. The computational reference remains `46c83b90f396e942452241fdbe09ae7f1546a7c7`. No new numerical run or test execution was performed for this cleanup.
