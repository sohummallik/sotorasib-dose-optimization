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

**240-mg calibration and projected overlap.** Two F parameters were calibrated to two ratios with the remainder fixed. Exact final-trial target timepoints, averaging convention, PK denominators and uncertainty remain unverified. The shared-subject steady-state ratio is algebraically determined by the assumed F ratio; simulation adds no independent corroboration of that ratio. The historical interval used assumed independent target perturbations and omitted most parameter/model uncertainty. The claimed induction-pattern consistency was false (19.18% versus the cited adjacent-dose range 7.5%–16.6%). These analyses and figures remain visible only as excluded historical work.

**Post hoc power and R.** Observed-effect post hoc power is removed from active evidence. The explicit prospective-style 35% versus 25% illustration remains with all assumptions. R outputs are archived: their incomplete covariance/SE step, diagonal random effects and reduced-scale synthetic run are not represented as successful clinical estimation or validation. No R rerun was required by the retained analysis.

## Dependency reconciliation

| Current input/method | Current results | Figure/document dependencies |
|---|---|---|
| `poppk_parameters_final.csv`, `covariate_targets.json`, corrected model/interval logic | `pk_covariate_checks.csv`, `pk_check_summary.json`, `reports/numerical_checks.json` | `pk_covariate_checks.png/.pdf`; README; methods and numerical report |
| `fda_figure26_quartiles.csv`, grouped likelihood and coordinate schemes | `grouped_coordinate_sensitivity.csv`, `grouped_exposure_response.json` | `grouped_coordinate_sensitivity.png/.pdf`; README; methods |
| Explicit hypothetical response probabilities and alpha/power | `design_scenario.json` | README; methods; no clinical recommendation figure |
| Supported runner and focused tests | Command logs and `execution_report.json` | Execution report and reproducibility instructions |

All prior results and figures have been removed from the active output directories and retained in the baseline archive. Current figure generation uses only current retained outputs. No preprint, formal release or application has been submitted by these corrections. No institutional, funding or conflict declaration is inferred from missing information.
