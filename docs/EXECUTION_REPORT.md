# Execution and numerical verification

## Recorded computational version

The results below belong to corrected snapshot [46c83b90f396e942452241fdbe09ae7f1546a7c7](https://github.com/sohummallik/sotorasib-dose-optimization/tree/46c83b90f396e942452241fdbe09ae7f1546a7c7). Its [execution manifest](https://github.com/sohummallik/sotorasib-dose-optimization/blob/46c83b90f396e942452241fdbe09ae7f1546a7c7/reports/execution_report.json) records the environment, commands, statuses and original hashes. The fourteen passing tests and all numerical values refer to that run. Later documentation revisions do not represent new test or numerical execution.

## Reproduction commands

The pipeline command is `python src/run_all.py` from the repository root. It executes:

1. `python -m unittest discover -s tests -v`
2. `python src/validate.py`
3. `python src/numerical_checks.py`
4. `python src/grouped_exposure_response.py`
5. `python src/design_scenario.py`
6. `python src/figures.py`

The recorded run used an existing Anaconda environment, not a newly provisioned isolated environment. Tested top-level versions are in [`requirements-tested.txt`](../requirements-tested.txt). The active pipeline has no random draws or seeds; archived seed-dependent analyses are excluded.

## Focused tests and tolerances

Fourteen tests cover closed endpoints, unaligned-boundary interpolation, refusal to extrapolate, rejection of ambiguous discontinuous duplicate samples, finite inputs, oral-dose central continuity, integral partition additivity, a 30-dose 0.05-h grid, terminal-versus-minimum semantics, analytical mass balance, grid/solver sensitivity, exact FDA field transcription, grouped likelihood equivalence and design arithmetic/null power.

Numerical error budgets were chosen to be small compared with the pre-existing 5% point-contrast threshold: relative AUC budget 0.0005 (0.05%), relative peak-grid budget 0.001 (0.1%), solver-relative budget 0.00001 (0.001%), and covariate-ratio grid budget 0.1%. The frozen-parameter model is checked against the independent identity AUC=dose×F/CL. Synthetic linear concentration curves check boundary interpolation analytically.

With the selected 0.03125-h grid, the frozen-parameter mass-balance relative error is 1.254e-7. Halving the grid to 0.015625 h changes no retained covariate ratio by more than 0.0161%. Against that finer grid and a stricter solver, typical day-8 differences are 0.00001009% for AUC, −0.007616% for Cmax and −0.000002356% for terminal concentration. Full numbers are in the [numerical-check record at the computational snapshot](https://github.com/sohummallik/sotorasib-dose-optimization/blob/46c83b90f396e942452241fdbe09ae7f1546a7c7/reports/numerical_checks.json).

All ten source-point contrasts pass the unchanged 5% threshold; the largest discrepancy is 1.9007% for the high-fat Cmax comparison. Selected point agreement does not validate absolute patient predictions, residual variation, unobserved doses or the original control stream. Source ambiguities remain in the [parameter notes](../data/PARAMETER_RESOLUTION.md).

## Statistical and figure checks

FDA bounds, coordinates, counts and denominators are tested as separate source-linked fields. A second likelihood implementation, repeating grouped counts at their fixed coordinates, agrees with the grouped-binomial coefficient; this is a mathematical check, not recovered patient data. The design calculation crosses its target power between 328 and 329 per arm and agrees with an independently expressed standard planning formula. No observed-effect post hoc power is active.

Two current scientific figures were regenerated as PNG and PDF and visually inspected for readable labels, complete intervals, units and clipping. The grouped figure's left intervals are Clopper–Pearson intervals calculated here from source counts; right intervals are conditional Wald model intervals. Captions retain these distinctions. No archived figure is used as a current result.

## Checks not performed

No R estimation/covariance rerun, patient-level trial fit, causal exposure-response validation, original control-stream replication, comprehensive fixed-effect uncertainty propagation or independent human specialist review was performed. The full trial article was inspected after the original correction and verifies the rounded day-1/day-8 ratios and overall PK analysis population. Remaining calibration-definition and uncertainty questions, together with the previously documented methodological limitations, prevent reinstatement of excluded projection claims. These questions do not affect the retained source-verified numerical inputs. Successful execution establishes only the tests and calculations described here.

## Baseline reproduction

The reviewed baseline is `043a80508e35306f2ddaac6edc9a621c12a25b21`. Before editing, the original validation, calibration, population and decision scripts all exited successfully in the correction checkout. Their logs and exact timings are in `reports/baseline_*.log` and `reports/baseline_defects.json`. The original endpoint omission, 0.05-h sampling failure, incorrect quartile-bound fields and toxicity calibration-description mismatch were reproduced. This confirms computational behavior, not scientific validity.

## Subsequent documentation checks

The [source-access update](https://github.com/sohummallik/sotorasib-dose-optimization/commit/1e46c8cbe37b6e9bbdf9e483908f8e1c591aaa46) and the [v0.3 verification record](../reports/documentation_review_v03.json) document later source review and unchanged code, tests, numerical inputs, dependencies, results, figures and baseline archive. Neither involved a new numerical run. Detailed chronology is retained in the [correction record](CORRECTION_LOG.md).

The original hashes for Markdown files under `data/` in [`reports/execution_report.json`](../reports/execution_report.json) identify their versions at execution. Its `post_execution_documentation_update` entry identifies their versions at the September 27 source-access update. These are historical hashes, not assertions that later documentation revisions have identical content. The immutable manifest above remains the reference for the reported computation.
