# Execution and numerical verification, v0.2

## Baseline reproduction

The reviewed baseline is `043a80508e35306f2ddaac6edc9a621c12a25b21`. Before editing, the original validation, calibration, population and decision scripts all exited successfully in the correction checkout. Their logs and exact timings are in `reports/baseline_*.log` and `reports/baseline_defects.json`. The original endpoint omission, 0.05-h sampling failure, incorrect quartile-bound fields and toxicity calibration-description mismatch were reproduced. This confirms computational behavior, not scientific validity.

## Supported corrected commands

The final command is `python src/run_all.py` from the repository root. It executes:

1. `python -m unittest discover -s tests -v`
2. `python src/validate.py`
3. `python src/numerical_checks.py`
4. `python src/grouped_exposure_response.py`
5. `python src/design_scenario.py`
6. `python src/figures.py`

The exact executable, dependency versions, command statuses, elapsed times, code/input hashes and output paths are recorded in `reports/execution_report.json`. This was a documented existing Anaconda environment, not a claim of a newly provisioned isolated environment. Tested top-level versions are in `requirements-tested.txt`. The active pipeline has no random draws or seeds; archived seed-dependent analyses are excluded.

## Focused tests and tolerances

Fourteen tests cover closed endpoints, unaligned-boundary interpolation, refusal to extrapolate, rejection of ambiguous discontinuous duplicate samples, finite inputs, oral-dose central continuity, integral partition additivity, a 30-dose 0.05-h grid, terminal-versus-minimum semantics, analytical mass balance, grid/solver sensitivity, exact FDA field transcription, grouped likelihood equivalence and design arithmetic/null power.

Numerical error budgets were chosen to be small compared with the pre-existing 5% point-contrast threshold: relative AUC budget 0.0005 (0.05%), relative peak-grid budget 0.001 (0.1%), solver-relative budget 0.00001 (0.001%), and covariate-ratio grid budget 0.1%. They were not weakened to obtain agreement. The frozen-parameter model is checked against the independent identity AUC=dose×F/CL. Synthetic linear concentration curves check boundary interpolation analytically. These checks address the actual defects, not merely similarity with archived output.

With the selected 0.03125-h grid, the frozen-parameter mass-balance relative error is 1.254e-7. Halving the grid to 0.015625 h changes no retained covariate ratio by more than 0.0161%. Against that finer grid and a stricter solver, typical day-8 differences are 0.00001009% for AUC, −0.007616% for Cmax and −0.000002356% for terminal concentration. Full numbers are in `reports/numerical_checks.json`.

All ten source-point contrasts pass the unchanged 5% threshold; the largest discrepancy is 1.9007% for the high-fat Cmax comparison. Selected point agreement does not validate absolute patient predictions, residual variation, unobserved doses or the original control stream. Source ambiguities remain in `data/PARAMETER_RESOLUTION.md`.

## Statistical and figure checks

FDA bounds, coordinates, counts and denominators are tested as separate source-linked fields. A second likelihood implementation, repeating grouped counts at their fixed coordinates, agrees with the grouped-binomial coefficient; this is a mathematical check, not recovered patient data. The design calculation crosses its target power between 328 and 329 per arm and agrees with an independently expressed standard planning formula. No observed-effect post hoc power is active.

Two current scientific figures were regenerated as PNG and PDF and visually inspected for readable labels, complete intervals, units and clipping. The grouped figure's left intervals are Clopper–Pearson intervals calculated here from source counts; right intervals are conditional Wald model intervals. Captions retain these distinctions. No archived figure is used as a current result.

## Checks not performed

No R estimation/covariance rerun, patient-level trial fit, causal exposure-response validation, original control-stream replication, comprehensive fixed-effect uncertainty propagation or independent human specialist review was performed. The full trial article was inspected after the original correction and verifies the rounded day-1/day-8 ratios and overall PK analysis population. Remaining calibration-definition and uncertainty questions, together with the previously documented methodological limitations, prevent reinstatement of excluded projection claims. These questions do not affect the retained source-verified numerical inputs. Successful execution establishes only the tests and calculations described here.


## Documentation-only source-access update

On September 27, 2026, the full Hochmair EJC article became available and was inspected, including visual review of Figure 3. README, source/parameter notes, analysis notes and the correction/execution reports were updated. Active code, numerical inputs, results and figures were unchanged, so no numerical rerun was required. The existing machine-readable execution report continues to describe the preceding computation. Its original hashes for the two Markdown documents under `data/` preserve their pre-addendum versions; the transparent `post_execution_documentation_update` entry records current document hashes and the unchanged-code/numerical-input check. Derived results and figures were also checked byte-for-byte against the preceding correction commit. Documentation changes do not alter any reported test or numerical result.
