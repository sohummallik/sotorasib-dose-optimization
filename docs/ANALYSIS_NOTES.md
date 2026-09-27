# Methods and limits of the retained analysis

## Scope

The active pipeline reproduces selected PK point contrasts, evaluates sensitivity of an aggregate exposure-response fit to coordinate choices, and calculates one hypothetical response-rate design scenario. Source review was targeted, not a systematic review. No preregistration, independent clinical validation or new patient-level model fit is claimed. Original analyses excluded after the audit remain intact in the baseline archive.

## PK implementation and numerical measurement

Structural parameter values come from Nagase Table II. Dose enters the first of three pre-central transit compartments, following Figure 1. Time-dependent apparent clearance and relative bioavailability share an induction coefficient. The active comparisons use deterministic typical-subject scenarios with random effects set to zero, 960 mg once daily and one covariate change at a time. Published variances/covariances remain in the provenance table, but the active pipeline does not sample a population or propagate parameter uncertainty.

The equations are integrated separately between oral dose events using LSODA with relative tolerance 1e-8 and absolute tolerance 1e-6. Doses affect the transit compartment, so central plasma concentration is continuous across the event. Each output interval includes both exact boundaries. Duplicate boundary concentrations must agree to relative tolerance 1e-10 before they are collapsed; inconsistent duplicates are rejected rather than choosing a silent pre/post convention. The helper does not support instantaneous central-compartment dosing.

For an externally supplied grid that does not contain a requested boundary, concentration is linearly interpolated between enclosing observations. The helper refuses to extrapolate. AUC is calculated by the trapezoid rule on the closed interval. Cmax/Tmax are sample-grid quantities. C_tau is the concentration at the interval endpoint, not necessarily the minimum within the interval. The historical floor-based interval label dropped each endpoint and has been removed.

The active PK grid is 0.03125 h. Selection is supported by refinement against 0.015625 h and a stricter solve, not by tuning to published ratios. The largest change in any retained covariate ratio when halving the grid was 0.0161%, small relative to the pre-existing 5% source-point comparison threshold. That threshold is retained unchanged; it is a limited implementation acceptance criterion, not a clinical validation threshold. Published uncertainty bands, full fixed-effect covariance, residual-error behavior and the original control stream have not been reproduced.

A second check freezes clearance and F at constant values and simulates to periodic steady state. In this linear system AUC over a complete dose interval equals dose×F/CL, independently of the numerical implementation. Tests compare against that identity. Additional tests cover exact and unaligned boundaries, no extrapolation, oral-event continuity, interval additivity and the previously failing 0.05-h output grid.

## Grouped exposure-coordinate sensitivity

FDA Figure 26, left panel, provides four modeled-AUC groups. The verified boundary pairs are 5000–16000, 16000–22000, 22000–35000 and 35000–85000 h·ng/mL. Displayed central coordinates are separately retained as 13000, 19000, 27000 and 46000. Their summary statistic is not explicitly named in the inspected figure. Counts are 25/57, 26/57, 21/57 and 12/57. The four denominators sum to 228, while nearby methods mention 248; this discrepancy is not filled by assumption.

For each coordinate scheme, logit(p_j)=a+b log(x_j) is fitted to the four binomial counts. The descriptive OR per coordinate doubling is exp(b log 2). Its Wald interval uses b ± 1.96 SE(b). Three schemes are reported: displayed labels, geometric range midpoints and arithmetic range midpoints. The alternatives do not estimate the true mean or median exposures within a bin. Expanding the counts at fixed coordinates into repeated binary observations is used only to check the grouped likelihood; it does not manufacture observed patient-level exposures.

The negative fitted association is not a treatment effect. Disease burden and albumin can affect exposure and prognosis, and these aggregate fits do not adjust for them. The conditional intervals omit within-bin exposure uncertainty and the decision to choose a particular coordinate. The analysis does not put efficacy on the same evidentiary footing as the FDA's patient-level continuous safety analysis.

## Hypothetical ORR design

The retained scenario assumes response probabilities 0.35 and 0.25, independent equal-sized arms, two-sided alpha 0.05 and desired power 0.80. The normal approximation uses a pooled-null variance at the rejection boundary and an unpooled variance under the specified alternative. Both rejection tails are included. A scalar root is solved for n and rounded up; 328.47083 becomes 329 per arm, 658 total. Power at 329 is 0.8006336; at 328 it is 0.7994348.

There is no continuity correction, attrition allowance or adjustment for unevaluable participants. Assessment time, clinically meaningful difference, efficacy objective, tolerability measures and operational feasibility need clinical justification before using such a design. The calculation does not reconstruct the sponsor's assumptions, implement noninferiority or identify a uniquely appropriate next trial. Observed-effect post hoc power is excluded because it does not add independent evidence about the observed trial result.

## Why previous simulation claims are excluded

The original low-dose calibration solved for two F parameters from two selected ratios. Exact trial timepoints, averaging conventions, PK denominators and uncertainty remain unverified. Agreement with those two targets is not goodness-of-fit validation. Most other model parameters were fixed. The sensitivity interval perturbed targets independently at an assumed log SD of 0.15; it was not an empirical clinical confidence interval and did not propagate all structural uncertainty.

With identical virtual subjects and clearance at both doses, steady-state AUC ratio is (240×F240)/(960×F960) for every subject. The original GMR therefore restated the calibrated F assumption. Distributional overlap characterized consequences of additional variability assumptions, not independent evidence that doses are equivalent. The original sensitivity also allowed an increase rather than decrease in low-dose F in 22.75% of draws, and its point estimate failed its own claimed adjacent-dose interpolation check. Those details are preserved in the audit record rather than silently repaired by selecting convenient constraints.

The historical toxicity function anchored a probability at a reference AUC and matched only the lower-dose population mean. Matching both means with an intercept and slope would fix the description numerically but would not justify causality, account for efficacy/safety cutoff differences or establish a meaningful utility weight. The active analysis therefore excludes the utility surface instead of retuning it toward a desired conclusion. The efficacy-only OR threshold was also strongly dependent on the reference anchor. Neither appears in current scientific figures.

## Unresolved source ambiguities and human responsibility

The active implementation follows the published three-transfer figure while documenting a conflicting transit-time formula. Albumin units and random-effect notation are treated as source ambiguities, not author-confirmed errors. The full trial report/supplement remains important for any future calibration extension; it is not required for the retained source-verified grouped data and implementation checks.

Original work involved reported Claude code assistance. This correction, test and documentation cycle involved Codex assistance. No independent human pharmacometric or clinical review is implied. The author must review and understand the code, claims, contribution record and limitations before dissemination.
