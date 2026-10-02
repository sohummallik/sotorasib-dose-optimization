# Methods

The PK checks and grouped exposure-response analysis use published summaries. The sample-size example uses hypothetical response rates. No model was fitted to individual patient records, and these calculations do not establish which sotorasib dose is preferable.

## PK model checks

The model uses the parameter estimates in Nagase Table II and the three transit compartments shown in Figure 1. Apparent clearance and relative bioavailability change over time using the same induction coefficient. Each comparison uses a typical subject receiving 960 mg once daily, with one covariate changed at a time. Random effects are set to zero, so the comparisons do not include between-patient variability. Uncertainty in the parameter estimates is also not simulated.

LSODA solves the model equations between oral doses, with relative tolerance 1e-8 and absolute tolerance 1e-6. AUC is calculated by trapezoidal integration over a complete dosing interval, including both endpoints. When an interval boundary falls between sampled times, its concentration is linearly interpolated. The calculation does not extrapolate beyond the simulated time range. Cmax and Tmax depend on the sampling times. C_tau is the concentration at the end of the interval; it is not necessarily the lowest concentration during that interval.

Simulation outputs are sampled every 0.03125 hours. Halving that spacing to 0.015625 hours changed the exposure ratios by at most 0.0161%. All 20 ratios were within 1.901% of the published point values. The 5% comparison threshold checks whether the code reproduces those values; it does not establish clinical validity. Tests also check interval boundaries, continuity of the central concentration across oral doses, and the steady-state relationship AUC = dose × F / CL when clearance and F are held constant. The results are in the [numerical checks](../reports/numerical_checks.json).

The checks do not reproduce the published uncertainty bands or residual variability. The original model code was not available. The albumin units and transit-time descriptions also require interpretation; those choices are explained in the [parameter notes](../data/PARAMETER_RESOLUTION.md).

## Grouped exposure-response analysis

FDA Figure 26 shows four predicted-AUC groups with response counts of 25/57, 26/57, 21/57, and 12/57. Together, the four groups contain 228 patients, while nearby methods describe 248. That difference remains unresolved. The figure also does not identify whether its displayed exposure coordinates are means, medians, or another summary.

A grouped binomial logistic regression relates the log odds of response to log exposure. The fit is repeated using the displayed coordinates, geometric range midpoints, and arithmetic range midpoints. The odds ratio per exposure doubling is exp(slope × log 2). Its Wald interval treats the assigned exposures as fixed; it does not account for uncertainty about those exposures, within-group exposure variation, or confounding.

The negative association does not show that higher exposure reduces efficacy. Disease burden and albumin may affect both exposure and prognosis, and the grouped data cannot adjust for those differences. A separate check expands the counts into binary response outcomes at the assigned group exposures to confirm agreement with the grouped fit. This does not recover actual patient-level exposures.

## Hypothetical response-rate comparison

The sample-size calculation assumes response rates of 35% and 25%, equal allocation, two-sided alpha 0.05, and 80% power. Under a normal approximation, pooled variance defines the rejection boundaries under the null hypothesis. Separate variances are used to calculate power under the alternative, including both rejection tails. Rounding up gives 329 participants per arm, or 658 total.

The calculation includes no continuity correction or allowance for attrition. The response rates are hypothetical and do not come from the sponsor's documented planning assumptions. This example addresses a formal superiority test. It does not set a minimum size for a dose-optimization trial or assess noninferiority. A descriptive randomized comparison can still inform dose selection, as discussed in Section III.B of [FDA's August 2024 guidance](https://www.fda.gov/media/164555/download). That guidance is nonbinding and was issued after the reported trial data cutoffs.

## R reproduction checks

[check_pk.R](../src/check_pk.R) implements the same deterministic 960 mg model in `rxode2`, using the source parameter table and covariate targets. It checks the complete day-30 interval at 0.03125-hour and 0.015625-hour spacing and compares the AUC and Cmax ratios with the saved Python results. These checks retain the parameter interpretations described above.

[check_statistics.R](../src/check_statistics.R) uses a binomial `glm` for the four response groups under each choice of exposure coordinates. It compares the estimates and conditional intervals with Python and uses `binom.test` for exact group intervals. The sample-size check uses `power.prop.test(strict = TRUE)` to include both rejection tails and gives the same result: 329 participants per arm. These checks use the same inputs and assumptions as Python; agreement does not provide new clinical evidence.

## Analyses not used in the conclusions

The earlier low-dose calibration fitted two bioavailability parameters to two rounded AUC ratios while holding the other parameters fixed. Matching two targets with two fitted parameters leaves no separate target to check the fit. Population-average results are also not necessarily comparable with predictions for a typical subject. The uncertainty assigned to the targets was assumed, rather than taken from a clinical confidence interval.

With the same simulated subjects and clearance at both doses, the steady-state AUC ratio reduces to (240 × F240) / (960 × F960). The resulting exposure ratio largely restates the bioavailability assumptions. Exposure overlap alone cannot establish clinical equivalence.

The toxicity and utility projections depended on unsupported causal assumptions and chosen weights. The toxicity calibration also did not reproduce both population means. Efficacy thresholds depended on the reference exposure, and the efficacy and safety results used different data cutoffs. These limitations are why the projections are excluded from the current findings. Post hoc power calculated from the observed effect and the incomplete R estimation example are also excluded.

[Source notes](../data/SOURCE_NOTES.md) identify the inputs and remaining source limitations. [Contributions](CONTRIBUTIONS.md) describes project roles and assistance.
