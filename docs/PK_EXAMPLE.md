# Worked PK example: albumin and exposure

This example compares albumin concentrations of 30 and 40 g/L in two otherwise identical typical subjects receiving sotorasib 960 mg once daily. Random effects are zero. The published inputs come from [Nagase et al., Table II and Figure 3](https://doi.org/10.1208/s12248-024-01013-6).

The [parameter table](../data/poppk_parameters_final.csv) gives an albumin coefficient of 0.0296. In [the implementation](../src/model.py), the clearance multiplier is:

```text
1 + 0.0296 × (albumin − 40)
At 40 g/L: 1.000
At 30 g/L: 0.704
```

Thus, the model assigns 29.6% lower apparent clearance at 30 g/L. This multiplier applies to both baseline and steady-state clearance; the reference estimates are 21.6 and 41.3 L/hour. Albumin does not change relative bioavailability, absorption, or volume in this comparison. The g/L interpretation is an assumption needed to reconcile inconsistent source descriptions; [parameter notes](../data/PARAMETER_RESOLUTION.md#albumin) explain why approximate agreement supports, but does not prove, it.

Clearance and volume are **apparent** parameters because oral data cannot separately identify absolute bioavailability. Relative bioavailability is scaled to the 960-mg steady-state reference, fixed at 1; this does not mean that 100% of an oral dose reaches the circulation.

With constant parameters at periodic steady state, AUC over one dosing interval is proportional to dose × relative bioavailability / apparent clearance. Because dose and relative bioavailability are the same in both subjects, a useful hand check is:

```text
AUC ratio ≈ 1 / 0.704 = 1.42045
```

The [comparison code](../src/validate.py) instead simulates time-varying clearance and relative bioavailability and integrates the full day-30 interval, from 696 to 720 hours after the first dose, including both endpoints. The [saved results](../results/pk_covariate_checks.csv) are:

| Exposure ratio, 30 versus 40 g/L | Reproduced | Published |
|---|---:|---:|
| AUC over 24 hours | 1.42057 | 1.414 |
| Maximum concentration, Cmax | 1.07412 | 1.075 |

The reciprocal is an approximation: the simulation uses a finite treatment duration, parameters still approaching their limits, and numerical integration. Its AUC ratio differs from the published point value by 0.46%. Rounded published parameters and reconstruction choices can contribute to disagreement; without the original model code, its exact source cannot be established. The [extracted comparison inputs](../data/covariate_targets.json) retain the published values.

Cmax measures the peak rather than total exposure. Absorption, distribution, elimination, and accumulation shape that peak, so it does not follow the reciprocal-clearance rule for AUC.

This checks one published model association. It does not show that changing albumin causes the exposure difference, establish an exposure target, or justify changing a patient's dose. Patient variability, parameter uncertainty, and clinical benefit and harm are not evaluated by this deterministic comparison.
