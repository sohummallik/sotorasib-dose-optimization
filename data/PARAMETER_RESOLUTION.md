# Resolution of open parameter questions

Primary source for parameter values: Nagase M, Houk B, Vuu I, Cardona P, Dutta S, Lin CW.
"Population Pharmacokinetics of Sotorasib in Healthy Subjects and Advanced Solid Tumor Patients
Harboring a KRAS-G12C Mutation from Phase 1 and Phase 2 Studies." AAPS J. 2025;27:26.
doi:10.1208/s12248-024-01013-6. PMID 39806205. Table II, Final Model Estimate column.

Cross-reference: FDA NDA 214665 Multi-Discipline Review, Table 35, PDF pp. 237-239.

The two sources agree to rounding on every parameter. Where they differ in the third significant
figure (V2 220 vs 219; Ka 7.87 vs 7.854; CLECOGBL1 0.151 vs 0.149; F1BS_DG2 3.43 vs 3.44),
the AAPS values are used, since FDA's Table 35 reproduces an earlier version of Amgen's
popPK report. Both values are retained in the CSV.

---

## Item 1: OMEGA units. RESOLVED.

The point estimate is the **variance on the log scale**. The parenthetical number is the
**%CV**, not the %RSE, computed for a lognormal distribution as sqrt(exp(omega^2) - 1) * 100.

Reproduced exactly:

| Parameter | omega^2 | sqrt(exp(w2)-1)*100 | Reported |
|---|---|---|---|
| V2/F | 0.285 | 57.4% | 57.4 |
| CL/F | 0.248 | 53.1% | 53.0 |
| Ka | 0.306 | 59.8% | 59.8 |

Note the naive sqrt(omega^2) gives 53.4 / 49.8 / 55.3, which does NOT match. The lognormal
formula does. This confirms both the units and the distributional assumption.

The off-diagonal rows are **covariances on the log scale**; the parenthetical is correlation x 100:

| Pair | Covariance | cov/(sd_a*sd_b)*100 | Reported |
|---|---|---|---|
| V2-CL | 0.178 | 66.95 | 67.0 |
| V2-Ka | -0.113 | -38.26 | -38.4 |
| CL-Ka | 0.021 | 7.62 | 7.64 |

**Consequence: the FDA Table 35 column header "%RSE" is wrong for every IIV row.** Those
numbers are %CV and correlations, not relative standard errors. Precision on the OMEGA
estimates is not reported in either source; use the bootstrap column as the stability check.

IIV is on V2, CL, and Ka only. There is **no IIV on F1**, which matters: the dose-group
bioavailability terms are fixed effects with no between-subject spread of their own.

Caution on notation: AAPS Methods writes the IIV model as log(p_i) = log(theta + eta_i),
which is non-standard and almost certainly a typesetting error for log(theta) + eta_i
(equivalently p_i = theta * exp(eta_i)). Implement the standard form. The %CV
reproduction above confirms the standard lognormal form is what was actually fitted.

---

## Item 2: CLALB centering and units. RESOLVED, and the published unit label is wrong.

Both FDA Table 35 and AAPS Table II label CLALB1 as "(L/hr)/(g/dL)", implying an additive
effect on clearance with albumin in g/dL. That cannot be right. Test it against the
published forest plot target (AAPS Fig. 3): low albumin (median 30 g/L) vs normal
albumin (median 40 g/L) gives an AUCtau,ss ratio of **1.414**.

| Interpretation | CL ratio | Implied AUC ratio | vs target 1.414 |
|---|---|---|---|
| Additive, g/dL | 0.9704 | 1.031 | fails badly |
| Proportional, g/dL | 0.9704 | 1.031 | fails badly |
| Exponential, g/L | 0.7438 | 1.345 | close but off |
| **Proportional-linear, g/L** | **0.7040** | **1.4205** | **matches** |

So the effect is a **proportional change per g/L of albumin**, centered at **40 g/L**
(the typical-subject value stated in AAPS Methods, and the median of the normal-albumin group):

    CL_typical = CLSS * (1 + 0.0296 * (ALB_gL - 40))

Direction: higher albumin, higher clearance, lower exposure. A patient at 30 g/L has 29.6%
lower CL and about 42% higher AUC than a patient at 40 g/L.

---

## Item 3: Dose-group lumping. RESOLVED, and it is worse than assumed.

AAPS Methods states the doses actually administered were 180, 360, 720, 960 mg QD and
480 mg BID in patients, and 360 or 960 mg single dose in healthy subjects.

**240 mg was never an assigned dose level in the popPK dataset.** The 120 mg and 240 mg
data inside dose group 1 can only have come from protocol-directed dose reductions.

Two consequences, and both must be stated before any exposure claim is made:

1. F1SS_DG1 = 4.58 is driven predominantly by 180 mg data. Applying it to a 240 mg
   assigned dose is an extrapolation across a dose the model never saw as a randomized level.

2. Dose-reduction data are informatively censored. Patients reduce dose because they had
   toxicity, and toxicity correlates with exposure. Any bioavailability estimated partly
   from reduced-dose observations inherits that selection. This works in the direction of
   inflating apparent exposure at low doses.

Naive exposure index (Dose x F1SS, relative to 960 mg), propagating the published 95% CI:

| | F1SS_DG1 | 240 mg / 960 mg exposure ratio |
|---|---|---|
| Point estimate | 4.58 | **1.145** |
| Lower 95% CI | 2.70 | 0.675 |
| Upper 95% CI | 6.46 | 1.615 |

The point estimate says 240 mg gives about 15% more exposure than 960 mg. The interval
spans 0.68 to 1.62, which does not exclude 240 mg being meaningfully underexposed.

**This is a placeholder calculation, not the analysis.** It ignores the covariance between
parameters, ignores IIV, and treats the dose-group term as if it were dose-specific. The real
number comes from Monte Carlo simulation with the full variance-covariance structure. But it
sets the expectation, and it is already the single most decision-relevant number in the project.

---

## Item 4: Absorption parameterization. RESOLVED.

Single first-order rate constant Ka = 7.87 /hr, shared between transit compartments AND from
the last transit compartment into the central compartment. Three transit compartments.

Mean transit time = (n+1)/Ka with n = 3: 4/7.87 = **0.508 h**. Paper reports 0.51 h. Confirmed.

"KAIND" in FDA Table 35 is not an induction term on absorption. It is the Ka random effect.
Kind acts on F1 and CL only. FDA's row labels are misleading here; AAPS Fig. 1 makes the
structure unambiguous:

    Dose --F1--> TR1 --Ka--> TR2 --Ka--> TR3 --Ka--> Central (V2/F, CL/V2) <--K24/K42--> Peripheral

Minor unreconciled discrepancy: applying KAHIGHFAT1 proportionally gives
Ka = 7.87 * (1 - 0.614) = 3.04 and MTT = 1.32 h, but the paper reports 1.23 h for the
high-fat condition, which implies a proportional change of -0.587 rather than -0.614.
Small, and it does not affect the fasted-state simulations that this project runs. Flagged
so it is not mistaken for an implementation bug later.

---

## Item 5 (new): residual error structure.

The error model is not a standard combined additive-plus-proportional. It follows Dosne 2016:

    log(Y_obs) = log(Y_pred) + sqrt( theta_x^2 + theta_y^2 / Y_pred^2 ) * eps_1

with eps_1 ~ N(0,1), theta_x = 0.625 (exponential term), theta_y = 3.8 ng/mL (additive term).
This is additive on the log scale with a concentration-dependent standard deviation. Do not
substitute a conventional prop+add error model without noting the deviation.

---

## Validation target suite

AAPS Fig. 3 gives model-predicted exposure ratios for a defined typical subject
(Caucasian male, NSCLC, ECOG 1, baseline tumor >70 mm, albumin >34 g/L with median 40 g/L,
960 mg QD, fasted, no PPI). These are a free acceptance test. The implementation is not
trusted until it reproduces them.

| Comparison | Cmax,ss ratio | AUCtau,ss ratio |
|---|---|---|
| Low albumin (30 g/L) vs normal | 1.075 | 1.414 |
| ECOG 0 vs ECOG 1 | 0.977 | 0.868 |
| ECOG 2 vs ECOG 1 | 1.017 | 1.083 |
| Asian vs Caucasian | 0.966 | 0.817 |
| Black or African American vs Caucasian | 0.988 | 0.934 |
| Female vs Male | 1.242 | 1.190 |
| Tumor <=70 mm vs >70 mm | 0.982 | 0.892 |
| High-fat meal vs fasted | 1.130 | 1.327 |
| PPI vs no PPI | 0.828 | 0.820 |
| Healthy vs patients | 0.898 | 0.534 |

Additional absolute-scale targets from AAPS Model Based Simulations, 960 mg QD, patients:

| | Observed | Model predicted |
|---|---|---|
| AUC Day 1 (h*ng/mL) | 57300 | 53600 |
| AUC Day 8 (h*ng/mL) | 37200 | 35300 |
| Cmax Day 1 (ng/mL) | 7820 | 6600 |
| Cmax Day 8 (ng/mL) | 5970 | 4970 |

Note AUC falls roughly 35% from Day 1 to Day 8 through autoinduction. Any exposure comparison
between doses must be made at the same time point, at steady state, or the induction
difference alone will dominate. Recall that the steady-state F1 decrease is itself
dose-dependent: 34.2% at 960 mg vs 7.5% at 180 mg, 16.6% at 360 mg, 31.3% at 720 mg.

**The high dose loses proportionally more bioavailability to autoinduction than the low dose.**
That mechanism, not absorption saturation alone, is part of why the exposures converge.
