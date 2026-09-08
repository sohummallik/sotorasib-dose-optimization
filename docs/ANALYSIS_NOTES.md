# Analysis notes after reading Hochmair 2024 and Popat/Ratain 2024

## 1. Scope relative to Popat and Ratain

**What they did.** Their Figure 1 is a direct reproduction of FDA's exposure-response
figures: ORR by exposure quartile (panels A, B) plus Kaplan-Meier curves for OS and PFS
stratified by exposure quartile (panels C through F), on both AUCtau,ss and Ctrough,ss.
They also assembled the regulatory and commercial argument: the patent filing, the
cost-effectiveness case, the label's failure to address chronic grade 1-2 GI toxicity,
and the burden-of-proof framing.

**What they did NOT do, verified by reading the full text:**

- No power calculation. They mention power once, in passing, regarding CodeBreaK 200's
  OS endpoint. They never compute the sample size CodeBreaK 100 part B would have needed,
  and they never state the power it had.
- No continuous exposure-response modeling. They reproduce FDA's quartile bins as-is.
  No logistic regression on continuous exposure, no covariate adjustment for ECOG,
  tumor size, or albumin.
- No popPK simulation. They use published summary exposures only.
- No analysis of the asymmetric dose-reduction rule (see section 3).

These comparisons define the scope relative to the cited papers, not a systematic claim of novelty. Hochmair reports observed exposure ratios; this repository additionally explores model-based exposure overlap.

---

## 2. Reframe of contribution 1: model versus observation

The naive popPK calculation and the observed data disagree, and the disagreement is
informative rather than embarrassing.

| Source | 240 mg / 960 mg exposure ratio |
|---|---|
| popPK Table II, naive Dose x F1SS | 1.15 |
| Hochmair observed, Day 1 | 0.67 |
| Hochmair observed, Day 8 | 0.77 |

The model overpredicts 240 mg exposure by roughly 1.5-fold. The explanation is the
dose-group lumping already documented in PARAMETER_RESOLUTION.md: F1SS_DG1 pools
120, 180, and 240 mg, is driven mainly by 180 mg data, and is contaminated by
dose-reduction records. Amgen's own discussion concedes the point: the study was designed
on the assumption that 240 mg and 960 mg exposures were similar, and the PK data did not
support that assumption.

**This repository quantifies the timing argument discussed by Popat and Ratain.**

Popat and Ratain note that CodeBreaK 100 part B sampled PK only on Days 1 and 8, and that
steady state is not reached until as late as Day 22. They assert this matters. They do not
quantify it.

The popPK model can. Autoinduction reduces steady-state relative bioavailability by 34.2%
at 960 mg but only about 7.5% at 180 mg. The high dose loses proportionally more. So the
exposure gap should keep narrowing after Day 8.

The observed data already show this happening:

| Timepoint | 960/240 AUC ratio |
|---|---|
| Day 1 | 1.5 |
| Day 8 | 1.3 |

A 13% narrowing in seven days, with roughly two more induction half-times still to run
(Kind = 0.00845/hr, t-half 82 hours).

**Deliverable:** simulate the 960/240 exposure ratio as a function of time from Day 1
through Day 30, calibrated so the Day 1 and Day 8 predictions match the observed 1.5 and
1.3. Report the model-projected steady-state ratio with uncertainty. 

---

## 3. The asymmetric dose-reduction rule: implications for interpretation

From Hochmair Methods, section 2.1, verbatim in substance:
**dose reductions were permitted for the 960 mg group only.**

Neither Amgen's discussion nor Popat and Ratain's commentary analyzes what this does to
the comparison.

**What was actually compared.** Not 960 mg versus 240 mg. Rather:

- Arm A: start at 960 mg, with tolerability-guided downward titration available
- Arm B: fixed 240 mg, interrupt or discontinue only

Arm A and Arm B were not identical dosing interventions from a treatment-management perspective. The 960 mg arm permitted dose reduction, whereas the 240 mg arm did not. This asymmetry complicates interpretation of treatment exposure, dose intensity, and tolerability between the arms. It does not by itself establish that the modification rule favored either arm for efficacy.

**Reported treatment-management and safety outcomes:**

| | 960 mg | 240 mg |
|---|---|---|
| TEAE leading to dose reduction | 17.3% | N/A (not permitted) |
| Hepatotoxicity leading to discontinuation | 5.8% | **10.6%** |
| ALT increased, any grade | 14.4% | **17.3%** |

Hepatotoxicity-driven discontinuation is nearly twice as high in the **low** dose arm.
Amgen's own text explains the mechanism: elevated transaminases were managed by dose
interruption in the 240 mg arm, but by interruption followed by dose reduction in the
960 mg arm. The low-dose arm had one fewer tool for staying on drug.

**Interpretation.** The different dose-modification options complicate comparison of delivered dose intensity, treatment persistence and tolerability. The reported rates do not establish that the rules caused a tolerability or efficacy difference.

**Potential extension (not performed here):** assess delivered dose intensity and compare symmetric versus asymmetric modification strategies under explicit assumptions. This would require additional data and would not by itself identify causal bias.

---

## 4. Baseline imbalance, quantified

Hochmair Table 1 states baseline characteristics were generally balanced with one
exception: history of liver metastasis was 9.6% in the 960 mg arm and 18.1% in the
240 mg arm. Nearly two-fold, in a poor-prognostic factor, favoring the high-dose arm.

This identifies a baseline imbalance that may affect interpretation; aggregate summaries cannot quantify its contribution to the observed efficacy difference.

Note the interaction with the subgroup analysis: Amgen reports that OS favored 960 mg in
patients with liver metastasis (HR 0.58, 95% CI 0.23 to 1.46), a subgroup that was
twice as prevalent in the arm that did worse overall. Wide CI, small n, and a confounded
starting point.

---

## 5. Power, and Amgen's own number

Hochmair Discussion states the trial was not powered for formal hypothesis testing
because that would have required enrollment of more than 600 patients.

This is Amgen's figure, in the primary publication, not a secondary characterization.
The repository estimates 658 total participants for 80% power at 35% versus 25% ORR, with equal allocation and two-sided alpha 0.05. Separately, treating the observed 34/104 versus 26/105 rates as true probabilities gives approximately 24% post hoc power. These are distinct calculations; the latter does not reconstruct Amgen's planning assumptions or add independent evidence about the null hypothesis.

Observed: 34/104 versus 26/105. Stratified ORR difference 6.9%, 90% CI -3.4 to 17.1.
The interval includes zero, and the entire difference is eight patients.

---

## 6. Additional items worth carrying

- **Ratain's absorption arithmetic:** a 300% dose increase produced roughly a 30% exposure
  increase, implying less than 10% of the incremental 720 mg is absorbed. This is a proposed explanation for GI toxicity, not a causal result of this repository.

- **The 240 mg BID proposal:** because absorption is saturable, splitting the dose could alter exposure. This is an untested extension here; the QD-derived model does not validate a BID regimen.

- **CodeBreaK 300 (colorectal) showed the opposite PK result:** 240 mg gave numerically
  higher exposure than 960 mg at Day 15. The cross-trial difference does not establish its cause; populations and settings also differ.

- **Real-world VA data:** dose-reduced patients had PFS HR 0.60 and OS HR 0.42. Confounded
  by immortal time bias, since patients must survive long enough to be reduced. State the
  bias in the same breath as the numbers.

- **Fatal TEAEs:** 6 versus 4 patients. The 960 mg causes listed include assisted suicide,
  hemorrhagic stroke, and vascular rupture.
