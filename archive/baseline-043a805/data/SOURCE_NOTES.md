# Source provenance: Table 35 transcription

**Source document:** FDA Center for Drug Evaluation and Research, NDA/BLA Multi-disciplinary
Review and Evaluation, NDA 214665, LUMAKRAS (sotorasib). Reference ID 4803204. 269 pages.
Retrieved from Drugs@FDA.

**Table 35 location:** PDF pages 237, 238, 239 (document-numbered pages 237, 238, 239).
Table caption appears on PDF page 236. Not redacted.

**FDA's stated source for Table 35:** "Applicant's PopPK report, Table 1, Page 4".
Table 35 is FDA's reproduction of Amgen's own parameter table, not an independent FDA fit.

**Table 33** (baseline covariate distributions): PDF pages 233 to 235. Source: Applicant's
PopPK report, Table 5, Page 41.

**Figure 18** (covariate map): PDF page 239. Confirms covariate-to-parameter assignment:
ALB, RACE, ECOG_BL, TUM_BS_CAT, SEX -> CL/F; SEX -> V2/F; PPI -> F1; HIGHFAT -> F1 and KA.

**Figure 26** (exposure vs ORR): PDF page 248.
**Figure 29** (exposure vs Grade 3+ AEs): PDF page 250.

---

## Structural model, as described on PDF page 231

- Two-compartment disposition
- Three transit compartments for absorption
- Relative bioavailability (F1) parameterized separately **by dose group**, with 960 mg as reference
- Enzyme induction modeled as an exponential function with first-order rate coefficient KIND_F
- Combined additive and proportional residual error

Autoinduction magnitude, stated in text on PDF page 231: reaches steady state in 2 to 3 weeks,
associated with a 35% decrease in relative bioavailability and a 91% increase in clearance.

Cross-check: CLSS / CLBS = 41.3 / 21.6 = 1.91. Consistent with the stated 91% increase.
Cross-check: F1SS_DG5 / F1BS_DG5 = 1.00 / 1.53 = 0.654, a 35% decrease. Consistent.

Induction half-time = ln(2) / KIND_F = 0.693 / 0.00845 = 82 hours = 3.4 days.
Approximately 4 to 5 induction half-lives = 14 to 17 days, consistent with "2 to 3 weeks".

---
