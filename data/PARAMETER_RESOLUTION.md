# Model parameter notes

The PK checks use the published estimates from [Nagase et al. AAPS J. 2025;27:26](https://doi.org/10.1208/s12248-024-01013-6), Table II and Figures 1 and 3. They use deterministic scenarios at 960 mg once daily. Several details need interpretation because the original model code was not available.

## Random effects

All random effects are set to zero in the current checks. The code supports lognormal effects: parameter_i = typical × exp(eta_i). This interpretation is consistent with the reported variance/CV and covariance/correlation relationships, with CV calculated as sqrt(exp(variance) − 1) × 100. The printed random-effect equation is ambiguous, so the original implementation cannot be confirmed from that equation alone.

## Albumin

The clearance multiplier is 1 + 0.0296 × (albumin − 40), with albumin in g/L. This approximately reproduces the published low-albumin contrast. The reported unit and other descriptions are not fully consistent. Agreement with that contrast supports this interpretation, but does not confirm an error in the paper or recover the original model code. The linear relationship can give physiologically implausible values outside its range and should not be extrapolated beyond the described patients.

## Absorption and transit time

Figure 1 shows dose entering TR1, then passing through TR2 and TR3 to the central compartment, with a shared Ka. The implementation follows that chain, giving a mean pre-central transit time of 3/Ka = 0.381 hours. The text instead gives MTT = (n + 1)/Ka with n = 3, implying 0.508 hours. This discrepancy remains unresolved. The high-fat transit-time description is also not fully reconciled; the high-fat Cmax contrast has the largest difference from its published point value.

## Dose groups

The low-dose relative bioavailability (F) estimates pool 120, 180, and 240 mg and do not provide a separate estimate for the randomized 240-mg arm. The available data do not show that every low-dose record followed a toxicity-related reduction, or whether any selection bias would raise or lower the estimate. The 480-mg terms came from a twice-daily setting. The current pipeline evaluates neither low-dose projections nor 480-mg once-daily dosing.

## Uncertainty

A coefficient interval crossing zero does not establish nonidentifiability. Small covariate subgroups warrant caution, but do not establish a particular artifact or selection effect. Residual-error parameters are retained in the source table but are not used in the deterministic checks. Marginal parameter intervals cannot replace the full fixed-effect covariance matrix. Agreement with 20 point ratios does not validate residual variability, absolute clinical predictions, or extrapolation to another dose.

The [source notes](SOURCE_NOTES.md) describe the separate randomized trial article and the PK details it leaves unresolved.
