# Operational parameter choices and unresolved reporting ambiguities

Primary source: [Nagase et al., AAPS Journal 2025;27:26](https://doi.org/10.1208/s12248-024-01013-6), Table II and Figures 1/3. Numeric final-estimate values remain unchanged. Retained analyses use 960 mg deterministic scenarios; they do not estimate new parameters or simulate clinical low-dose equivalence.

## Random effects

The operational model uses lognormal effects: parameter_i = typical × exp(eta_i). Table II's reported variance/CV and covariance/correlation relationships are consistent with this interpretation. The variance-to-CV relation is sqrt(exp(variance)-1)×100. The paper's printed random-effect equation is not fully consistent with the customary expression. This is documented as reporting ambiguity; no original control stream was obtained.

Nagase explicitly describes variability/correlation presentation. The FDA's broad table heading should not be used to interpret all parenthetical random-effect values as relative standard errors. Selected numerical agreement does not establish all hidden model conventions. Random effects are zero in the active deterministic checks; none are drawn.

## Albumin

The implemented clearance multiplier is 1+0.0296×(albumin_g/L−40). This reproduces the low-albumin point contrast approximately. The table's reported unit and other descriptions are not fully consistent. It is an operational interpretation supported by that check, not confirmation of an author error or the original model code. It should not be extrapolated beyond the described patient range. A proportional-linear form can become nonphysical outside its domain.

## Absorption/transit time

Figure 1 shows dose entering TR1 followed by TR2, TR3 and the central compartment, with a shared Ka. The implementation follows this three-transfer chain. Its mean pre-central transit time is 3/Ka=0.381 h. The text separately states MTT=(n+1)/Ka with n=3, implying 0.508 h. This internal source inconsistency remains unresolved. A fourth compartment was not added merely to reproduce the reported time. The high-fat transit-time description is likewise not fully reconciled; the largest retained source-point discrepancy remains the high-fat Cmax contrast.

## Dose-group pooling

The low-dose F values pool 120, 180 and 240 mg. The listed assigned levels do not identify a distinct randomized 240-mg estimate in this PK dataset. The data cannot support the old assertion that every low-dose record must be a toxicity-selected reduction, or establish the direction of resulting bias. No low-dose projection remains active. The 480-mg terms arose in a BID setting and are not evaluated as a QD regimen by the active pipeline.

## Sparse covariates, residual error and uncertainty

A coefficient CI crossing zero does not prove nonidentifiability. Very small subgroups warrant caution without claiming a specific artifact or causal selection effect. Residual-error parameters are retained for source provenance but are not used by deterministic point checks. A single marginal parameter CI is not a substitute for full fixed-effect covariance or structural uncertainty. Matching twenty point ratios does not validate residual variability, absolute clinical predictions or dose extrapolation.
