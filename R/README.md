# R / rxode2 / nlmixr2 component

The primary analysis was implemented in Python. An independent R implementation using rxode2 and nlmixr2 was subsequently executed in Posit Cloud in September 2026.

## Reported results and execution status

`01_model_validate.R` reproduced the ten published covariate contrasts, with R and Python agreeing within 0.2%. `03_population_240_vs_960.R` produced a GMR of 0.818, compared with Python's 0.819. These checks support the relative-exposure implementation; they do not establish global model validity.

The currently reported estimation results come from `02c_simulate_and_recover_diag.R`: 60 synthetic subjects, diagonal omega, 100 burn-in and 100 estimation iterations, with Day 1/8 sampling. Ka, V2, k24 and k42 recovered within 7%; CLBS, CLSS and Kind had approximately 17–31% negative bias. A single reduced-scale run cannot attribute those biases to sampling design alone; sample size, iteration count, numerical behavior and model assumptions may contribute.

The saved recovery table contains point estimates recorded from the final printed SAEM iteration. The session exited during the covariance/SE calculation, before the fit object and downstream exports persisted. Formal RSE, shrinkage and saved omega recovery are therefore unavailable from that run. See `results/parameter_recovery_nlmixr2_NOTES.txt`. A successful rerun writes a full-precision recovery table and fit object; exact reproduction of the archived rounded values is not guaranteed across package versions.

The full 150-subject implementation, `02_simulate_and_recover.R`, is retained as a higher-complexity alternative with correlated random effects and 300/300 iterations. It is not the source of the reported recovery table. `02b_simulate_and_recover_fast.R` is the 60-subject correlated-omega alternative, also not the source of that table.

## Run order (from the repository root)

```r
source("R/00_install.R") # once, to install dependencies
source("R/run_all.R")
```

The runner executes these scripts in order:

| Script | Purpose |
|---|---|
| `01_model_validate.R` | Check the rxode2 implementation against ten published contrasts and Python results. |
| `03_population_240_vs_960.R` | Simulate 1,000 subjects at each dose using calibrated F1; compare GMR, overlap and CV. |
| `02c_simulate_and_recover_diag.R` | Demonstrate SAEM estimation on 60 synthetic subjects with diagonal omega; export recovery and diagnostics if the fit returns. |

Runtime depends on the machine and package versions. The SAEM estimation and covariance step can take substantially longer than the exposure calculations. The alternative estimation scripts write to the same result paths; use a separate checkout to compare runs without replacing the reported outputs.

## Scope and limitations

This component demonstrates implementation of the published structural model and a synthetic estimation workflow. It does not fit patient-level trial data, which are not public. The reduced-scale estimation uses diagonal random effects rather than the published covariance structure, and a lognormal residual approximation described in the script.

Setting `RICH <- TRUE` adds Day 15 and Day 22 troughs as an exploratory sampling-design comparison. No completed paired sparse-versus-rich experiment is reported here. Repeated simulations with comparable estimation settings would be needed to assess effects on bias, precision and convergence.
