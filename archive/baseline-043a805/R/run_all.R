## Runs the reported R configuration in order, from the repository root.
## The reduced SAEM run uses 60 subjects and 100/100 iterations; runtime varies.
stopifnot(dir.exists("results"))
source("R/01_model_validate.R")
source("R/03_population_240_vs_960.R")
source("R/02c_simulate_and_recover_diag.R")
cat("\nAll R components complete.\n")
