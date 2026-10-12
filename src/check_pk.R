# Check the ten published PK contrasts in R against the current Python results.
# These are deterministic point comparisons, not a fit to patient data.
suppressPackageStartupMessages({
  library(rxode2)
  library(jsonlite)
})

script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (length(script) != 1L) stop("Run with Rscript src/check_pk.R")
root <- dirname(dirname(normalizePath(script)))
parameters <- read.csv(file.path(root, "data/poppk_parameters_final.csv"))
p <- setNames(parameters$final_estimate, parameters$parameter)
targets <- fromJSON(file.path(root, "data/covariate_targets.json"),
                    simplifyVector = FALSE)$contrasts
python <- read.csv(file.path(root, "results/pk_covariate_checks.csv"))

# Keep compiled model files outside the repository.
model_dir <- tempfile("sotorasib-rxode2-")
dir.create(model_dir)
model <- rxode2({
  cl <- clss + (clbs - clss) * exp(-kind * t)
  f(tr1) <- fss + (fbs - fss) * exp(-kind * t)
  d/dt(tr1) <- -ka * tr1
  d/dt(tr2) <- ka * tr1 - ka * tr2
  d/dt(tr3) <- ka * tr2 - ka * tr3
  d/dt(centr) <- ka * tr3 - (cl / v2) * centr - k24 * centr + k42 * peri
  d/dt(peri) <- k24 * centr - k42 * peri
  cp <- centr / v2 / 1000
}, modName = "sotorasib_pk_check", wd = model_dir)

source(file.path(root, "src/pk_parameters.R"))

simulate_subject <- function(overrides = list(), grid_h = 0.03125,
                             n_doses = 30L, freeze_induction = FALSE) {
  # Amounts are ng, volume is L, and cp is ng/mL, matching the Python units.
  individual <- subject_parameters(p, overrides)
  if (freeze_induction) {
    individual[c("clbs", "fbs")] <- individual[c("clss", "fss")]
  }
  start <- (n_doses - 1L) * 24
  end <- n_doses * 24
  events <- et(amt = 960 * 1e6, ii = 24, addl = n_doses - 1L, cmt = "tr1")
  events <- et(events, seq(0, end, by = grid_h))
  result <- as.data.frame(rxSolve(model, individual, events,
                                 method = "liblsoda", rtol = 1e-8, atol = 1e-6,
                                 addDosing = TRUE))
  if (any(!is.finite(result$cp))) stop("Nonfinite concentration")
  # Oral dosing changes the transit compartment; central concentration is continuous.
  duplicate_times <- unique(result$time[duplicated(result$time)])
  if (!all(seq(0, start, by = 24) %in% duplicate_times)) {
    stop("Missing before/after dose records")
  }
  for (time in duplicate_times) {
    values <- result$cp[result$time == time]
    if (diff(range(values)) > 1e-8 * max(1, abs(values))) {
      stop("Central concentration changed across an oral dose")
    }
  }
  result <- result[!duplicated(result$time), ]
  interval <- result[result$time >= start & result$time <= end, ]
  if (interval$time[1] != start || tail(interval$time, 1) != end ||
      any(diff(interval$time) <= 0) ||
      abs(sum(diff(interval$time)) - 24) > 1e-10) {
    stop("Incomplete final dosing interval")
  }
  c(AUC_tau = sum(diff(interval$time) *
                   (head(interval$cp, -1) + tail(interval$cp, -1)) / 2),
    Cmax = max(interval$cp))
}

# With induction held fixed and accumulation complete, AUC_tau = dose * F / CL.
# The factor of 1000 converts mg/L to ng/mL; this checks numerical units as well.
frozen <- simulate_subject(n_doses = 60L, freeze_induction = TRUE)
reference_parameters <- subject_parameters(p)
exact_auc <- 960 * reference_parameters[["fss"]] / reference_parameters[["clss"]] * 1000
mass_balance <- data.frame(calculated_AUC_ng_h_mL = frozen[["AUC_tau"]],
                           exact_AUC_ng_h_mL = exact_auc,
                           relative_error = frozen[["AUC_tau"]] / exact_auc - 1)
if (abs(mass_balance$relative_error) >= 0.0005) {
  stop("Fixed-parameter AUC error exceeds 0.05%.")
}

reference <- simulate_subject()
reference_fine <- simulate_subject(grid_h = 0.015625)
rows <- lapply(targets, function(target) {
  estimate <- simulate_subject(target$overrides) / reference
  fine <- simulate_subject(target$overrides, grid_h = 0.015625) / reference_fine
  match <- python[python$contrast == target$contrast, ]
  if (nrow(match) != 1) stop("Missing or duplicate Python contrast: ", target$contrast)
  row <- data.frame(contrast = target$contrast)
  for (metric in c("AUC_tau", "Cmax")) {
    python_value <- match[[paste0(metric, "_ratio_model")]]
    if (match[[paste0(metric, "_ratio_published")]] != target[[metric]]) {
      stop("Python and R source targets disagree")
    }
    row[[paste0(metric, "_ratio_R")]] <- estimate[[metric]]
    row[[paste0(metric, "_ratio_published")]] <- target[[metric]]
    row[[paste0(metric, "_ratio_Python")]] <- python_value
    row[[paste0(metric, "_published_error_pct")]] <- 100 * (estimate[[metric]] / target[[metric]] - 1)
    row[[paste0(metric, "_Python_difference_pct")]] <- 100 * (estimate[[metric]] / python_value - 1)
    row[[paste0(metric, "_grid_change_pct")]] <- 100 * (fine[[metric]] / estimate[[metric]] - 1)
  }
  row
})
output <- do.call(rbind, rows)
max_abs <- function(suffix) max(abs(as.matrix(output[endsWith(names(output), suffix)])))
published_error <- max_abs("_published_error_pct")
python_difference <- max_abs("_Python_difference_pct")
grid_change <- max_abs("_grid_change_pct")

# Numerical tolerances are separate from the 5% published-point comparison.
if (published_error > 5) stop("Published-point discrepancy exceeds 5%")
if (python_difference > 0.01) stop("R/Python discrepancy exceeds 0.01%")
if (grid_change > 0.1) stop("Grid-refinement change exceeds 0.1%")
write.csv(output, file.path(root, "results/r_pk_covariate_checks.csv"), row.names = FALSE)
write.csv(mass_balance, file.path(root, "results/r_pk_numerical_checks.csv"), row.names = FALSE)
cat(sprintf("Ten contrasts, twenty ratios checked at 960 mg daily on day 30.\n"))
cat(sprintf("Maximum published-point error: %.6f%%\n", published_error))
cat(sprintf("Maximum R/Python difference: %.6f%%\n", python_difference))
cat(sprintf("Maximum grid-refinement change: %.6f%%\n", grid_change))
cat("Full interval endpoints and central concentration continuity checked.\n")
cat(sprintf("Fixed-parameter AUC relative error: %.6f%%\n", 100 * mass_balance$relative_error))
cat("R ", as.character(getRversion()), "; rxode2 ", as.character(packageVersion("rxode2")),
    "; jsonlite ", as.character(packageVersion("jsonlite")), "\n", sep = "")
