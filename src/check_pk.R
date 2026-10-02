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

subject_parameters <- function(overrides = list()) {
  subject <- modifyList(list(ecog = 1, sex = "M", tumor = ">70",
                             race = "Caucasian", alb_gL = 40,
                             ppi = FALSE, highfat = FALSE), overrides)
  cl_multiplier <- 1 + p[["CLALB1"]] * (subject$alb_gL - 40)
  if (subject$ecog == 0) cl_multiplier <- cl_multiplier * (1 + p[["CLECOGBL1"]])
  if (subject$ecog == 2) cl_multiplier <- cl_multiplier * (1 + p[["CLECOGBL2"]])
  if (subject$sex == "F") cl_multiplier <- cl_multiplier * (1 + p[["CLSEX1"]])
  if (subject$tumor == "<=70") cl_multiplier <- cl_multiplier * (1 + p[["CLTUM_BS_CAT1"]])
  if (subject$tumor == "healthy") cl_multiplier <- cl_multiplier * (1 + p[["CLTUM_BS_CAT2"]])
  race_parameter <- c(Asian = "CLRACE1", Black = "CLRACE2", Other = "CLRACE3",
                      NHPI = "CLRACE4", AIAN = "CLRACE5")
  if (subject$race %in% names(race_parameter)) {
    cl_multiplier <- cl_multiplier * (1 + p[[race_parameter[[subject$race]]]])
  }
  f_multiplier <- 1
  if (subject$ppi) f_multiplier <- f_multiplier * (1 + p[["F1PPI1"]])
  if (subject$highfat) f_multiplier <- f_multiplier * (1 + p[["F1HIGHFAT1"]])
  c(ka = p[["Ka"]] * if (subject$highfat) 1 + p[["KAHIGHFAT1"]] else 1,
    v2 = p[["V2"]] * if (subject$sex == "F") 1 + p[["S2SEX1"]] else 1,
    clbs = p[["CLBS"]] * cl_multiplier, clss = p[["CLSS"]] * cl_multiplier,
    kind = p[["KIND_F"]], k24 = p[["K24"]], k42 = p[["K42"]],
    fbs = p[["F1BS_DG5"]] * f_multiplier, fss = p[["F1SS_DG5"]] * f_multiplier)
}

simulate_subject <- function(overrides = list(), grid_h = 0.03125) {
  # Amounts are ng, volume is L, and cp is ng/mL, matching the Python units.
  events <- et(amt = 960 * 1e6, ii = 24, addl = 29, cmt = "tr1")
  events <- et(events, seq(0, 720, by = grid_h))
  result <- as.data.frame(rxSolve(model, subject_parameters(overrides), events,
                                 method = "liblsoda", rtol = 1e-8, atol = 1e-6,
                                 addDosing = TRUE))
  if (any(!is.finite(result$cp))) stop("Nonfinite concentration")
  # Oral dosing changes the transit compartment; central concentration is continuous.
  duplicate_times <- unique(result$time[duplicated(result$time)])
  if (!all(seq(0, 696, by = 24) %in% duplicate_times)) {
    stop("Missing before/after dose records")
  }
  for (time in duplicate_times) {
    values <- result$cp[result$time == time]
    if (diff(range(values)) > 1e-8 * max(1, abs(values))) {
      stop("Central concentration changed across an oral dose")
    }
  }
  result <- result[!duplicated(result$time), ]
  interval <- result[result$time >= 696 & result$time <= 720, ]
  if (interval$time[1] != 696 || tail(interval$time, 1) != 720 ||
      any(diff(interval$time) <= 0) ||
      abs(sum(diff(interval$time)) - 24) > 1e-10) {
    stop("Incomplete final dosing interval")
  }
  c(AUC_tau = sum(diff(interval$time) *
                   (head(interval$cp, -1) + tail(interval$cp, -1)) / 2),
    Cmax = max(interval$cp))
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
cat(sprintf("Ten contrasts, twenty ratios checked at 960 mg daily on day 30.\n"))
cat(sprintf("Maximum published-point error: %.6f%%\n", published_error))
cat(sprintf("Maximum R/Python difference: %.6f%%\n", python_difference))
cat(sprintf("Maximum grid-refinement change: %.6f%%\n", grid_change))
cat("Full interval endpoints and central concentration continuity checked.\n")
cat("R ", as.character(getRversion()), "; rxode2 ", as.character(packageVersion("rxode2")),
    "; jsonlite ", as.character(packageVersion("jsonlite")), "\n", sep = "")
