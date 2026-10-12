# Run: Rscript src/check_statistics.R [data_dir] [results_dir]
# FDA Figure 26 provides four grouped counts, not individual exposures.
# These fits cannot establish a causal exposure-response relationship.

script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
root <- dirname(dirname(normalizePath(script)))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2) stop("Supply at most data_dir and results_dir.")
data_dir <- if (length(args) >= 1) args[1] else file.path(root, "data")
results_dir <- if (length(args) >= 2) args[2] else file.path(root, "results")

data <- read.csv(file.path(data_dir, "fda_figure26_quartiles.csv"))
required <- c("quartile", "auc_lower", "auc_upper", "auc_label", "responders", "total")
if (!all(required %in% names(data))) stop("Missing Figure 26 data columns.")
if (nrow(data) != 4 || !all(vapply(data[required], is.numeric, logical(1))) ||
    any(!is.finite(as.matrix(data[required])))) {
  stop("Expected four groups with finite numeric data.")
}
data <- data[order(data$quartile), ]
counts <- as.matrix(data[c("quartile", "responders", "total")])
if (!identical(as.integer(data$quartile), 1:4) || any(counts != floor(counts)) ||
    any(data$total <= 0) || any(data$responders < 0 | data$responders > data$total)) {
  stop("Invalid quartile numbers or binomial counts.")
}
if (any(data$auc_lower <= 0 | data$auc_upper <= data$auc_lower) ||
    any(data$auc_label < data$auc_lower | data$auc_label > data$auc_upper) ||
    any(diff(data$auc_label) <= 0) ||
    any(data$auc_lower[-1] != data$auc_upper[-nrow(data)])) {
  stop("Invalid AUC ranges or displayed coordinates.")
}

coordinates <- list(
  displayed_labels = data$auc_label,
  geometric_range_midpoints = sqrt(as.numeric(data$auc_lower) * data$auc_upper),
  arithmetic_range_midpoints = (data$auc_lower + data$auc_upper) / 2
)
fit_grouped <- function(auc) {
  fit <- glm(cbind(responders, total - responders) ~ log(auc),
             data = data, family = binomial(),
             control = glm.control(epsilon = 1e-12, maxit = 100))
  if (!fit$converged) stop("Grouped binomial fit did not converge.")
  coefficients <- coef(summary(fit))
  slope <- coefficients[2, "Estimate"]
  se <- coefficients[2, "Std. Error"]
  # Match the Python calculation's 1.96 Wald multiplier. This interval does
  # not include coordinate uncertainty, within-group variation or confounding.
  interval <- exp((slope + c(-1.96, 1.96) * se) * log(2))
  c(intercept = coefficients[1, "Estimate"], slope_per_log_AUC = slope,
    slope_SE = se, OR_per_doubling = exp(slope * log(2)),
    conditional_95CI_lower = interval[1], conditional_95CI_upper = interval[2],
    setNames(auc, paste0("coordinate_q", 1:4)))
}
fits <- as.data.frame(do.call(rbind, lapply(coordinates, fit_grouped)))
fits <- data.frame(coordinate_scheme = names(coordinates), fits, row.names = NULL)

# Compare saved Python results before writing any R outputs.
python_file <- file.path(results_dir, "grouped_coordinate_sensitivity.csv")
if (!file.exists(python_file)) stop("Run python src/run_all.py first: missing ", python_file)
python <- read.csv(python_file)
if (!setequal(names(python), names(fits)) || nrow(python) != nrow(fits) ||
    anyDuplicated(python$coordinate_scheme) ||
    !setequal(python$coordinate_scheme, fits$coordinate_scheme)) {
  stop("Python grouped results have unexpected columns or coordinate schemes.")
}
numeric_columns <- setdiff(names(fits), "coordinate_scheme")
python <- python[match(fits$coordinate_scheme, python$coordinate_scheme), ]
if (!all(vapply(python[numeric_columns], is.numeric, logical(1))) ||
    any(!is.finite(as.matrix(python[numeric_columns])))) {
  stop("Python grouped results contain invalid numeric values.")
}
difference <- abs(as.matrix(fits[numeric_columns]) - as.matrix(python[numeric_columns]))
if (any(!is.finite(difference)) || any(difference > 1e-6)) {
  stop("R/Python grouped results disagree (absolute tolerance 1e-6); regenerate or inspect Python results.")
}

# Exact binomial intervals describe response in each reported group.
intervals <- t(vapply(seq_len(nrow(data)), function(i) {
  unname(binom.test(data$responders[i], data$total[i], conf.level = 0.95)$conf.int)
}, numeric(2)))
response <- data.frame(data[c("quartile", "responders", "total")],
                       response_proportion = data$responders / data$total,
                       exact_95CI_lower = intervals[, 1], exact_95CI_upper = intervals[, 2])

# Hypothetical superiority design; not an estimate of the original trial plan.
# Equal allocation, two-sided normal approximation, no continuity correction
# or attrition allowance. strict=TRUE counts rejection in either tail.
p1 <- 0.35
p2 <- 0.25
alpha <- 0.05
target <- 0.80
design <- power.prop.test(p1 = p1, p2 = p2, sig.level = alpha, power = target,
                          alternative = "two.sided", strict = TRUE, tol = 1e-10)
n <- ceiling(design$n)
power_at <- function(size) power.prop.test(n = size, p1 = p1, p2 = p2,
                                          sig.level = alpha, alternative = "two.sided",
                                          strict = TRUE)$power
if (power_at(n) < target || power_at(n - 1) >= target) stop("Sample-size rounding failed.")
scenario <- data.frame(assumed_ORR_arm_1 = p1, assumed_ORR_arm_2 = p2,
                       alpha_two_sided = alpha, target_power = target,
                       continuous_n_per_arm = design$n, n_per_arm = n, n_total = 2 * n,
                       power_at_selected_n = power_at(n),
                       power_at_one_fewer_per_arm = power_at(n - 1))

write.csv(fits, file.path(results_dir, "r_grouped_coordinate_sensitivity.csv"), row.names = FALSE)
write.csv(response, file.path(results_dir, "r_group_response_intervals.csv"), row.names = FALSE)
write.csv(scenario, file.path(results_dir, "r_design_scenario.csv"), row.names = FALSE)
cat(sprintf("R/Python grouped fits agree: maximum absolute difference %.3g.\n", max(difference)))
cat(sprintf("Exact response intervals calculated for %d patients in four groups.\n", sum(data$total)))
cat(sprintf("Hypothetical design: %d per arm, %d total; power %.6f.\n", n, 2 * n, power_at(n)))
