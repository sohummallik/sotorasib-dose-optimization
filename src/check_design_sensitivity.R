# Recalculate the hypothetical design scenarios with base R.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
root <- dirname(dirname(normalizePath(script)))
assumed <- read.csv(file.path(root, "data", "design_scenarios.csv"))
python <- read.csv(file.path(root, "results", "design_sensitivity.csv"))
if (anyDuplicated(assumed$scenario_id) || anyDuplicated(python$scenario_id) ||
    !setequal(assumed$scenario_id, python$scenario_id)) stop("Scenario identifiers differ.")
python <- python[match(assumed$scenario_id, python$scenario_id), ]
input_columns <- setdiff(names(assumed), "scenario_id")
if (!isTRUE(all.equal(unname(as.matrix(assumed[input_columns])),
                     unname(as.matrix(python[input_columns])), tolerance = 1e-12))) {
  stop("Python output does not match the stated design assumptions.")
}
rows <- lapply(seq_len(nrow(assumed)), function(i) {
  s <- assumed[i, ]
  design <- power.prop.test(p1 = s$comparison_orr, p2 = s$reference_orr,
                           sig.level = s$alpha_two_sided, power = s$target_power,
                           alternative = "two.sided", strict = TRUE, tol = 1e-10)
  n <- ceiling(design$n)
  power_at <- function(size) power.prop.test(n = size, p1 = s$comparison_orr,
      p2 = s$reference_orr, sig.level = s$alpha_two_sided,
      alternative = "two.sided", strict = TRUE)$power
  selected <- power_at(n)
  lower <- power_at(n - 1)
  if (selected < s$target_power || lower >= s$target_power) stop("Invalid integer sample size.")
  if (n != python$n_per_arm[i] || 2 * n != python$n_total[i] ||
      abs(selected - python$power_at_selected_n[i]) > 1e-8 ||
      abs(lower - python$power_at_one_fewer_per_arm[i]) > 1e-8) {
    stop("R/Python design results differ for ", s$scenario_id)
  }
  data.frame(scenario_id = s$scenario_id, n_per_arm = n, n_total = 2 * n,
             power_at_selected_n = selected, power_at_one_fewer_per_arm = lower)
})
precision <- read.csv(file.path(root, "results", "design_precision.csv"))
half_width <- qnorm((1 + precision$confidence_level) / 2) * sqrt(
    (precision$comparison_orr * (1 - precision$comparison_orr) +
     precision$reference_orr * (1 - precision$reference_orr)) / precision$n_per_arm) * 100
if (any(!is.finite(half_width)) ||
    max(abs(half_width - precision$approximate_half_width_percentage_points)) > 1e-10) {
  stop("R/Python precision approximations differ.")
}
write.csv(do.call(rbind, rows), file.path(root, "results", "r_design_sensitivity.csv"), row.names = FALSE)
cat(sprintf("R/Python agreement: %d hypothetical design scenarios and %d precision examples.\n",
            length(rows), nrow(precision)))
