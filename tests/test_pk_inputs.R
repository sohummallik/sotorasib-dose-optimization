# Run from any directory: Rscript tests/test_pk_inputs.R
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
root <- dirname(dirname(normalizePath(script)))
source(file.path(root, "src/pk_parameters.R"))
parameters <- read.csv(file.path(root, "data/poppk_parameters_final.csv"))
p <- setNames(parameters$final_estimate, parameters$parameter)
checks <- 0L
check <- function(condition) {
  if (!isTRUE(condition)) stop("PK input check failed.")
  checks <<- checks + 1L
}
reject <- function(overrides, pattern, parameters = p) {
  error <- tryCatch({ subject_parameters(parameters, overrides); NULL },
                    error = function(e) conditionMessage(e))
  check(!is.null(error) && grepl(pattern, error, fixed = TRUE))
}

# Typos, missing values, and numeric indicators must not silently become reference values.
for (value in list(99, TRUE, "1", NA, NULL, c(0, 1))) reject(list(ecog = value), "ecog")
for (value in list("female", NA, NULL, c("M", "F"))) reject(list(sex = value), "sex")
reject(list(tumor = "small"), "tumor")
reject(list(race = "unknown"), "race")
for (name in c("ppi", "highfat")) {
  for (value in list(0, 1, "FALSE", NA, NULL, c(TRUE, FALSE))) {
    reject(setNames(list(value), name), name)
  }
}
for (value in list(0, -1, NA, Inf, TRUE, "40", NULL, c(30, 40))) {
  reject(list(alb_gL = value), "alb_gL")
}
reject(list(alb_gL = 1), "Calculated PK parameters")
reject(list(sexx = "F"), "supported fields")
reject(list("F"), "supported fields")
reject(list(sex = "F", sex = "M"), "supported fields")
reject(list(), "Calculated PK parameters", replace(p, "CLSS", -1))
reject(list(), "Calculated PK parameters", replace(p, "V2", 0))
reject(list(), "Calculated PK parameters", replace(p, "K24", Inf))

# Reference values and every supported category remain usable.
reference <- subject_parameters(p)
check(identical(unname(reference[c("clbs", "clss", "v2")]),
                unname(p[c("CLBS", "CLSS", "V2")])))
accepted <- c(lapply(0:2, function(value) list(ecog = value)),
              lapply(c("M", "F"), function(value) list(sex = value)),
              lapply(c(">70", "<=70", "healthy"), function(value) list(tumor = value)),
              lapply(c("Caucasian", "Asian", "Black", "Other", "NHPI", "AIAN"),
                     function(value) list(race = value)),
              list(list(ppi = TRUE), list(highfat = TRUE),
                   list(tumor = "healthy", ecog = 0)))
for (overrides in accepted) {
  values <- subject_parameters(p, overrides)
  check(length(values) == length(reference) && all(is.finite(values) & values > 0))
}
# A lower albumin value changes clearance, leaving relative bioavailability unchanged.
low_albumin <- subject_parameters(p, list(alb_gL = 30))
check(abs(low_albumin[["clss"]] / reference[["clss"]] - 0.704) < 1e-12)
check(identical(low_albumin[c("fbs", "fss")], reference[c("fbs", "fss")]))
cat(sprintf("R PK input checks passed: %d assertions.\n", checks))
