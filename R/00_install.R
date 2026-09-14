## 00_install.R -- install dependencies once before running the R pipeline.
## Requires internet access and build tools for compiled packages.
## Installation time depends on the platform.

pkgs <- c("rxode2", "nlmixr2", "lotri", "dplyr", "tibble", "jsonlite", "ggplot2", "vpc")
install.packages(pkgs, repos = "https://cloud.r-project.org")

ok <- vapply(c("rxode2", "nlmixr2", "lotri", "jsonlite"), requireNamespace, logical(1), quietly = TRUE)
if (all(ok)) {
  cat("\nInstalled. rxode2", as.character(packageVersion("rxode2")),
      " nlmixr2", as.character(packageVersion("nlmixr2")), "\n")
  cat("Next:  source('R/01_model_validate.R')  from the repository root.\n")
} else stop("Missing: ", paste(names(ok)[!ok], collapse = ", "))
