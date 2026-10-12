# Deterministic covariate parameters at the published 960 mg dose group.
subject_parameters <- function(p, overrides = list()) {
  reference <- list(ecog = 1, sex = "M", tumor = ">70", race = "Caucasian",
                    alb_gL = 40, ppi = FALSE, highfat = FALSE)
  if (!is.list(overrides) ||
      (length(overrides) > 0 && (is.null(names(overrides)) ||
       anyNA(names(overrides)) || anyDuplicated(names(overrides)) ||
       !all(names(overrides) %in% names(reference))))) {
    stop("Covariate overrides must be a named list of supported fields.")
  }
  subject <- modifyList(reference, overrides, keep.null = TRUE)
  if (!is.numeric(subject$ecog) || length(subject$ecog) != 1L ||
      is.na(subject$ecog) || !subject$ecog %in% 0:2) {
    stop("ecog must be 0, 1, or 2.")
  }
  allowed <- list(sex = c("M", "F"), tumor = c(">70", "<=70", "healthy"),
                  race = c("Caucasian", "Asian", "Black", "Other", "NHPI", "AIAN"))
  for (name in names(allowed)) {
    value <- subject[[name]]
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !value %in% allowed[[name]]) stop("Unsupported ", name, " category.")
  }
  for (name in c("ppi", "highfat")) {
    value <- subject[[name]]
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(name, " must be TRUE or FALSE.")
    }
  }
  if (!is.numeric(subject$alb_gL) || length(subject$alb_gL) != 1L ||
      !is.finite(subject$alb_gL) || subject$alb_gL <= 0) {
    stop("alb_gL must be finite and positive.")
  }
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
  result <- c(ka = p[["Ka"]] * if (subject$highfat) 1 + p[["KAHIGHFAT1"]] else 1,
    v2 = p[["V2"]] * if (subject$sex == "F") 1 + p[["S2SEX1"]] else 1,
    clbs = p[["CLBS"]] * cl_multiplier, clss = p[["CLSS"]] * cl_multiplier,
    kind = p[["KIND_F"]], k24 = p[["K24"]], k42 = p[["K42"]],
    fbs = p[["F1BS_DG5"]] * f_multiplier, fss = p[["F1SS_DG5"]] * f_multiplier)
  if (any(!is.finite(result)) || any(result <= 0)) {
    stop("Calculated PK parameters must be finite and positive.")
  }
  result
}

