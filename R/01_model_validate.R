## 01_model_validate.R
## rxode2 implementation of the sotorasib popPK model, validated against
## (a) the ten published covariate exposure ratios in Nagase et al. AAPS J
##     2025;27:26 Fig. 3, and
## (b) the independent Python implementation in src/model.py.
##
## RUN FROM THE REPOSITORY ROOT:   source("R/01_model_validate.R")
## Requires: R/00_install.R run once on a machine with internet access.
##
## Executed September 2026; see R/README.md for verification results.

suppressPackageStartupMessages({ library(rxode2); library(dplyr) })
if (!dir.exists("results")) stop("Run from the repository root (setwd to sotorasib-dose-optimization).")

## ---- Parameters: single source of truth is data/poppk_parameters_final.csv
prm <- read.csv("data/poppk_parameters_final.csv", stringsAsFactors = FALSE)
P   <- setNames(prm$final_estimate, prm$parameter)

## ---- Structural model ----------------------------------------------------
## Units: dose in mg. V2 in L. centr/v2 = mg/L = ug/mL; x1000 -> ng/mL.
## Covariate multipliers (CLMULT, V2MULT, F1MULT, KAMULT) and the dose-group
## bioavailability terms (F1BS_IN, F1SS_IN) are computed in R per subject and
## passed as covariates, exactly as a NONMEM dataset would carry them.
mod <- rxode2({
  ka   <- TVKA   * KAMULT * exp(eta.ka)
  v2   <- TVV2   * V2MULT * exp(eta.v2)
  clbs <- TVCLBS * CLMULT * exp(eta.cl)
  clss <- TVCLSS * CLMULT * exp(eta.cl)
  cl   <- clss + (clbs - clss) * exp(-TVKIND * t)
  f1   <- (F1SS_IN + (F1BS_IN - F1SS_IN) * exp(-TVKIND * t)) * F1MULT
  f(tr1) <- f1
  d/dt(tr1)   <- -ka * tr1
  d/dt(tr2)   <-  ka * tr1 - ka * tr2
  d/dt(tr3)   <-  ka * tr2 - ka * tr3
  d/dt(centr) <-  ka * tr3 - (cl / v2) * centr - TVK24 * centr + TVK42 * peri
  d/dt(peri)  <-  TVK24 * centr - TVK42 * peri
  cp <- centr / v2 * 1000
})

theta <- c(TVKA = P["Ka"], TVV2 = P["V2"], TVCLBS = P["CLBS"], TVCLSS = P["CLSS"],
           TVK24 = P["K24"], TVK42 = P["K42"], TVKIND = P["KIND_F"])
names(theta) <- c("TVKA","TVV2","TVCLBS","TVCLSS","TVK24","TVK42","TVKIND")

dose_group <- function(d) {
  if (d %in% c(120,180,240)) "DG1" else if (d == 360) "DG2" else if (d == 480) "DG3"
  else if (d %in% c(600,720,840)) "DG4" else if (d == 960) "DG5"
  else stop("dose not in a published dose group")
}

## Mirrors individual_params() in src/model.py exactly.
build_subject <- function(dose_mg, ecog = 1, sex = "M", tumor = ">70", race = "Caucasian",
                          alb_gL = 40, ppi = FALSE, highfat = FALSE,
                          f1bs_override = NA, f1ss_override = NA) {
  clm <- 1
  if (ecog == 0) clm <- clm * (1 + P["CLECOGBL1"])
  if (ecog == 2) clm <- clm * (1 + P["CLECOGBL2"])
  if (sex == "F") clm <- clm * (1 + P["CLSEX1"])
  if (tumor == "<=70")    clm <- clm * (1 + P["CLTUM_BS_CAT1"])
  if (tumor == "healthy") clm <- clm * (1 + P["CLTUM_BS_CAT2"])
  race_par <- c(Asian="CLRACE1", Black="CLRACE2", Other="CLRACE3", NHPI="CLRACE4", AIAN="CLRACE5")
  if (race %in% names(race_par)) clm <- clm * (1 + P[race_par[[race]]])
  clm <- clm * (1 + P["CLALB1"] * (alb_gL - 40))        # proportional per g/L, ref 40 g/L
  v2m <- if (sex == "F") 1 + P["S2SEX1"] else 1
  f1m <- 1
  if (ppi)     f1m <- f1m * (1 + P["F1PPI1"])
  if (highfat) f1m <- f1m * (1 + P["F1HIGHFAT1"])
  kam <- if (highfat) 1 + P["KAHIGHFAT1"] else 1
  dg  <- dose_group(dose_mg)
  c(CLMULT = unname(clm), V2MULT = unname(v2m), F1MULT = unname(f1m), KAMULT = unname(kam),
    F1BS_IN = if (is.na(f1bs_override)) unname(P[paste0("F1BS_", dg)]) else f1bs_override,
    F1SS_IN = if (is.na(f1ss_override)) unname(P[paste0("F1SS_", dg)]) else f1ss_override)
}

simulate_typical <- function(dose_mg, subj, n_doses = 30, dt = 0.25) {
  ev <- et(amt = dose_mg, ii = 24, addl = n_doses - 1, cmt = "tr1") %>%
        et(seq(0, n_doses * 24, by = dt))
  par <- c(theta, subj, eta.ka = 0, eta.v2 = 0, eta.cl = 0)
  as.data.frame(rxSolve(mod, par, ev))
}

interval_metrics <- function(sim, dose_no, tau = 24) {
  d <- sim[sim$time >= (dose_no - 1) * tau & sim$time <= dose_no * tau, ]
  auc <- sum(diff(d$time) * (head(d$cp, -1) + tail(d$cp, -1)) / 2)
  c(AUC_tau = auc, Cmax = max(d$cp), Cmin = tail(d$cp, 1))
}

## ---- Ten published contrasts (AAPS Fig. 3) + Python cross-check ----------
targets <- tribble(
  ~contrast,                            ~ecog, ~sex, ~tumor,   ~race,       ~alb, ~ppi,  ~hf,   ~cmax_pub, ~auc_pub, ~cmax_py, ~auc_py,
  "Low albumin (30 g/L) vs normal (40)",   1, "M", ">70",     "Caucasian", 30,   FALSE, FALSE, 1.075, 1.414, 1.079, 1.419,
  "ECOG 0 vs ECOG 1",                      0, "M", ">70",     "Caucasian", 40,   FALSE, FALSE, 0.977, 0.868, 0.974, 0.869,
  "ECOG 2 vs ECOG 1",                      2, "M", ">70",     "Caucasian", 40,   FALSE, FALSE, 1.017, 1.083, 1.016, 1.081,
  "Asian vs Caucasian",                    1, "M", ">70",     "Asian",     40,   FALSE, FALSE, 0.966, 0.817, 0.966, 0.820,
  "Black vs Caucasian",                    1, "M", ">70",     "Black",     40,   FALSE, FALSE, 0.988, 0.934, 0.986, 0.932,
  "Female vs Male",                        1, "F", ">70",     "Caucasian", 40,   FALSE, FALSE, 1.242, 1.190, 1.237, 1.192,
  "Tumor <=70 mm vs >70 mm",               1, "M", "<=70",    "Caucasian", 40,   FALSE, FALSE, 0.982, 0.892, 0.979, 0.894,
  "High-fat meal vs fasted",               1, "M", ">70",     "Caucasian", 40,   FALSE, TRUE,  1.130, 1.327, 1.157, 1.328,
  "PPI vs no PPI",                         1, "M", ">70",     "Caucasian", 40,   TRUE,  FALSE, 0.828, 0.820, 0.822, 0.822,
  "Healthy (0 mm, ECOG 0) vs patients",    0, "M", "healthy", "Caucasian", 40,   FALSE, FALSE, 0.898, 0.534, 0.907, 0.532
)

cat("\nSimulating reference subject (960 mg QD, 30 days) ...\n")
ref <- interval_metrics(simulate_typical(960, build_subject(960)), 30)

out <- lapply(seq_len(nrow(targets)), function(i) {
  r <- targets[i, ]
  s <- build_subject(960, r$ecog, r$sex, r$tumor, r$race, r$alb, r$ppi, r$hf)
  m <- interval_metrics(simulate_typical(960, s), 30)
  data.frame(contrast = r$contrast,
             Cmax_R = round(m["Cmax"] / ref["Cmax"], 3), Cmax_pub = r$cmax_pub, Cmax_py = r$cmax_py,
             AUC_R  = round(m["AUC_tau"] / ref["AUC_tau"], 3), AUC_pub = r$auc_pub, AUC_py = r$auc_py)
}) %>% bind_rows() %>%
  mutate(Cmax_err_vs_pub_pct = round((Cmax_R / Cmax_pub - 1) * 100, 1),
         AUC_err_vs_pub_pct  = round((AUC_R  / AUC_pub  - 1) * 100, 1),
         R_vs_Python_AUC_pct = round((AUC_R  / AUC_py   - 1) * 100, 2),
         pass = abs(Cmax_err_vs_pub_pct) <= 5 & abs(AUC_err_vs_pub_pct) <= 5)

cat("\n=== rxode2 vs published (AAPS Fig. 3) vs Python (src/validate.py) ===\n")
print(out, row.names = FALSE)
cat(sprintf("\nContrasts within 5%% of published: %d / %d\n", sum(out$pass), nrow(out)))
cat(sprintf("Max |R - Python| AUC ratio discrepancy: %.2f%%  (should be well under 1%%)\n",
            max(abs(out$R_vs_Python_AUC_pct))))

## Absolute-scale check, typical subject, 960 mg
sim8 <- simulate_typical(960, build_subject(960), n_doses = 8)
d1 <- interval_metrics(sim8, 1); d8 <- interval_metrics(sim8, 8)
cat(sprintf("\nDay 1 AUC0-24 %.0f h*ng/mL (Python 52084; AAPS model-pred 53600; observed 57300)\n", d1["AUC_tau"]))
cat(sprintf("Day 8 AUC0-24 %.0f h*ng/mL (Python 29749; AAPS model-pred 35300; observed 37200)\n", d8["AUC_tau"]))
cat(sprintf("Day 8 Tmax %.2f h (observed median 1.1 h)\n", sim8$time[sim8$time >= 168 & sim8$time <= 192][which.max(sim8$cp[sim8$time >= 168 & sim8$time <= 192])] - 168))

write.csv(out, "results/validation_covariate_ratios_R.csv", row.names = FALSE)
cat("\nWritten: results/validation_covariate_ratios_R.csv\n")
