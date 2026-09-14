## 03_population_240_vs_960.R
## Component C in R: steady-state exposure overlap between 240 mg and 960 mg
## in a virtual CodeBreaK 100 part B population, with the published omega
## matrix and the 240 mg bioavailability calibrated to the observed Day 1 /
## Day 8 exposure ratios (results/calibrated_240mg_F1.json).
##
## Parity target from the Python pipeline (results/exposure_overlap.json,
## scenario B):  GMR 240/960 = 0.819, overlap coefficient 0.833, BSV CV ~59%.
##
## RUN FROM THE REPOSITORY ROOT:   source("R/03_population_240_vs_960.R")
## Requires 01_model_validate.R to have been sourced (uses mod, theta, build_subject).

suppressPackageStartupMessages({ library(rxode2); library(dplyr); library(jsonlite) })
if (!exists("mod")) source("R/01_model_validate.R")
set.seed(20260901)

cal <- fromJSON("results/calibrated_240mg_F1.json")
N <- 1000; n_doses <- 30

## Covariate distributions: Hochmair 2024 Table 1 (pooled arms) and AAPS Table I (NSCLC)
covs <- tibble(
  id    = 1:N,
  ecog  = sample(c(0,1,2), N, TRUE, prob = c(0.354, 0.569, 0.077)),
  sex   = sample(c("M","F"), N, TRUE, prob = c(0.545, 0.455)),
  race  = sample(c("Caucasian","Asian"), N, TRUE, prob = c(0.828, 0.172)),
  tumor = sample(c(">70","<=70"), N, TRUE, prob = c(0.475, 0.525)),
  alb   = pmin(pmax(rnorm(N, 37.9, 5.0), 21), 51),   # SD 5 g/L is an ASSUMPTION (source gives median/range only)
  ppi   = runif(N) < 0.37
)

## Same eta draws for both arms (each virtual patient receives both doses)
omega <- lotri({ eta.v2 + eta.cl + eta.ka ~ c(P[["OMEGA_1_1"]], P[["OMEGA_2_1"]], P[["OMEGA_2_2"]],
                                               P[["OMEGA_3_1"]], P[["OMEGA_3_2"]], P[["OMEGA_3_3"]]) })
etas <- as.data.frame(rxRmvn(N, rep(0, 3), omega)); names(etas) <- c("eta.v2","eta.cl","eta.ka")

make_icov <- function(dose_mg, f1bs = NA, f1ss = NA) {
  bind_cols(covs["id"], etas,
            bind_rows(lapply(seq_len(N), function(i) {
              as.data.frame(t(build_subject(dose_mg, covs$ecog[i], covs$sex[i], covs$tumor[i], covs$race[i],
                                            covs$alb[i], covs$ppi[i], FALSE, f1bs, f1ss)))
            })))
}

run_arm <- function(dose_mg, f1bs = NA, f1ss = NA) {
  ev <- et(amt = dose_mg, ii = 24, addl = n_doses - 1, cmt = "tr1") %>%
        et(seq((n_doses - 1) * 24, n_doses * 24, by = 0.25)) %>% et(id = 1:N)
  s <- rxSolve(mod, theta, ev, iCov = make_icov(dose_mg, f1bs, f1ss), returnType = "data.frame")
  s %>% group_by(id) %>%
    summarise(AUC_tau = sum(diff(time) * (head(cp, -1) + tail(cp, -1)) / 2),
              Cmax = max(cp), Cmin = last(cp), .groups = "drop") %>%
    mutate(dose_mg = dose_mg)
}

cat("Simulating 960 mg arm ...\n"); a960 <- run_arm(960)
cat("Simulating 240 mg arm (calibrated F1) ...\n"); a240 <- run_arm(240, cal$F1BS_240, cal$F1SS_240)

overlap <- function(a, b) {
  la <- log(a); lb <- log(b)
  br <- seq(min(c(la, lb)), max(c(la, lb)), length.out = 61)
  ha <- hist(la, breaks = br, plot = FALSE)$density; hb <- hist(lb, breaks = br, plot = FALSE)$density
  c(GMR_240_over_960 = exp(mean(la) - mean(lb)),
    overlap_coefficient = sum(pmin(ha, hb) * diff(br)),
    P_240_exceeds_960 = mean(outer(a, b, ">")),
    frac_240_below_960_p10 = mean(a < quantile(b, 0.10)),
    CV_pct_960 = 100 * sqrt(exp(var(lb)) - 1),
    CV_pct_240 = 100 * sqrt(exp(var(la)) - 1))
}
res <- overlap(a240$AUC_tau, a960$AUC_tau)

cat("\n=== rxode2 exposure overlap, scenario B (calibrated 240 mg) ===\n")
print(round(res, 3))
cat(sprintf("\nPython parity: GMR 0.819 (R %.3f), overlap 0.833 (R %.3f), CV 59.3%% (R %.1f%%)\n",
            res["GMR_240_over_960"], res["overlap_coefficient"], res["CV_pct_960"]))
cat(sprintf("Cmin GMR 240/960: %.3f   Cmax GMR 240/960: %.3f\n",
            exp(mean(log(a240$Cmin)) - mean(log(a960$Cmin))), exp(mean(log(a240$Cmax)) - mean(log(a960$Cmax)))))

write.csv(bind_rows(a960, a240), "results/population_exposures_R.csv", row.names = FALSE)
write_json(as.list(res), "results/exposure_overlap_R.json", auto_unbox = TRUE, pretty = TRUE)

## Figure
png("figures/fig2_exposure_overlap_R.png", width = 1300, height = 800, res = 200)
br <- 10^seq(log10(3), log10(200), length.out = 50)
hist(a960$AUC_tau / 1000, breaks = br, col = rgb(0.12, 0.31, 0.47, 0.55), border = NA,
     xlab = "Steady-state AUCtau (h*ug/mL), log scale", main = sprintf("rxode2: GMR 240/960 = %.2f, overlap %.2f", res[1], res[2]), log = "x")
hist(a240$AUC_tau / 1000, breaks = br, col = rgb(0.75, 0.22, 0.17, 0.55), border = NA, add = TRUE)
legend("topright", c("960 mg", "240 mg (calibrated F1)"), fill = c(rgb(0.12,0.31,0.47,0.55), rgb(0.75,0.22,0.17,0.55)), bty = "n")
dev.off()
cat("\nWritten: results/population_exposures_R.csv, exposure_overlap_R.json, figures/fig2_exposure_overlap_R.png\n")
