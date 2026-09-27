## 02_simulate_and_recover.R
## Alternative full-scale estimation: 150 synthetic subjects, correlated omega,
## 300 burn-in and 300 estimation iterations. Not the source of the
## reported recovery table; see R/README.md for the executed reduced-scale run.
## Synthetic estimation demonstration; patient-level trial data are not public.
## RUN FROM THE REPOSITORY ROOT: source("R/02_simulate_and_recover.R")
## Runtime depends on the machine and package versions.

suppressPackageStartupMessages({ library(nlmixr2); library(rxode2); library(dplyr) })
if (!dir.exists("results")) stop("Run from the repository root.")
set.seed(20260906)

prm <- read.csv("data/poppk_parameters_final.csv", stringsAsFactors = FALSE)
P   <- setNames(prm$final_estimate, prm$parameter)

## ---- TRUE parameters (960 mg dose group, no covariates) ------------------
## F1 for 960 mg is fixed and known (F1SS_DG5 = 1 reference; F1BS_DG5 = 1.52).
## Bioavailability is not separately identifiable from oral data, so it is
## not estimated; this mirrors how the published model fixes F1SS_DG5.
truth <- c(ka = P[["Ka"]], v2 = P[["V2"]], clbs = P[["CLBS"]], clss = P[["CLSS"]],
           kind = P[["KIND_F"]], k24 = P[["K24"]], k42 = P[["K42"]],
           f1bs = P[["F1BS_DG5"]], f1ss = P[["F1SS_DG5"]])

## Published omega (log-scale variances and covariances) for (V2, CL, Ka)
omega_true <- lotri({
  eta.v2 + eta.cl + eta.ka ~ c(P[["OMEGA_1_1"]],
                               P[["OMEGA_2_1"]], P[["OMEGA_2_2"]],
                               P[["OMEGA_3_1"]], P[["OMEGA_3_2"]], P[["OMEGA_3_3"]])
})

gen <- rxode2({
  ka   <- TKA   * exp(eta.ka)
  v2   <- TV2   * exp(eta.v2)
  clbs <- TCLBS * exp(eta.cl)
  clss <- TCLSS * exp(eta.cl)
  cl   <- clss + (clbs - clss) * exp(-TKIND * t)
  f(tr1) <- TF1SS + (TF1BS - TF1SS) * exp(-TKIND * t)
  d/dt(tr1)   <- -ka * tr1
  d/dt(tr2)   <-  ka * tr1 - ka * tr2
  d/dt(tr3)   <-  ka * tr2 - ka * tr3
  d/dt(centr) <-  ka * tr3 - (cl / v2) * centr - TK24 * centr + TK42 * peri
  d/dt(peri)  <-  TK24 * centr - TK42 * peri
  cp <- centr / v2 * 1000        # dose in mg, v2 in L -> ng/mL
})
gen_par <- c(TKA = truth[["ka"]], TV2 = truth[["v2"]], TCLBS = truth[["clbs"]], TCLSS = truth[["clss"]],
             TKIND = truth[["kind"]], TK24 = truth[["k24"]], TK42 = truth[["k42"]],
             TF1BS = truth[["f1bs"]], TF1SS = truth[["f1ss"]])

## ---- Sampling design ----------------------------------------------------
## Mirrors the CodeBreaK 100 part B design criticized by Popat & Ratain:
## intensive sampling on Day 1 and Day 8 only. Set RICH <- TRUE to add
## Day 15 and Day 22 troughs for an exploratory design comparison.
RICH <- FALSE
n_subj <- 150
samp <- c(0.5, 1, 2, 4, 6, 8, 12, 24,                  # Day 1
          168.5, 169, 170, 172, 174, 176, 180, 192)    # Day 8
if (RICH) samp <- c(samp, 336, 504)                    # Day 15, Day 22 troughs
n_doses <- if (RICH) 22 else 8

ev <- et(amt = 960, ii = 24, addl = n_doses - 1, cmt = "tr1") %>% et(samp) %>% et(id = 1:n_subj)
sim <- rxSolve(gen, gen_par, ev, omega = omega_true, returnType = "data.frame")

## ---- Residual error: the published Dosne-type structure ------------------
## log(Yobs) = log(Ypred) + sqrt(theta_x^2 + theta_y^2 / Ypred^2) * eps,  eps ~ N(0,1)
## theta_x = ERR_EXP = 0.625 (log-scale), theta_y = ERR_ADD = 3.8 ng/mL.
## Simulating on the log scale (as published) avoids the negative
## concentrations a natural-scale 62.5% proportional error would produce.
tx <- P[["ERR_EXP"]]; ty <- P[["ERR_ADD"]]
obs <- sim %>%
  filter(cp > 0) %>%
  mutate(sd_log = sqrt(tx^2 + ty^2 / cp^2),
         DV = exp(log(cp) + sd_log * rnorm(n()))) %>%
  transmute(ID = id, TIME = time, DV, AMT = 0, EVID = 0, CMT = "cp")

dose_rows <- expand.grid(ID = 1:n_subj, TIME = 24 * (0:(n_doses - 1))) %>%
  transmute(ID, TIME, DV = NA_real_, AMT = 960, EVID = 1, CMT = "tr1")

dat <- bind_rows(dose_rows, obs) %>% arrange(ID, TIME, desc(EVID))
write.csv(dat, "results/synthetic_popPK_dataset.csv", row.names = FALSE)
cat(sprintf("Synthetic dataset: %d subjects, %d observations, design = %s\n",
            n_subj, nrow(obs), if (RICH) "Day 1/8/15/22" else "Day 1/8 only (as in the trial)"))

## ---- Estimation model: same structure, initial values deliberately off ----
est <- function() {
  ini({
    tka   <- log(4.0)      ; label("Ka (1/h); truth 7.87")
    tv2   <- log(150)      ; label("V2/F (L); truth 220")
    tclbs <- log(15)       ; label("CL/F baseline (L/h); truth 21.6")
    tclss <- log(30)       ; label("CL/F steady state (L/h); truth 41.3")
    tkind <- log(0.005)    ; label("Induction rate (1/h); truth 0.00845")
    tk24  <- log(0.02)     ; label("k24 (1/h); truth 0.0269")
    tk42  <- log(0.06)     ; label("k42 (1/h); truth 0.0421")
    eta.v2 + eta.cl + eta.ka ~ c(0.15,
                                 0.02, 0.15,
                                 0.00, 0.00, 0.15)
    lsd <- 0.5             ; label("log-scale residual SD; truth ~0.625")
  })
  model({
    ka   <- exp(tka + eta.ka)
    v2   <- exp(tv2 + eta.v2)
    clbs <- exp(tclbs + eta.cl)
    clss <- exp(tclss + eta.cl)
    kind <- exp(tkind)
    k24  <- exp(tk24)
    k42  <- exp(tk42)
    cl   <- clss + (clbs - clss) * exp(-kind * t)
    f(tr1) <- 1.00 + (1.52 - 1.00) * exp(-kind * t)   # 960 mg F1, fixed as published
    d/dt(tr1)   <- -ka * tr1
    d/dt(tr2)   <-  ka * tr1 - ka * tr2
    d/dt(tr3)   <-  ka * tr2 - ka * tr3
    d/dt(centr) <-  ka * tr3 - (cl / v2) * centr - k24 * centr + k42 * peri
    d/dt(peri)  <-  k24 * centr - k42 * peri
    cp <- centr / v2 * 1000
    cp ~ lnorm(lsd)
  })
}
## Note on the error model: the published structure adds a small additive
## term (3.8 ng/mL) inside the log-scale SD. At the concentrations sampled
## here it is negligible relative to 0.625, so a pure lognormal residual is
## used for estimation. State this in the write-up.

cat("\nRunning nlmixr2 SAEM ...\n")
fit <- nlmixr2(est, dat, est = "saem",
               control = saemControl(nBurn = 300, nEm = 300, print = 50, seed = 1))

## ---- Recovery table ------------------------------------------------------
th <- fit$theta
rec <- data.frame(
  parameter = c("Ka", "V2/F", "CLBS/F", "CLSS/F", "Kind", "k24", "k42"),
  truth     = c(truth[["ka"]], truth[["v2"]], truth[["clbs"]], truth[["clss"]], truth[["kind"]], truth[["k24"]], truth[["k42"]]),
  estimate  = exp(c(th[["tka"]], th[["tv2"]], th[["tclbs"]], th[["tclss"]], th[["tkind"]], th[["tk24"]], th[["tk42"]]))
) %>% mutate(pct_bias = round((estimate / truth - 1) * 100, 1))

om <- fit$omega
rec_omega <- data.frame(
  element = c("omega V2", "omega CL", "omega Ka", "cov V2-CL", "cov V2-Ka", "cov CL-Ka"),
  truth = c(P[["OMEGA_1_1"]], P[["OMEGA_2_2"]], P[["OMEGA_3_3"]], P[["OMEGA_2_1"]], P[["OMEGA_3_1"]], P[["OMEGA_3_2"]]),
  estimate = c(om["eta.v2","eta.v2"], om["eta.cl","eta.cl"], om["eta.ka","eta.ka"],
               om["eta.v2","eta.cl"], om["eta.v2","eta.ka"], om["eta.cl","eta.ka"])
)

cat("\n=== FIXED-EFFECT RECOVERY ===\n"); print(rec, row.names = FALSE)
cat("\n=== RANDOM-EFFECT RECOVERY (log-scale variance / covariance) ===\n"); print(rec_omega, row.names = FALSE)
cat("\n=== nlmixr2 fit summary (RSE and shrinkage if available) ===\n"); print(fit)

shr <- tryCatch(fit$shrink, error = function(e) NULL)
if (!is.null(shr)) { cat("\nEta shrinkage (%):\n"); print(shr) }
cat("\nA single recovery run does not isolate the effect of sampling design.\n")
cat("Compare repeated sparse and rich runs to assess bias and precision.\n")

write.csv(rec,       "results/parameter_recovery_nlmixr2.csv", row.names = FALSE)
write.csv(rec_omega, "results/omega_recovery_nlmixr2.csv",     row.names = FALSE)
saveRDS(fit, "results/nlmixr2_fit.rds")
cat("\nWritten: results/parameter_recovery_nlmixr2.csv, omega_recovery_nlmixr2.csv, nlmixr2_fit.rds\n")

## Standard diagnostics (uncomment after the fit succeeds):
## library(ggplot2)
## p1 <- plot(fit)                       # GOF panels
## vpc <- nlmixr2::vpcPlot(fit, n = 500)  # VPC (requires vpc package)
## ggsave("figures/nlmixr2_gof.png", p1[[1]], width = 8, height = 6, dpi = 150)
