"""
Component C: Steady-state exposure overlap between 240 mg QD and 960 mg QD.

Virtual population covariates are drawn from the CodeBreaK 100 part B baseline table
(Hochmair 2024, Table 1) where available, and from the NSCLC column of the popPK
analysis dataset (AAPS Table I) where part B does not report the distribution.

  Source: Hochmair Table 1 (pooled arms, n=209)
    ECOG 0/1/2      : 35.4% / 56.9% / 7.7%
    Male            : 54.5%
    Asian           : 17.2%   (no Black / NHPI / AIAN patients enrolled)
  Source: AAPS Table I, NSCLC column (n=258)
    Tumor >70 mm    : 47.5%  (>70 vs <=70; healthy category not applicable)
    Albumin g/L     : median 37.9, range 21-51  -> Normal(37.9, 5.0) truncated [21, 51]
                      SD is NOT published; 5.0 g/L is an assumption chosen so that
                      ~20% of subjects fall below 34 g/L (hypoalbuminemia threshold).
    PPI use         : 37% (AAPS text: 160/431 patients)

Variability and uncertainty are represented separately:
  1. Between-subject variability (IIV): the published 3x3 omega matrix on V2, CL, Ka.
  2. Parameter uncertainty on the dose-group F1 terms: sampled lognormally from the
     published 95% CI (AAPS Table II). This is the dominant source of uncertainty in
     published-DG1 scenario. For calibrated scenario B, the interval instead propagates
     assumed uncertainty in the Day 1/Day 8 calibration targets (log-scale SD 0.15).

Two scenarios for 240 mg bioavailability:
  A. "Published DG1": use F1BS_DG1 / F1SS_DG1 as published (pooled 120/180/240 mg).
  B. "Calibrated 240": use F1 values for 240 mg calibrated to the observed Day 1 and
     Day 8 exposure ratios in Hochmair Fig. 3 (see calibrate_240.py). This is the
     scenario used for the headline result, because scenario A is known to misfit
     (Hochmair Discussion; ANALYSIS_NOTES.md section 2).
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from model import load_params, individual_params, simulate, interval_metrics, omega_matrix

ROOT = Path(__file__).resolve().parents[1]
RNG = np.random.default_rng(20260901)
N_PER_ARM = 1000
N_DOSES = 30       # steady state
N_PARAM_DRAWS = 200  # for parameter-uncertainty layer (cheap analytic AUC)


def draw_covariates(n: int) -> list[dict]:
    ecog = RNG.choice([0, 1, 2], size=n, p=[0.354, 0.569, 0.077])
    sex = RNG.choice(["M", "F"], size=n, p=[0.545, 0.455])
    race = RNG.choice(["Caucasian", "Asian"], size=n, p=[0.828, 0.172])
    tumor = RNG.choice([">70", "<=70"], size=n, p=[0.475, 0.525])
    alb = np.clip(RNG.normal(37.9, 5.0, size=n), 21, 51)
    ppi = RNG.random(n) < 0.37
    return [dict(ecog=int(ecog[i]), sex=sex[i], race=race[i], tumor=tumor[i],
                 alb_gL=float(alb[i]), ppi=bool(ppi[i]), highfat=False) for i in range(n)]


def ss_auc_analytic(ip: dict, dose_mg: float) -> float:
    """AUC0-tau at steady state for a linear system: F*Dose/CL. Units h*ng/mL."""
    return dose_mg * ip["F1SS"] / ip["CLSS"] * 1000.0


def analytic_arm(p, covs, etas, dose_mg, f1bs=None, f1ss=None) -> pd.DataFrame:
    rows = []
    for i, (cov, eta) in enumerate(zip(covs, etas)):
        c = dict(cov, dose_mg=dose_mg)
        ip = individual_params(p, c, eta, f1bs_override=f1bs, f1ss_override=f1ss)
        rows.append(dict(id=i, dose_mg=dose_mg, AUC_tau=ss_auc_analytic(ip, dose_mg),
                         CLSS_i=ip["CLSS"], F1SS_i=ip["F1SS"], **cov))
    return pd.DataFrame(rows)


def ode_subsample(p, covs, etas, dose_mg, n=150, f1bs=None, f1ss=None) -> pd.DataFrame:
    rows = []
    for i in range(n):
        c = dict(covs[i], dose_mg=dose_mg)
        ip = individual_params(p, c, etas[i], f1bs_override=f1bs, f1ss_override=f1ss)
        df = simulate(ip, dose_mg, N_DOSES, dt=0.5)
        m = interval_metrics(df, N_DOSES)
        rows.append(dict(id=i, dose_mg=dose_mg, **m))
    return pd.DataFrame(rows)


def overlap_stats(a: np.ndarray, b: np.ndarray) -> dict:
    """a = 240 mg AUCs, b = 960 mg AUCs."""
    gmr = np.exp(np.log(a).mean() - np.log(b).mean())
    med_ratio = np.median(a) / np.median(b)
    lo, hi = min(a.min(), b.min()), max(a.max(), b.max())
    bins = np.linspace(np.log(lo), np.log(hi), 60)
    ha, _ = np.histogram(np.log(a), bins=bins, density=True)
    hb, _ = np.histogram(np.log(b), bins=bins, density=True)
    ovl = np.sum(np.minimum(ha, hb) * np.diff(bins))
    p_sup = np.mean(a[:, None] > b[None, :])
    p10_960 = np.percentile(b, 10)
    tail = np.mean(a < p10_960)
    return dict(GMR_240_over_960=gmr, median_ratio=med_ratio, overlap_coefficient=ovl,
                P_240_exceeds_960=p_sup, frac_240_below_960_p10=tail,
                p10_960=p10_960, p10_240=np.percentile(a, 10),
                p50_240=np.median(a), p50_960=np.median(b),
                p90_240=np.percentile(a, 90), p90_960=np.percentile(b, 90),
                CV_pct_240=100 * np.sqrt(np.exp(np.var(np.log(a))) - 1),
                CV_pct_960=100 * np.sqrt(np.exp(np.var(np.log(b))) - 1))


def main():
    p = load_params()
    covs = draw_covariates(N_PER_ARM)
    etas = RNG.multivariate_normal(np.zeros(3), omega_matrix(p), size=N_PER_ARM)
    calib = json.loads((ROOT / "results" / "calibrated_240mg_F1.json").read_text())

    arm960 = analytic_arm(p, covs, etas, 960).assign(scenario="960")
    armA = analytic_arm(p, covs, etas, 240).assign(scenario="A_published_DG1")
    armB = analytic_arm(p, covs, etas, 240, f1bs=calib["F1BS_240"], f1ss=calib["F1SS_240"]).assign(scenario="B_calibrated_240")

    results = {
        "scenario_A_published_DG1": overlap_stats(armA.AUC_tau.values, arm960.AUC_tau.values),
        "scenario_B_calibrated_240": overlap_stats(armB.AUC_tau.values, arm960.AUC_tau.values),
    }

    # Parameter uncertainty on the population GMR (scenario A): depends only on F1SS_DG1
    lo, hi = p["_ci"]["F1SS_DG1"]
    sd_log = (np.log(float(hi)) - np.log(float(lo))) / (2 * 1.96)
    f1_draws = np.exp(RNG.normal(np.log(p["F1SS_DG1"]), sd_log, size=5000))
    g = 240 * f1_draws / 960.0
    results["scenario_A_GMR_param_uncertainty"] = dict(median=float(np.median(g)),
        ci95=[float(np.percentile(g, 2.5)), float(np.percentile(g, 97.5))])
    results["scenario_B_GMR_param_uncertainty"] = dict(
        median=calib["implied_240_over_960_at_SS"], ci95=calib["implied_240_over_960_at_SS_ci95"],
        interval_interpretation="95% sensitivity interval under independent calibration-target lognormal perturbations with assumed log-scale SD 0.15; legacy ci95 key retained")

    # ODE subsample for Cmax / Cmin (scenario B and 960)
    print("ODE subsample for Cmax/Cmin ...")
    o960 = ode_subsample(p, covs, etas, 960)
    oB = ode_subsample(p, covs, etas, 240, f1bs=calib["F1BS_240"], f1ss=calib["F1SS_240"])
    results["cmax_cmin_subsample_n150"] = dict(
        Cmax_GMR_240_over_960=float(np.exp(np.log(oB.Cmax).mean() - np.log(o960.Cmax).mean())),
        Cmin_GMR_240_over_960=float(np.exp(np.log(oB.Cmin).mean() - np.log(o960.Cmin).mean())),
        Cmin_median_240=float(oB.Cmin.median()), Cmin_median_960=float(o960.Cmin.median()),
        AUC_GMR_check=float(np.exp(np.log(oB.AUC_tau).mean() - np.log(o960.AUC_tau).mean())))

    pd.concat([arm960, armA, armB]).to_csv(ROOT / "results" / "population_exposures.csv", index=False)
    (ROOT / "results" / "exposure_overlap.json").write_text(json.dumps(results, indent=2, default=float))
    print(json.dumps(results, indent=2, default=float))


if __name__ == "__main__":
    main()
