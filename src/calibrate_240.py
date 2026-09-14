"""
The timing argument, quantified.

Hochmair 2024 Fig. 3 reports the 960/240 ratio of mean AUC0-24h:
    Day 1: 1.5
    Day 8: 1.3
Popat & Ratain note PK sampling stopped at Day 8 while steady state is not reached
until ~Day 22, and argue the gap would continue to narrow. They do not quantify it.

The published popPK pools 240 mg into dose group 1 with 120 and 180 mg (F1BS 4.95,
F1SS 4.58). Those values predict a Day 1 ratio of ~1.23 and a steady-state ratio of
~0.87, i.e. they overpredict 240 mg exposure relative to what was observed.

Approach:
  Treat F1BS_240 and F1SS_240 as unknown. Keep every other parameter (Kind, CLBS, CLSS,
  V2, K24, K42, Ka, and the 960 mg F1 terms) at the published values. Solve for the two
  240 mg F1 values that reproduce the observed Day 1 and Day 8 ratios in the typical
  patient. Then simulate forward to Day 30 to project the model-based steady-state ratio.

  Uncertainty: the observed ratios are ratios of arithmetic means over ~104 patients per
  arm with large SDs (Fig. 3 error bars). Fig. 3 does not report the SDs numerically.
  We independently perturb each ratio lognormally with log-scale SD 0.15
  (approximately 15% relative uncertainty). This is a sensitivity assumption, not
  an estimated standard error or a confidence interval derived from trial data.

Two exactly-identified unknowns from two observations: no goodness-of-fit is possible.
The value of the exercise is the *projection*, and the internal consistency check that
the implied 240 mg autoinduction loss falls where the published dose-dependence pattern
says it should (180 mg: 7.5%, 360 mg: 16.6%).
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from model import load_params, individual_params, simulate, interval_metrics, typical_subject

ROOT = Path(__file__).resolve().parents[1]
OBS_D1, OBS_D8 = 1.5, 1.3     # Hochmair Fig. 3, 960/240 AUC0-24 ratio
REL_UNC = 0.15                # assumed relative uncertainty on each observed ratio
N_DOSES = 30


def auc_by_day(ip, dose):
    df = simulate(ip, dose, N_DOSES)
    return np.array([interval_metrics(df, k)["AUC_tau"] for k in range(1, N_DOSES + 1)])


def f1_t(f1bs, f1ss, kind, t):
    return f1ss + (f1bs - f1ss) * np.exp(-kind * t)


def ratio_curve_analytic(p, f1bs_240, f1ss_240, n_days=N_DOSES):
    """
    960/240 ratio of AUC0-24 on each day. CL, V2, Ka, K24, K42 are identical between
    arms for the same subject and cancel in the ratio; the residual carryover
    difference is second order (verified against the ODE below, <1%).
    """
    t = 24.0 * np.arange(n_days)
    num = 960 * f1_t(p["F1BS_DG5"], p["F1SS_DG5"], p["KIND_F"], t)
    den = 240 * f1_t(f1bs_240, f1ss_240, p["KIND_F"], t)
    return num / den


def ratio_curve_ode(p, f1bs_240, f1ss_240):
    ip960 = individual_params(p, typical_subject(960))
    ip240 = individual_params(p, typical_subject(240), f1bs_override=f1bs_240, f1ss_override=f1ss_240)
    return auc_by_day(ip960, 960) / auc_by_day(ip240, 240)


def fit(p, obs_d1, obs_d8):
    def resid(x):
        r = ratio_curve_analytic(p, x[0], x[1])
        return [r[0] - obs_d1, r[7] - obs_d8]
    sol = least_squares(resid, x0=[4.0, 3.3], bounds=([0.5, 0.5], [10, 10]))
    return sol.x


def main():
    p = load_params()

    # Published DG1 baseline
    r_pub = ratio_curve_analytic(p, p["F1BS_DG1"], p["F1SS_DG1"])

    # Point calibration
    f1bs, f1ss = fit(p, OBS_D1, OBS_D8)
    r_cal = ratio_curve_analytic(p, f1bs, f1ss)
    r_cal_ode = ratio_curve_ode(p, f1bs, f1ss)
    ode_check = dict(day1=[float(r_cal[0]), float(r_cal_ode[0])], day8=[float(r_cal[7]), float(r_cal_ode[7])], day30=[float(r_cal[29]), float(r_cal_ode[29])])

    # Uncertainty: sample observed ratios, refit
    rng = np.random.default_rng(1)
    draws = []
    for _ in range(2000):
        d1 = OBS_D1 * np.exp(rng.normal(0, REL_UNC))
        d8 = OBS_D8 * np.exp(rng.normal(0, REL_UNC))
        try:
            b, s = fit(p, d1, d8)
            draws.append(ratio_curve_analytic(p, b, s))
        except Exception:
            pass
    draws = np.array(draws)
    lo, hi = np.percentile(draws, [2.5, 97.5], axis=0)

    days = np.arange(1, N_DOSES + 1)
    out = pd.DataFrame(dict(day=days, ratio_published_DG1=r_pub, ratio_calibrated=r_cal,
                            ratio_calibrated_lo95=lo, ratio_calibrated_hi95=hi))
    out.to_csv(ROOT / "results" / "exposure_ratio_timecourse.csv", index=False)

    induction_loss_240 = 1 - f1ss / f1bs
    induction_loss_960 = 1 - p["F1SS_DG5"] / p["F1BS_DG5"]
    summary = dict(
        F1BS_240=float(f1bs), F1SS_240=float(f1ss),
        implied_240mg_autoinduction_F_loss_pct=float(100 * induction_loss_240),
        published_960mg_autoinduction_F_loss_pct=float(100 * induction_loss_960),
        published_pattern_pct={"180mg": 7.5, "360mg": 16.6, "720mg": 31.3, "960mg": 34.2},
        ratio_day1=float(r_cal[0]), ratio_day8=float(r_cal[7]),
        ratio_day15=float(r_cal[14]), ratio_day22=float(r_cal[21]), ratio_day30=float(r_cal[29]),
        ratio_day30_ci95=[float(lo[29]), float(hi[29])],
        implied_240_over_960_at_SS=float(1 / r_cal[29]),
        implied_240_over_960_at_SS_ci95=[float(1 / hi[29]), float(1 / lo[29])],
        published_DG1_ratio_day1=float(r_pub[0]), published_DG1_ratio_day8=float(r_pub[7]),
        published_DG1_ratio_day30=float(r_pub[29]),
        analytic_vs_ode_check=ode_check,
        assumptions="Kind, CL, V2, K24, K42, Ka and 960 mg F1 at published values; independent lognormal perturbations of observed ratios with assumed log-scale SD 0.15",
        interval_interpretation="95% sensitivity intervals (2.5th/97.5th percentiles); legacy ci95 keys are retained for compatibility, not conventional confidence intervals",
    )
    (ROOT / "results" / "calibrated_240mg_F1.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
