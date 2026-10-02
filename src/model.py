"""
Sotorasib population PK model.

Structure (Nagase et al., AAPS J 2025;27:26, Fig. 1):
  Dose --F1(t)--> TR1 --Ka--> TR2 --Ka--> TR3 --Ka--> Central (V2) <--K24/K42--> Peripheral
                                                        |
                                                      CL(t)/V2

Time-dependent autoinduction (same Kind modulates both F1 and CL):
  F1(t) = F1SS + (F1BS - F1SS) * exp(-Kind * t)
  CL(t) = CLSS + (CLBS - CLSS) * exp(-Kind * t)
where t is time after first dose (TAFD), in hours. This is an implementation of
selected published point estimates, not an independently fitted clinical model.

Relative bioavailability is parameterized by dose group with 960 mg at steady state
fixed to 1. V2 and CL are therefore *apparent* (V2/F, CL/F) relative to that reference.

Covariates (AAPS Table II; operational assumptions in data/PARAMETER_RESOLUTION.md):
  CL:  ECOG, sex, tumor size category, race, albumin (proportional per g/L, ref 40 g/L)
  V2:  sex
  F1:  PPI, high-fat meal
  Ka:  high-fat meal
Covariate effects are applied multiplicatively as (1 + theta). This form reproduces the
published forest-plot ratios (see validate.py).

IIV: lognormal on V2, CL, Ka with the published 3x3 covariance matrix (log scale).
No IIV on F1 (none was estimated).

All numeric values are loaded from data/poppk_parameters_final.csv.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
PARAM_CSV = ROOT / "data" / "poppk_parameters_final.csv"


def load_params() -> dict:
    df = pd.read_csv(PARAM_CSV)
    p = {r.parameter: float(r.final_estimate) for r in df.itertuples()}
    # Published confidence limits retained with the source parameters.
    ci = {r.parameter: (r.ci_lower, r.ci_upper) for r in df.itertuples()}
    p["_ci"] = ci
    return p


# Dose-group assignment (AAPS Table II row labels)
def dose_group(dose_mg: float) -> str:
    if dose_mg in (120, 180, 240):
        return "DG1"
    if dose_mg == 360:
        return "DG2"
    if dose_mg == 480:
        return "DG3"
    if dose_mg in (600, 720, 840):
        return "DG4"
    if dose_mg == 960:
        return "DG5"
    raise ValueError(f"Dose {dose_mg} mg is not in any published dose group")


def omega_matrix(p: dict) -> np.ndarray:
    """3x3 log-scale variance-covariance for (V2, CL, Ka), AAPS Table II."""
    return np.array([
        [p["OMEGA_1_1"], p["OMEGA_2_1"], p["OMEGA_3_1"]],
        [p["OMEGA_2_1"], p["OMEGA_2_2"], p["OMEGA_3_2"]],
        [p["OMEGA_3_1"], p["OMEGA_3_2"], p["OMEGA_3_3"]],
    ])


def individual_params(p: dict, cov: dict, eta: np.ndarray | None = None,
                      f1bs_override: float | None = None,
                      f1ss_override: float | None = None) -> dict:
    """
    Build the parameter set for one subject.

    cov keys (defaults = AAPS 'typical subject'):
      ecog: 0,1,2           (ref 1)
      sex: 'M'/'F'          (ref M)
      tumor: '>70', '<=70', 'healthy'   (ref '>70')
      race: 'Caucasian','Asian','Black','Other','NHPI','AIAN'  (ref Caucasian)
      alb_gL: albumin g/L   (ref 40)
      ppi: bool             (ref False)
      highfat: bool         (ref False)
      dose_mg: dose level, used only to pick the F1 dose group
    eta: length-3 array of log-scale random effects on (V2, CL, Ka); None = typical
    f1bs_override / f1ss_override: replace the dose-group F1 values (used for the
  exploratory scenarios only; no calibrated-dose analysis is in the active pipeline)
    """
    ecog = cov.get("ecog", 1)
    sex = cov.get("sex", "M")
    tumor = cov.get("tumor", ">70")
    race = cov.get("race", "Caucasian")
    alb = cov.get("alb_gL", 40.0)
    ppi = cov.get("ppi", False)
    highfat = cov.get("highfat", False)
    dose = cov.get("dose_mg", 960)

    # ---- CL covariate multiplier ----
    cl_mult = 1.0
    if ecog == 0:
        cl_mult *= 1 + p["CLECOGBL1"]
    elif ecog == 2:
        cl_mult *= 1 + p["CLECOGBL2"]
    if sex == "F":
        cl_mult *= 1 + p["CLSEX1"]
    if tumor == "<=70":
        cl_mult *= 1 + p["CLTUM_BS_CAT1"]
    elif tumor == "healthy":
        cl_mult *= 1 + p["CLTUM_BS_CAT2"]
    race_map = {"Asian": "CLRACE1", "Black": "CLRACE2", "Other": "CLRACE3",
                "NHPI": "CLRACE4", "AIAN": "CLRACE5"}
    if race in race_map:
        cl_mult *= 1 + p[race_map[race]]
    # Albumin: operational proportional-linear g/L assumption; source ambiguity documented
    cl_mult *= 1 + p["CLALB1"] * (alb - 40.0)

    # ---- V2 covariate multiplier ----
    v2_mult = 1.0
    if sex == "F":
        v2_mult *= 1 + p["S2SEX1"]

    # ---- F1 covariate multiplier ----
    f1_mult = 1.0
    if ppi:
        f1_mult *= 1 + p["F1PPI1"]
    if highfat:
        f1_mult *= 1 + p["F1HIGHFAT1"]

    # ---- Ka covariate multiplier ----
    ka_mult = 1.0
    if highfat:
        ka_mult *= 1 + p["KAHIGHFAT1"]

    # ---- IIV ----
    if eta is None:
        eta = np.zeros(3)
    v2 = p["V2"] * v2_mult * np.exp(eta[0])
    clbs = p["CLBS"] * cl_mult * np.exp(eta[1])
    clss = p["CLSS"] * cl_mult * np.exp(eta[1])
    ka = p["Ka"] * ka_mult * np.exp(eta[2])

    dg = dose_group(dose)
    f1bs = p[f"F1BS_{dg}"] if f1bs_override is None else f1bs_override
    f1ss = p[f"F1SS_{dg}"] if f1ss_override is None else f1ss_override

    return dict(V2=v2, CLBS=clbs, CLSS=clss, Ka=ka, K24=p["K24"], K42=p["K42"],
                Kind=p["KIND_F"], F1BS=f1bs * f1_mult, F1SS=f1ss * f1_mult)


def _rhs(t, y, ip):
    tr1, tr2, tr3, a2, a4 = y
    ka = ip["Ka"]
    cl = ip["CLSS"] + (ip["CLBS"] - ip["CLSS"]) * np.exp(-ip["Kind"] * t)
    return [
        -ka * tr1,
        ka * tr1 - ka * tr2,
        ka * tr2 - ka * tr3,
        ka * tr3 - (cl / ip["V2"]) * a2 - ip["K24"] * a2 + ip["K42"] * a4,
        ip["K24"] * a2 - ip["K42"] * a4,
    ]


def simulate(ip: dict, dose_mg: float, n_doses: int, tau: float = 24.0,
             dt: float = 0.25, rtol: float = 1e-8, atol: float = 1e-6) -> pd.DataFrame:
    """
    Simulate once-daily dosing. Returns a DataFrame with time (h) and concentration
    (ng/mL) in the central compartment. Each interval includes its exact endpoints.
    A shared boundary can occur twice, immediately before/after an oral dose.
    Drug is added to TR1, so central concentration is continuous at that event.

    Dose amount in mg -> ng: x 1e6. V2 in L, so conc = A2(ng) / V2(L) / 1000 = ng/mL.
    """
    if not isinstance(n_doses, (int, np.integer)) or n_doses < 1:
        raise ValueError("n_doses must be a positive integer")
    if not all(np.isfinite(v) and v > 0 for v in [dose_mg, tau, dt, rtol, atol]):
        raise ValueError("dose, tau, dt and solver tolerances must be positive and finite")
    y = np.zeros(5)
    out_t, out_c = [], []
    for k in range(n_doses):
        t0 = k * tau
        f1 = ip["F1SS"] + (ip["F1BS"] - ip["F1SS"]) * np.exp(-ip["Kind"] * t0)
        y[0] += dose_mg * 1e6 * f1  # ng into TR1
        # Construct relative times; append the exact endpoint to prevent floating
        # arange overshoot, including when dt does not divide tau (e.g. 0.7 h).
        offsets = np.arange(0.0, tau, dt)
        offsets = offsets[offsets < tau]
        t_eval = np.unique(np.minimum(t0 + np.r_[offsets, tau], t0 + tau))
        sol = solve_ivp(_rhs, (t0, t0 + tau), y, args=(ip,), t_eval=t_eval,
                        method="LSODA", rtol=rtol, atol=atol)
        if not sol.success:
            raise RuntimeError(sol.message)
        y = sol.y[:, -1].copy()
        out_t.append(sol.t)
        out_c.append(sol.y[3] / ip["V2"] / 1000.0)  # ng/mL
    t = np.concatenate(out_t)
    c = np.concatenate(out_c)
    df = pd.DataFrame({"time_h": t, "conc_ng_mL": c})
    return df


def interval_metrics(df: pd.DataFrame, dose_no: int, tau: float = 24.0) -> dict:
    """Closed-interval trapezoidal AUC, sampled Cmax/Tmax, and terminal C_tau.

    If a requested boundary lies between samples, its concentration is linearly
    interpolated. Extrapolation is forbidden. Consistent duplicate oral-dose
    boundary samples are collapsed. Unequal duplicate values are rejected: this
    helper does not define a pre/post convention for discontinuous central doses.
    C_tau is the final concentration, not necessarily the minimum in the interval.
    """
    if not isinstance(dose_no, (int, np.integer)) or dose_no < 1 or not np.isfinite(tau) or tau <= 0:
        raise ValueError("dose_no and tau must be positive")
    d = df[["time_h", "conc_ng_mL"]].copy().sort_values("time_h")
    if len(d) < 2 or not np.isfinite(d.to_numpy()).all():
        raise ValueError("Need at least two finite time/concentration samples")
    grouped = d.groupby("time_h")["conc_ng_mL"]
    spreads = grouped.max() - grouped.min()
    scales = np.maximum(1.0, grouped.max().abs())
    if np.any(spreads > 1e-10 * scales):
        raise ValueError("Unequal duplicate concentrations require an explicit event convention")
    d = d.drop_duplicates("time_h")
    times, concentrations = d.time_h.to_numpy(), d.conc_ng_mL.to_numpy()
    start, end = (dose_no - 1) * tau, dose_no * tau
    if start < times[0] or end > times[-1]:
        raise ValueError("Requested interval is outside the sampled time range; no extrapolation")
    interior = (times > start) & (times < end)
    t = np.r_[start, times[interior], end]
    c = np.r_[np.interp(start, times, concentrations), concentrations[interior],
              np.interp(end, times, concentrations)]
    t = t - start
    auc = np.trapezoid(c, t)
    return dict(AUC_tau=float(auc), Cmax=float(c.max()), Tmax=float(t[c.argmax()]), C_tau=float(c[-1]))


def typical_subject(dose_mg: float = 960) -> dict:
    """AAPS Methods definition of the reference subject."""
    return dict(ecog=1, sex="M", tumor=">70", race="Caucasian", alb_gL=40.0,
                ppi=False, highfat=False, dose_mg=dose_mg)
