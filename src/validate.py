"""
Validation of the model implementation against published targets.

Targets (AAPS J 2025;27:26):
  Fig. 3: steady-state Cmax and AUCtau ratios for covariate contrasts vs the typical subject
  Model Based Simulations: Day 1 and Day 8 AUC0-24 and Cmax for 960 mg QD in patients

Pass criterion: within 5% of the published ratio. Absolute Day 1 / Day 8 values are
reported for transparency; they are population-level means over a covariate mix the
paper does not fully specify, so agreement is expected to be approximate.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from model import load_params, individual_params, simulate, interval_metrics, typical_subject

ROOT = Path(__file__).resolve().parents[1]
N_DOSES_SS = 30  # ~ 8.6 induction half-lives; steady state

# AAPS Fig. 3 published targets (Cmax,ss ratio, AUCtau,ss ratio)
TARGETS = {
    "Low albumin (30 g/L) vs normal (40)":  (dict(alb_gL=30.0), 1.075, 1.414),
    "ECOG 0 vs ECOG 1":                     (dict(ecog=0), 0.977, 0.868),
    "ECOG 2 vs ECOG 1":                     (dict(ecog=2), 1.017, 1.083),
    "Asian vs Caucasian":                   (dict(race="Asian"), 0.966, 0.817),
    "Black vs Caucasian":                   (dict(race="Black"), 0.988, 0.934),
    "Female vs Male":                       (dict(sex="F"), 1.242, 1.190),
    "Tumor <=70 mm vs >70 mm":              (dict(tumor="<=70"), 0.982, 0.892),
    "High-fat meal vs fasted":              (dict(highfat=True), 1.130, 1.327),
    "PPI vs no PPI":                        (dict(ppi=True), 0.828, 0.820),
    "Healthy (0 mm, ECOG 0) vs patients":   (dict(tumor="healthy", ecog=0), 0.898, 0.534),
}


def ss_metrics(cov_overrides: dict) -> dict:
    p = load_params()
    cov = typical_subject(960)
    cov.update(cov_overrides)
    ip = individual_params(p, cov)
    df = simulate(ip, 960, N_DOSES_SS)
    return interval_metrics(df, N_DOSES_SS)


def main():
    ref = ss_metrics({})
    rows = []
    for label, (ov, t_cmax, t_auc) in TARGETS.items():
        m = ss_metrics(ov)
        r_cmax = m["Cmax"] / ref["Cmax"]
        r_auc = m["AUC_tau"] / ref["AUC_tau"]
        rows.append(dict(contrast=label,
                         Cmax_ratio_model=round(r_cmax, 3), Cmax_ratio_published=t_cmax,
                         Cmax_pct_err=round((r_cmax / t_cmax - 1) * 100, 1),
                         AUC_ratio_model=round(r_auc, 3), AUC_ratio_published=t_auc,
                         AUC_pct_err=round((r_auc / t_auc - 1) * 100, 1)))
    out = pd.DataFrame(rows)
    out["pass_5pct"] = (out.Cmax_pct_err.abs() <= 5) & (out.AUC_pct_err.abs() <= 5)

    # Absolute-scale check, typical subject, 960 mg
    p = load_params()
    ip = individual_params(p, typical_subject(960))
    df = simulate(ip, 960, 8)
    d1, d8 = interval_metrics(df, 1), interval_metrics(df, 8)
    abs_rows = pd.DataFrame([
        dict(metric="AUC0-24 Day 1 (h*ng/mL)", model_typical=round(d1["AUC_tau"]), published_model_pred=53600, published_observed=57300),
        dict(metric="AUC0-24 Day 8 (h*ng/mL)", model_typical=round(d8["AUC_tau"]), published_model_pred=35300, published_observed=37200),
        dict(metric="Cmax Day 1 (ng/mL)",      model_typical=round(d1["Cmax"]),    published_model_pred=6600,  published_observed=7820),
        dict(metric="Cmax Day 8 (ng/mL)",      model_typical=round(d8["Cmax"]),    published_model_pred=4970,  published_observed=5970),
        dict(metric="Tmax Day 8 (h)",          model_typical=round(d8["Tmax"], 2), published_model_pred=None,  published_observed=1.1),
        dict(metric="Day8/Day1 AUC ratio",     model_typical=round(d8["AUC_tau"] / d1["AUC_tau"], 3), published_model_pred=round(35300 / 53600, 3), published_observed=round(37200 / 57300, 3)),
    ])

    (ROOT / "results").mkdir(exist_ok=True)
    out.to_csv(ROOT / "results" / "validation_covariate_ratios.csv", index=False)
    abs_rows.to_csv(ROOT / "results" / "validation_absolute_scale.csv", index=False)
    pd.set_option("display.width", 200)
    print(out.to_string(index=False))
    print()
    print(abs_rows.to_string(index=False))
    print(f"\nCovariate ratio contrasts passing 5% tolerance: {out.pass_5pct.sum()}/{len(out)}")


if __name__ == "__main__":
    main()
