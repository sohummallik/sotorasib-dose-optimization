"""Figures for the write-up. All inputs are results/*.csv and *.json produced by the pipeline."""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from model import load_params, individual_params, simulate, typical_subject

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
C960, C240 = "#1f4e79", "#c0392b"


def fig1_dose_exposure():
    """FDA p.94 observed steady-state AUC vs dose, with the model's typical-subject prediction."""
    obs = pd.DataFrame(dict(dose=[180, 360, 720, 960], auc=[31.7, 38.9, 42.1, 32.4], n=[6, 24, 11, 24]))
    p = load_params()
    pred = []
    for d in [180, 360, 480, 720, 960]:
        ip = individual_params(p, typical_subject(d))
        pred.append(d * ip["F1SS"] / ip["CLSS"])  # mg*h/L = ug*h/mL
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot([180, 360, 480, 720, 960], pred, "o-", color="grey", label="popPK model, typical subject (steady state)")
    ax.errorbar(obs.dose, obs.auc, fmt="s", color=C960, ms=8, label="Observed Day 8 geometric mean (FDA MDR p.94)")
    for _, r in obs.iterrows():
        ax.annotate(f"n={r.n}", (r.dose, r.auc), textcoords="offset points", xytext=(0, 8), ha="center", fontsize=8)
    ax.plot([180, 960], [31.7, 31.7 * 960 / 180], "--", color="lightgrey", label="Dose-proportional reference")
    ax.set_ylim(0, 60); ax.set_xlabel("Sotorasib dose (mg once daily)"); ax.set_ylabel("AUC0-24h (h*ug/mL)")
    ax.set_title("A 5.3-fold dose range produces a 1.02-fold exposure range")
    ax.legend(fontsize=8, loc="upper left"); fig.tight_layout()
    fig.savefig(FIG / "fig1_dose_vs_exposure.png", dpi=200); plt.close(fig)


def fig2_exposure_overlap():
    pop = pd.read_csv(ROOT / "results" / "population_exposures.csv")
    ov = json.loads((ROOT / "results" / "exposure_overlap.json").read_text())["scenario_B_calibrated_240"]
    a = pop[pop.scenario == "B_calibrated_240"].AUC_tau / 1000
    b = pop[pop.scenario == "960"].AUC_tau / 1000
    fig, ax = plt.subplots(figsize=(6.5, 4))
    bins = np.logspace(np.log10(3), np.log10(200), 50)
    ax.hist(b, bins=bins, alpha=0.55, color=C960, label=f"960 mg  (median {ov['p50_960']/1000:.1f})")
    ax.hist(a, bins=bins, alpha=0.55, color=C240, label=f"240 mg  (median {ov['p50_240']/1000:.1f})")
    ax.set_xscale("log"); ax.set_xlabel("Steady-state AUCtau (h*ug/mL), 1000 virtual patients per arm")
    ax.set_ylabel("Patients")
    ax.set_title(f"GMR 240/960 = {ov['GMR_240_over_960']:.2f}; overlap coefficient {ov['overlap_coefficient']:.2f}; BSV CV {ov['CV_pct_960']:.0f}%")
    ax.axvline(ov["p10_960"] / 1000, color=C960, ls=":", lw=1)
    ax.text(ov["p10_960"] / 1000, ax.get_ylim()[1] * 0.9, " 960 mg 10th pct", fontsize=8, color=C960)
    ax.legend(); fig.tight_layout()
    fig.savefig(FIG / "fig2_exposure_overlap.png", dpi=200); plt.close(fig)


def fig3_ratio_timecourse():
    d = pd.read_csv(ROOT / "results" / "exposure_ratio_timecourse.csv")
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.fill_between(d.day, d.ratio_calibrated_lo95, d.ratio_calibrated_hi95, color=C240, alpha=0.15, label="95% interval (assumed 15% error on observed ratios)")
    ax.plot(d.day, d.ratio_calibrated, color=C240, lw=2, label="Calibrated to observed Day 1 and Day 8")
    ax.plot(d.day, d.ratio_published_DG1, color="grey", ls="--", label="Published pooled dose-group 1 parameters")
    ax.plot([1, 8], [1.5, 1.3], "ko", ms=8, label="Observed (Hochmair 2024, Fig. 3)")
    ax.axhline(1.0, color="lightgrey", lw=1)
    ax.axvline(8, color="lightgrey", ls=":", lw=1); ax.text(8.3, 1.9, "last PK sample\nin the trial", fontsize=8, color="grey")
    ax.set_xlabel("Day of once-daily dosing"); ax.set_ylabel("AUC0-24h ratio, 960 mg / 240 mg")
    ax.set_ylim(0.6, 2.1); ax.set_title("Autoinduction narrows the exposure gap; most of it by Day 8")
    ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(FIG / "fig3_exposure_ratio_timecourse.png", dpi=200); plt.close(fig)


def fig4_er_efficacy():
    er = json.loads((ROOT / "results" / "decision_analysis.json").read_text())["er_efficacy"]
    q = pd.DataFrame(dict(auc=[13000, 19000, 27000, 46000], r=[25, 26, 21, 12], n=[57] * 4))
    q["orr"] = q.r / q.n
    from scipy.stats import beta as B
    lo = B.ppf(0.025, q.r, q.n - q.r + 1); hi = B.ppf(0.975, q.r + 1, q.n - q.r)
    x = np.logspace(np.log10(5000), np.log10(80000), 200)
    y = 1 / (1 + np.exp(-(er["intercept"] + er["logistic_slope_per_log_auc"] * np.log(x))))
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.errorbar(q.auc, q.orr, yerr=[q.orr - lo, hi - q.orr], fmt="s", color=C960, ms=8, capsize=4, label="FDA Fig. 26 quartiles (Clopper-Pearson 95% CI)")
    ax.plot(x, y, color=C240, lw=2, label=f"Continuous logistic fit: OR {er['OR_per_doubling_of_AUC']:.2f} per doubling (95% CI {er['OR_per_doubling_95CI'][0]:.2f} to {er['OR_per_doubling_95CI'][1]:.2f})")
    ax.set_xscale("log"); ax.set_ylim(0, 0.8)
    ax.set_xlabel("Model-predicted AUCtau,ss (h*ng/mL)"); ax.set_ylabel("Objective response rate")
    ax.set_title("Observed exposure-response is inverted, consistent with confounding by disease burden")    
    ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(FIG / "fig4_exposure_response_efficacy.png", dpi=200); plt.close(fig)


def fig5_tipping_point():
    d = pd.read_csv(ROOT / "results" / "orr_diff_vs_er_slope.csv")
    tp = json.loads((ROOT / "results" / "decision_analysis.json").read_text())["tipping_point"]
    s = pd.read_csv(ROOT / "results" / "utility_surface_240_minus_960.csv", index_col=0)
    betas = np.array(s.index, dtype=float); ws = np.array(s.columns, dtype=float)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
    ax1.plot(np.exp(d.beta * np.log(2)), 100 * d.orr_diff_960_minus_240, color=C960, lw=2)
    ax1.axhline(7.9, color=C240, ls="--", label="Observed ORR gap, 7.9 points")
    ax1.axvline(tp["OR_per_doubling_required"], color=C240, ls=":", label=f"Required OR = {tp['OR_per_doubling_required']:.1f} per doubling")
    ax1.axvline(0.55, color="grey", ls=":", label="Binned FDA exposure-response fit: OR = 0.55 per AUC doubling")
    ax1.set_xscale("log"); ax1.set_xlim(0.4, 8); ax1.set_xlabel("Assumed true exposure-efficacy OR per doubling of AUC")
    ax1.set_ylabel("Predicted ORR difference, 960 minus 240 (points)")
    ax1.set_title("How steep would E-R need to be to explain the gap?"); ax1.legend(fontsize=8)
    im = ax2.contourf(ws, np.exp(betas * np.log(2)), s.values, levels=[-0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2], cmap="RdBu")
    ax2.contour(ws, np.exp(betas * np.log(2)), s.values, levels=[0], colors="k")
    ax2.set_yscale("log"); ax2.set_xlabel("Toxicity weight w  (utility = P(response) - w * P(Grade>=3 TEAE))")
    ax2.set_ylabel("Assumed true E-R OR per doubling"); ax2.set_title("Where 240 mg wins (blue) vs 960 mg (red)")
    fig.colorbar(im, ax=ax2, label="Utility(240) - Utility(960)")
    fig.tight_layout(); fig.savefig(FIG / "fig5_tipping_point.png", dpi=200); plt.close(fig)


def fig6_conc_time():
    p = load_params()
    cal = json.loads((ROOT / "results" / "calibrated_240mg_F1.json").read_text())
    ip960 = individual_params(p, typical_subject(960))
    ip240 = individual_params(p, typical_subject(240), f1bs_override=cal["F1BS_240"], f1ss_override=cal["F1SS_240"])
    d960 = simulate(ip960, 960, 22); d240 = simulate(ip240, 240, 22)
    fig, ax = plt.subplots(figsize=(8, 3.6))
    ax.plot(d960.time_h / 24, d960.conc_ng_mL, color=C960, lw=1.2, label="960 mg QD")
    ax.plot(d240.time_h / 24, d240.conc_ng_mL, color=C240, lw=1.2, label="240 mg QD (calibrated F1)")
    ax.set_yscale("log"); ax.set_ylim(20, 15000); ax.set_xlabel("Day"); ax.set_ylabel("Plasma sotorasib (ng/mL)")
    ax.set_title("Typical patient, 22 days: autoinduction lowers both curves, the high dose more")
    ax.legend(); fig.tight_layout(); fig.savefig(FIG / "fig6_conc_time_typical.png", dpi=200); plt.close(fig)


if __name__ == "__main__":
    for f in [fig1_dose_exposure, fig2_exposure_overlap, fig3_ratio_timecourse, fig4_er_efficacy, fig5_tipping_point, fig6_conc_time]:
        f(); print("done", f.__name__)
