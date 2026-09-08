"""
Component D: decision analysis.

D1. Post hoc power for the observed CodeBreaK 100 part B ORR difference, and prospective sample-size scenarios.
    Inputs (Hochmair 2024 Table 2): 34/104 vs 26/105.
    Amgen (Discussion): formal testing "would require enrollment of over 600 patients."
    FDA (Singh/Vellanki/Pazdur 2025, per public summaries): >500 under OS HR<=0.75 or
    ORR difference >=10%.
    Method: two-sided two-proportion z-test (pooled), alpha 0.05, plus Monte Carlo simulation of the same z-test.

D2. Exposure-efficacy reanalysis on continuous exposure.
    FDA Fig. 26 gives ORR by AUCtau,ss quartile: 25/57, 26/57, 21/57, 12/57 with bin
    ranges 13000-19000, 19000-27000, 27000-46000, 46000-80000 h*ng/mL.
    Individual data are not public. What CAN be done: fit a logistic model on log(AUC)
    to the binned data at the retained representative AUC values (13k/19k/27k/46k), with the observed counts. This puts the
    efficacy analysis on the same functional footing as FDA's safety analysis (continuous
    logistic), and yields a slope with a CI rather than four bars. It CANNOT adjust for
    covariates without patient-level data; that limitation is stated.

D3. Tipping-point analysis.
    Given the simulated exposure distributions for 240 and 960 mg (scenario B) and a
    hypothesized TRUE exposure-efficacy logistic slope (beta, per log-unit AUC), what ORR
    difference would the two doses produce? At what beta does the difference reach the
    observed 7.9 points? And at what beta does a benefit-risk utility
        U = P(response) - w * P(Grade>=3 TEAE)
    favor 240 mg? Grade>=3 exposure-toxicity slope is calibrated so that the simulated
    arms reproduce the observed 61.5% vs 49.0% (Hochmair Table 3).
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import brentq
import statsmodels.api as sm

ROOT = Path(__file__).resolve().parents[1]
RNG = np.random.default_rng(7)


# ---------------------------------------------------------------- D1: power
def power_two_prop(p1, p2, n1, n2, alpha=0.05):
    pbar = (p1 * n1 + p2 * n2) / (n1 + n2)
    se0 = np.sqrt(pbar * (1 - pbar) * (1 / n1 + 1 / n2))
    se1 = np.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
    z = stats.norm.ppf(1 - alpha / 2)
    d = abs(p1 - p2)
    return stats.norm.cdf((d - z * se0) / se1) + stats.norm.cdf((-d - z * se0) / se1)


def n_per_arm_for_power(p1, p2, power=0.8, alpha=0.05):
    f = lambda n: power_two_prop(p1, p2, n, n, alpha) - power
    return int(np.ceil(brentq(f, 10, 100000)))


def power_sim(p1, p2, n1, n2, alpha=0.05, nsim=20000):
    x1 = RNG.binomial(n1, p1, nsim)
    x2 = RNG.binomial(n2, p2, nsim)
    ph1, ph2 = x1 / n1, x2 / n2
    pbar = (x1 + x2) / (n1 + n2)
    se = np.sqrt(pbar * (1 - pbar) * (1 / n1 + 1 / n2))
    z = (ph1 - ph2) / np.where(se > 0, se, np.inf)
    return np.mean(np.abs(z) > stats.norm.ppf(1 - alpha / 2))


def run_power():
    p960, p240 = 34 / 104, 26 / 105
    out = {
        "interpretation": "Post hoc power treats the observed ORRs as true probabilities; it adds no independent evidence beyond the observed effect and interval. Sample-size scenarios assume equal allocation and two-sided alpha 0.05.",
        "observed": dict(p960=p960, p240=p240, diff=p960 - p240, n960=104, n240=105),
        "power_as_conducted_for_observed_diff": dict(
            normal_approx=float(power_two_prop(p960, p240, 104, 105)),
            simulation=float(power_sim(p960, p240, 104, 105))),
        "power_as_conducted_for_10pt_diff_35_vs_25": float(power_two_prop(0.35, 0.25, 104, 105)),
        "n_per_arm_80pct_power_observed_diff": n_per_arm_for_power(p960, p240),
        "n_per_arm_80pct_power_10pt_diff_35_vs_25": n_per_arm_for_power(0.35, 0.25),
        "n_per_arm_80pct_power_10pt_diff_37_vs_27": n_per_arm_for_power(0.37, 0.27),
        "amgen_stated_total_required": ">600",
        "fda_stated_total_required": ">500 (ORR diff >=10% or OS HR <=0.75)",
    }
    out["n_total_80pct_power_observed_diff"] = 2 * out["n_per_arm_80pct_power_observed_diff"]
    out["n_total_80pct_power_10pt_diff_35_vs_25"] = 2 * out["n_per_arm_80pct_power_10pt_diff_35_vs_25"]
    # Risk difference CI (Wald, unstratified) for comparison with published 7.9 (-2.3, 18.2) at 90%
    d = p960 - p240
    se = np.sqrt(p960 * (1 - p960) / 104 + p240 * (1 - p240) / 105)
    out["risk_diff_90CI_unstratified_reproduced"] = [float(d - 1.645 * se), float(d + 1.645 * se)]
    out["risk_diff_90CI_published_unstratified"] = [-2.3 / 100, 18.2 / 100]
    out["risk_diff_95CI"] = [float(d - 1.96 * se), float(d + 1.96 * se)]
    return out


# ------------------------------------------------ D2: exposure-efficacy, continuous
FDA_QUARTILES = pd.DataFrame(dict(
    q=[1, 2, 3, 4],
    auc_lo=[13000, 19000, 27000, 46000], auc_hi=[19000, 27000, 46000, 80000],
    auc_label=[13000, 19000, 27000, 46000],   # retained regression coordinates; not geometric bin midpoints
    responders=[25, 26, 21, 12], n=[57, 57, 57, 57]))


def run_er_efficacy():
    d = FDA_QUARTILES.copy()
    d["log_auc"] = np.log(d.auc_label)
    X = sm.add_constant(d.log_auc)
    m = sm.GLM(np.column_stack([d.responders, d.n - d.responders]), X,
               family=sm.families.Binomial()).fit()
    beta, se = m.params.iloc[1], m.bse.iloc[1]
    # Odds ratio per doubling of exposure
    or_per_doubling = np.exp(beta * np.log(2))
    ci = np.exp((beta + np.array([-1.96, 1.96]) * se) * np.log(2))
    return dict(logistic_slope_per_log_auc=float(beta), slope_se=float(se), slope_p=float(m.pvalues.iloc[1]),
                OR_per_doubling_of_AUC=float(or_per_doubling), OR_per_doubling_95CI=[float(ci[0]), float(ci[1])],
                intercept=float(m.params.iloc[0]),
                note="Binned data using representative AUC values 13000/19000/27000/46000 (not geometric bin midpoints); no patient-level covariate adjustment. The negative association may reflect confounding and does not establish a causal exposure effect.")


# ---------------------------------------------------- D3: tipping point
def run_tipping(er):
    pop = pd.read_csv(ROOT / "results" / "population_exposures.csv")
    a240 = pop[pop.scenario == "B_calibrated_240"].AUC_tau.values
    a960 = pop[pop.scenario == "960"].AUC_tau.values
    la240, la960 = np.log(a240), np.log(a960)
    la_ref = np.log(23000.0)  # anchor: 960 mg typical AUC; keeps overall ORR ~ observed

    def orr(la, beta, p_ref=0.327):
        a = np.log(p_ref / (1 - p_ref))
        return np.mean(1 / (1 + np.exp(-(a + beta * (la - la_ref)))))

    # Grade>=3 toxicity slope calibrated to observed 61.5% (960) vs 49.0% (240)
    def tox(la, gamma, p_ref=0.615):
        a = np.log(p_ref / (1 - p_ref))
        return np.mean(1 / (1 + np.exp(-(a + gamma * (la - la_ref)))))
    gamma = brentq(lambda g: tox(la240, g) - 0.49, 0.0, 10.0)

    betas = np.linspace(-1.5, 3.0, 91)
    diff = np.array([orr(la960, b) - orr(la240, b) for b in betas])
    # beta at which the dose difference reproduces the observed 7.9-point ORR gap
    beta_obs = brentq(lambda b: (orr(la960, b) - orr(la240, b)) - 0.079, 0.0, 5.0)
    # OR per doubling implied by that beta
    or_needed = np.exp(beta_obs * np.log(2))

    # Utility surface: U = P(resp) - w * P(G3+); find w threshold where 240 wins, for each beta
    ws = np.linspace(0, 1.0, 51)
    surface = np.zeros((len(betas), len(ws)))
    for i, b in enumerate(betas):
        e960, e240 = orr(la960, b), orr(la240, b)
        t960, t240 = tox(la960, gamma), tox(la240, gamma)
        for j, w in enumerate(ws):
            surface[i, j] = (e240 - w * t240) - (e960 - w * t960)   # >0 => 240 preferred
    pd.DataFrame(surface, index=np.round(betas, 3), columns=np.round(ws, 3)).to_csv(
        ROOT / "results" / "utility_surface_240_minus_960.csv")
    pd.DataFrame(dict(beta=betas, orr_diff_960_minus_240=diff)).to_csv(
        ROOT / "results" / "orr_diff_vs_er_slope.csv", index=False)

    # For the FDA-fitted (confounded) slope and for a null slope, what does the model predict?
    return dict(
        tox_slope_gamma_calibrated=float(gamma),
        tox_960=float(tox(la960, gamma)), tox_240=float(tox(la240, gamma)),
        beta_required_for_observed_7p9pt_ORR_gap=float(beta_obs),
        OR_per_doubling_required=float(or_needed),
        FDA_fitted_slope=er["logistic_slope_per_log_auc"],
        orr_diff_at_FDA_fitted_slope=float(orr(la960, er["logistic_slope_per_log_auc"]) - orr(la240, er["logistic_slope_per_log_auc"])),
        orr_diff_at_null_slope=float(orr(la960, 0.0) - orr(la240, 0.0)),
        toxicity_weight_at_which_240_wins_if_beta_is_required_value=float(
            (orr(la960, beta_obs) - orr(la240, beta_obs)) / (tox(la960, gamma) - tox(la240, gamma))),
        note="Exposure distributions from scenario B (calibrated). ORR anchored to 32.7% at 960 mg typical exposure.")


def main():
    out = dict(power=run_power(), er_efficacy=run_er_efficacy())
    out["tipping_point"] = run_tipping(out["er_efficacy"])
    (ROOT / "results" / "decision_analysis.json").write_text(json.dumps(out, indent=2, default=float))
    print(json.dumps(out, indent=2, default=float))


if __name__ == "__main__":
    main()
