"""Prospective, hypothetical ORR planning scenarios and approximate precision."""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import norm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from design_scenario import power_normal_approximation

ROOT = Path(__file__).resolve().parents[1]


def sample_size_per_arm(p1, p2, alpha, target_power):
    """Equal allocation; two-sided pooled-null normal test; no attrition allowance."""
    values = np.array([p1, p2, alpha, target_power], dtype=float)
    if not np.isfinite(values).all() or not np.all((values > 0) & (values < 1)):
        raise ValueError("Probabilities, alpha, and target power must lie between 0 and 1.")
    if p1 == p2 or target_power <= alpha:
        raise ValueError("Specify different response rates and target power above alpha.")
    upper = 2.0
    while power_normal_approximation(p1, p2, upper, alpha) < target_power:
        upper *= 2
        if upper > 10_000_000:
            raise ValueError("Scenario requires more than 10 million participants per arm.")
    if upper == 2:
        return 2
    root = brentq(lambda n: power_normal_approximation(p1, p2, n, alpha) - target_power,
                  1, upper, xtol=1e-9)
    return int(np.ceil(root))


def precision_half_width(p1, p2, n, confidence=0.95):
    """Planning approximation for a risk-difference interval, using assumed rates.

    This is not an observed interval or a guarantee of realized interval width.
    """
    if not all(np.isfinite(v) for v in [p1, p2, n, confidence]):
        raise ValueError("Inputs must be finite.")
    if not 0 < p1 < 1 or not 0 < p2 < 1 or n < 2 or not 0 < confidence < 1:
        raise ValueError("Invalid planning rates, per-arm sample size, or confidence level.")
    return float(norm.ppf((1 + confidence) / 2) *
                 np.sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / n))


def build_tables():
    scenarios = pd.read_csv(ROOT / "data/design_scenarios.csv")
    if scenarios.scenario_id.duplicated().any():
        raise ValueError("Scenario identifiers must be unique.")
    rows = []
    for s in scenarios.itertuples(index=False):
        n = sample_size_per_arm(s.comparison_orr, s.reference_orr,
                                s.alpha_two_sided, s.target_power)
        rows.append({
            **s._asdict(),
            "difference_percentage_points": 100 * (s.comparison_orr - s.reference_orr),
            "n_per_arm": n, "n_total": 2 * n,
            "power_at_selected_n": power_normal_approximation(
                s.comparison_orr, s.reference_orr, n, s.alpha_two_sided),
            "power_at_one_fewer_per_arm": power_normal_approximation(
                s.comparison_orr, s.reference_orr, n - 1, s.alpha_two_sided),
        })
    sensitivity = pd.DataFrame(rows)
    # Fixed assumed rates isolate the effect of sample size on precision.
    # They are educational assumptions, not estimates from the sotorasib trial.
    precision = pd.DataFrame([
        {"reference_orr": 0.25, "comparison_orr": 0.35,
         "n_per_arm": n, "n_total": 2 * n, "confidence_level": 0.95,
         "approximate_half_width_percentage_points":
             100 * precision_half_width(0.35, 0.25, n)}
        for n in [50, 100, 200, 329, 500]
    ])
    return sensitivity, precision


def plot(sensitivity, precision):
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.8))
    example = sensitivity[np.isclose(sensitivity.reference_orr, 0.25)]
    for power, color in [(0.8, "#244665"), (0.9, "#a65c27")]:
        rows = example[np.isclose(example.target_power, power)].sort_values(
            "difference_percentage_points")
        axes[0].plot(rows.difference_percentage_points, rows.n_per_arm,
                     "o-", color=color, label=f"{power:.0%} power")
    axes[0].set(xlabel="Assumed ORR difference (percentage points)",
                ylabel="Participants per arm", title="Sample size depends on the assumed difference",
                xticks=[5, 10, 15], ylim=(0, None))
    axes[0].legend(frameon=False)
    axes[0].grid(axis="y", alpha=0.2)
    x = np.arange(len(precision))
    axes[1].bar(x, precision.approximate_half_width_percentage_points, color="#244665", width=.65)
    axes[1].set(xticks=x, xticklabels=precision.n_per_arm.astype(str),
                xlabel="Participants per arm", ylabel="Approximate 95% interval half-width (points)",
                title="Planning precision for a 35% vs 25% scenario", ylim=(0, 21))
    for i, row in enumerate(precision.itertuples(index=False)):
        axes[1].text(i, row.approximate_half_width_percentage_points + .4,
                     f"±{row.approximate_half_width_percentage_points:.1f}", ha="center", fontsize=9)
    fig.suptitle("Hypothetical response-rate design scenarios", fontsize=14)
    fig.text(.5, .015,
             "Equal allocation; normal approximations. Left: reference ORR 25%, two-sided alpha 0.05.\n"
             "Right: planning half-widths use assumed rates, not observed trial results. No attrition allowance.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=[0, .10, 1, .94], w_pad=2)
    for ext in ["png", "pdf"]:
        fig.savefig(ROOT / "figures" / f"design_sensitivity.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "figures").mkdir(exist_ok=True)
    sensitivity, precision = build_tables()
    sensitivity.to_csv(ROOT / "results/design_sensitivity.csv", index=False)
    precision.to_csv(ROOT / "results/design_precision.csv", index=False)
    plot(sensitivity, precision)
    print(f"Calculated {len(sensitivity)} hypothetical design scenarios and {len(precision)} precision examples.")


if __name__ == "__main__":
    main()
