"""Rebuild paper Figure 3 from LMArena and Amazon summary results.

The original four-panel layout and method styling are retained, with two
additional EBCI baselines. As in the original figure, curves show means
without error bars. Run this script from the repository root.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullFormatter
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd
from eb_rebiasing.data_paths import LMARENA_RESULTS


# label, LMArena key, Amazon key, original color, marker
METHODS = [
    ("Classical", "mle_ci", "classic", "#1f77b4", "o"),
    ("Pred Mean", "pred_mean_ci", "classic_bias", "#bcbd22", "v"),
    ("PT", "pt_ci", "pt", "#d62728", "D"),
    ("EBCI Classical", "robust_eb_classical_ci", "classic_ebci", "#8c564b", "s"),
    ("EBCI PT", "robust_eb_pt_ci", "pt_ebci", "#829372", "^"),
    ("RB Normal (ours)", "eb_pt_ci", "pt_param", "#9467bd", "P"),
    ("RB NPMLE (ours)", "eb_npmle_pt_ci", "pt_npmle", "#e377c2", "*"),
]
ALPHAS = [0.05, 0.1, 0.2, 0.3]
MISCOV_TICKS = [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5]


def load_plot_data(lmarena_path: Path, amazon_path: Path) -> pd.DataFrame:
    lmarena = pd.read_csv(lmarena_path)
    amazon = pd.read_csv(amazon_path)
    if lmarena.duplicated(["alpha", "method"]).any():
        raise ValueError("Duplicate LMArena (alpha, method) rows")
    records = []
    for alpha in ALPHAS:
        for label, lm_key, am_key, _, _ in METHODS:
            lm_result = lmarena[
                np.isclose(lmarena["alpha"], alpha)
                & (lmarena["method"] == lm_key)
            ]
            if len(lm_result) != 1:
                raise ValueError(
                    f"Expected one LMArena row for {lm_key}, alpha={alpha}"
                )
            lm_row = lm_result.iloc[0]
            records.append(dict(
                dataset="LMArena", method=label, alpha=alpha,
                mean_coverage=lm_row["mean_coverage"],
                se_coverage=lm_row["se_coverage"],
                mean_width=lm_row["mean_width"],
                se_width=lm_row["se_width"],
                trials=np.nan,
            ))
            result = amazon[np.isclose(amazon["alpha"], alpha)
                            & (amazon["Method"] == am_key)]
            if len(result) != 1:
                raise ValueError(f"Expected one Amazon row for {am_key}, alpha={alpha}")
            r = result.iloc[0]
            records.append(dict(
                dataset="Amazon", method=label, alpha=alpha,
                mean_coverage=r["Avg_Coverage"], se_coverage=r["SE_Coverage"],
                mean_width=r["Avg_Length"], se_width=r["SE_Length"],
                trials=np.nan,  # Input summary does not encode trial count.
            ))
    result = pd.DataFrame(records)
    result["miscoverage"] = 1 - result["mean_coverage"]
    metrics = result[["mean_coverage", "se_coverage", "mean_width", "se_width"]]
    if not np.isfinite(metrics.to_numpy()).all():
        raise ValueError("Non-finite plotting metrics")
    if not result["mean_coverage"].between(0, 1).all() or (result["mean_width"] <= 0).any():
        raise ValueError("Invalid coverage or width")
    return result


def make_plot(data: pd.DataFrame, out: Path):
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 7,
        "axes.linewidth": 0.6, "pdf.fonttype": 42, "ps.fonttype": 42,
        "mathtext.fontset": "dejavusans",
    })
    fig, axes = plt.subplots(1, 4, figsize=(7.4, 1.65))
    fig.subplots_adjust(left=0.065, right=0.995, bottom=0.28, top=0.78, wspace=0.35)
    for dataset, panel_axes in zip(["LMArena", "Amazon"], [axes[:2], axes[2:]]):
        for label, _, _, color, marker in METHODS:
            d = data[(data["dataset"] == dataset) & (data["method"] == label)].sort_values("alpha")
            marker_size = 2.8 if label == "RB NPMLE (ours)" else 3.5
            style = dict(color=color, marker=marker, markersize=marker_size,
                         linewidth=0.8, markeredgewidth=0.45)
            panel_axes[0].plot(d["alpha"], d["mean_width"], label=label, **style)
            panel_axes[1].plot(d["alpha"], d["miscoverage"], **style)
        panel_axes[0].set_ylabel("Average Width", labelpad=2)
        panel_axes[0].set_ylim(0.04, 0.6)
        panel_axes[0].set_yticks([0.1, 0.3, 0.5])
        miscov_ax = panel_axes[1]
        miscov_ax.plot([0.04, 0.35], [0.04, 0.35], color="black",
                       linestyle="--", linewidth=0.6, label=r"Nominal $\alpha$")
        miscov_ax.set_ylabel("Miscoverage Rate", labelpad=2)
        miscov_ax.set_yscale("log")
        miscov_ax.set_ylim(0.01, 0.99)
        miscov_ax.set_yticks(MISCOV_TICKS, labels=[f"{v:g}" for v in MISCOV_TICKS])
        miscov_ax.yaxis.set_minor_formatter(NullFormatter())
        for ax in panel_axes:
            ax.set_xlabel(r"$\alpha$", labelpad=0)
            ax.set_xscale("log")
            ax.set_xlim(0.04, 0.35)
            ax.set_xticks(ALPHAS, labels=[f"{a:g}" for a in ALPHAS])
            # The final log-spaced ticks are close; give their labels room
            # while retaining the original horizontal tick text.
            ax.get_xticklabels()[-2].set_horizontalalignment("right")
            ax.xaxis.set_minor_locator(FixedLocator([0.06, 0.08, 0.15, 0.25]))
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.tick_params(axis="both", which="major", labelsize=6, length=2, pad=1.5)
            ax.tick_params(which="minor", length=0)
            ax.grid(which="major", color="#d9d9d9", linewidth=0.5)
            ax.grid(which="minor", color="#eeeeee", linewidth=0.35)
            ax.set_axisbelow(True)
    handles, labels = axes[0].get_legend_handles_labels()
    nominal_handle = axes[1].get_legend_handles_labels()[0][0]
    fig.legend(handles + [nominal_handle], labels + [r"Nominal $\alpha$"],
               loc="lower center", bbox_to_anchor=(0.52, 0.015), ncol=8,
               frameon=False, fontsize=5.7, handlelength=1.6,
               handletextpad=0.3, columnspacing=0.75)
    for i, label in enumerate(["(a) LMArena", "(b) Amazon"]):
        pair = axes[i * 2:i * 2 + 2]
        center = (pair[0].get_position().x0 + pair[1].get_position().x1) / 2
        fig.text(center, 0.88, label, ha="center", va="top", fontsize=7.5, fontweight="bold")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    content_bbox = fig.get_tightbbox(fig.canvas.get_renderer())
    output_bbox = Bbox.from_extents(
        0, 0, fig.get_figwidth(), min(fig.get_figheight(), content_bbox.y1 + 0.06)
    )
    fig.savefig(out, dpi=300, bbox_inches=output_bbox)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--lmarena-summary",
        type=Path,
        default=LMARENA_RESULTS / "summary.csv",
    )
    parser.add_argument(
        "--amazon-summary",
        type=Path,
        default=Path("amazon/results/amazon_real.csv"),
    )
    parser.add_argument("--out", type=Path, default=Path("paper_plots/figure3.pdf"))
    args = parser.parse_args()
    data = load_plot_data(args.lmarena_summary, args.amazon_summary)
    make_plot(data, args.out)
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
