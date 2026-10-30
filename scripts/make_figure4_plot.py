"""Plot Figure 4 from the supplied synthetic Amazon summaries.

Run from the repository root. Uses the CSV's average per-task width ratios
directly and plots the Monte Carlo means without error bars.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd


# Keep baseline/reference methods first and the proposed methods last.
METHODS = [
    ("Classical", "classic", "#1f77b4", "o"),
    ("Pred Mean", "classic_bias", "#bcbd22", "v"),
    ("PT", "pt", "#d62728", "D"),
    ("EBCI Classical", "classic_ebci", "#8c564b", "s"),
    ("EBCI PT", "pt_ebci", "#829372", "^"),
    ("Oracle (ours)", "pt_oracle", "#7f7f7f", "x"),
    ("RB Normal (ours)", "pt_param", "#9467bd", "P"),
    ("RB NPMLE (ours)", "pt_npmle", "#e377c2", "*"),
]
PANELS = [
    ("Normal", "sigma", [0.01, 0.03, 0.05, 0.1], r"$A^{1/2}$", (0.009, 0.11)),
    ("Two-point", "b0", [0.05, 0.1, 0.2, 0.5], r"$b_0$", (0.045, 0.56)),
]
MISCOV_TICKS = [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5]


def load_panel(path: Path, key: str, grid: list[float]) -> pd.DataFrame:
    data = pd.read_csv(path)
    metrics = ["Avg_Coverage", "Avg_Length_Ratio"]
    required = [key, "Method", *metrics]
    missing = set(required) - set(data.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    data = data[data["Method"].isin([m[1] for m in METHODS])].copy()
    if not np.isfinite(data[[key, *metrics]].to_numpy()).all():
        raise ValueError(f"{path}: non-finite values")
    if not data["Avg_Coverage"].between(0, 1).all() or (data["Avg_Length_Ratio"] <= 0).any():
        raise ValueError(f"{path}: invalid coverage or width ratio")
    if data.duplicated([key, "Method"]).any():
        raise ValueError(f"{path}: duplicate setting/method rows")
    if len(data) != len(grid) * len(METHODS):
        raise ValueError(f"{path}: unexpected number of plotted rows")
    for value in grid:
        for _, method, _, _ in METHODS:
            if len(data[np.isclose(data[key], value) & (data["Method"] == method)]) != 1:
                raise ValueError(f"{path}: expected one {method} row at {key}={value}")
    data["miscoverage"] = 1 - data["Avg_Coverage"]
    if (data["miscoverage"] <= 0).any():
        raise ValueError(f"{path}: zero miscoverage cannot be plotted on a log axis")
    return data


def make_plot(normal: pd.DataFrame, twopoint: pd.DataFrame, out: Path):
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 7,
        "axes.linewidth": 0.6, "pdf.fonttype": 42, "ps.fonttype": 42,
        "mathtext.fontset": "dejavusans",
    })
    fig, axes = plt.subplots(1, 4, figsize=(7.4, 1.65))
    fig.subplots_adjust(left=0.07, right=0.995, bottom=0.28, top=0.78, wspace=0.38)
    for i, (data, panel) in enumerate(zip([normal, twopoint], PANELS)):
        _, key, grid, xlabel, xlim = panel
        width_ax, miscov_ax = axes[i * 2:i * 2 + 2]
        for label, method, color, marker in METHODS:
            d = data[data["Method"] == method].sort_values(key)
            marker_size = 2.8 if label == "RB NPMLE (ours)" else 3.5
            style = dict(color=color, marker=marker, markersize=marker_size,
                         linewidth=0.8, markeredgewidth=0.45)
            width_ax.plot(d[key], d["Avg_Length_Ratio"], label=label, **style)
            miscov_ax.plot(d[key], d["miscoverage"], **style)
        width_ax.set_ylabel("Average Width Ratio", labelpad=2)
        width_ax.set_ylim(0.45, 1.025)
        width_ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        miscov_ax.axhline(0.05, color="black", linestyle="--", linewidth=0.6,
                          zorder=1, label=r"Nominal $\alpha$")
        miscov_ax.set_ylabel("Miscoverage Rate", labelpad=2)
        miscov_ax.set_yscale("log")
        miscov_ax.set_ylim(0.01, 0.99)
        miscov_ax.set_yticks(MISCOV_TICKS, labels=[f"{v:g}" for v in MISCOV_TICKS])
        miscov_ax.yaxis.set_minor_formatter(NullFormatter())
        for ax in [width_ax, miscov_ax]:
            ax.set_xlabel(xlabel, labelpad=3)
            ax.set_xscale("log")
            ax.set_xlim(xlim)
            ax.set_xticks(grid, labels=[f"{v:g}" for v in grid])
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.tick_params(axis="both", which="major", labelsize=6, length=2, pad=1.5)
            ax.tick_params(which="minor", length=0)
            ax.grid(which="major", color="#d9d9d9", linewidth=0.5)
            ax.grid(which="minor", color="#eeeeee", linewidth=0.35)
            ax.set_axisbelow(True)
    handles, labels = axes[0].get_legend_handles_labels()
    nominal_handle = axes[1].get_legend_handles_labels()[0][0]
    fig.legend(handles + [nominal_handle], labels + [r"Nominal $\alpha$"],
               loc="lower center", bbox_to_anchor=(0.52, 0.015), ncol=9,
               frameon=False, fontsize=5.3, handlelength=1.5,
               handletextpad=0.3, columnspacing=0.65)
    for i, label in enumerate(["(a) Normal", "(b) Two-point"]):
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
        "--normal-csv",
        type=Path,
        default=Path("amazon/results/synthetic_normal.csv"),
    )
    parser.add_argument(
        "--twopoint-csv",
        type=Path,
        default=Path("amazon/results/synthetic_twopoint.csv"),
    )
    parser.add_argument("--out", type=Path, default=Path("paper_plots/figure4.pdf"))
    args = parser.parse_args()
    normal = load_panel(args.normal_csv, PANELS[0][1], PANELS[0][2])
    twopoint = load_panel(args.twopoint_csv, PANELS[1][1], PANELS[1][2])
    make_plot(normal, twopoint, args.out)
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
