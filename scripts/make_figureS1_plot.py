"""Plot Figure S1 from the synthetic Amazon theta summary.

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
import numpy as np
import pandas as pd


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
KEY = "theta"
GRID = [0.01, 0.1, 0.5, 1.0]
XLABEL = r"$\sigma_0$"
XLIM = (0.008, 1.15)
MISCOV_TICKS = [0.0001, 0.001, 0.01, 0.05, 0.1, 0.5]


def load_data(path: Path) -> pd.DataFrame:
    data = pd.read_csv(path)
    metrics = ["Avg_Coverage", "Avg_Length_Ratio"]
    required = [KEY, "Method", *metrics]
    missing = set(required) - set(data.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")

    data = data[data["Method"].isin([method for _, method, _, _ in METHODS])].copy()
    if not np.isfinite(data[[KEY, *metrics]].to_numpy()).all():
        raise ValueError(f"{path}: non-finite values")
    if not data["Avg_Coverage"].between(0, 1).all() or (data["Avg_Length_Ratio"] <= 0).any():
        raise ValueError(f"{path}: invalid coverage or width ratio")
    if data.duplicated([KEY, "Method"]).any():
        raise ValueError(f"{path}: duplicate setting/method rows")
    if len(data) != len(GRID) * len(METHODS):
        raise ValueError(f"{path}: unexpected number of plotted rows")

    for value in GRID:
        for _, method, _, _ in METHODS:
            n = len(data[np.isclose(data[KEY], value) & (data["Method"] == method)])
            if n != 1:
                raise ValueError(f"{path}: expected one {method} row at sigma_0={value}")

    data["miscoverage"] = 1 - data["Avg_Coverage"]
    if (data["miscoverage"] <= 0).any():
        raise ValueError(f"{path}: zero miscoverage cannot be plotted on a log axis")
    return data


def make_plot(data: pd.DataFrame, out: Path) -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 7,
        "axes.linewidth": 0.6, "pdf.fonttype": 42, "ps.fonttype": 42,
        "mathtext.fontset": "dejavusans",
    })
    fig, axes = plt.subplots(1, 2, figsize=(3.7, 2))
    fig.subplots_adjust(left=0.12, right=0.99, bottom=0.25, top=0.72, wspace=0.38)
    width_ax, miscov_ax = axes

    for label, method, color, marker in METHODS:
        d = data[data["Method"] == method].sort_values(KEY)
        marker_size = 2.8 if label == "RB NPMLE (ours)" else 3.5
        style = dict(color=color, marker=marker, markersize=marker_size,
                     linewidth=0.8, markeredgewidth=0.45)
        width_ax.plot(d[KEY], d["Avg_Length_Ratio"], label=label, **style)
        miscov_ax.plot(d[KEY], d["miscoverage"], **style)

    width_ax.set_ylabel("Average Width Ratio", labelpad=2)
    width_ax.set_ylim(0.15, 1.025)
    width_ax.set_yticks([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

    miscov_ax.axhline(0.05, color="black", linestyle="--", linewidth=0.6,
                      zorder=1, label=r"Nominal $\alpha$")
    miscov_ax.set_ylabel("Miscoverage Rate", labelpad=2)
    miscov_ax.set_yscale("log")
    miscov_ax.set_ylim(0.00004, 0.6)
    miscov_ax.set_yticks(MISCOV_TICKS, labels=[f"{value:g}" for value in MISCOV_TICKS])
    miscov_ax.yaxis.set_minor_formatter(NullFormatter())

    for ax in (width_ax, miscov_ax):
        ax.set_xlabel(XLABEL, labelpad=3)
        ax.set_xscale("log")
        ax.set_xlim(XLIM)
        ax.set_xticks(GRID, labels=[f"{value:g}" for value in GRID])
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.tick_params(axis="both", which="major", labelsize=6, length=2, pad=1.5)
        ax.tick_params(which="minor", length=0)
        ax.grid(which="major", color="#d9d9d9", linewidth=0.5)
        ax.grid(which="minor", color="#eeeeee", linewidth=0.35)
        ax.set_axisbelow(True)

    handles, labels = width_ax.get_legend_handles_labels()
    nominal_handle = miscov_ax.get_legend_handles_labels()[0][0]
    fig.legend(handles + [nominal_handle], labels + [r"Nominal $\alpha$"],
               loc="upper center", bbox_to_anchor=(0.52, 0.95), ncol=5,
               frameon=False, fontsize=5.3, handlelength=1.5,
               handletextpad=0.3, columnspacing=0.65)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=Path("amazon/results/synthetic_theta.csv"))
    parser.add_argument("--out", type=Path, default=Path("paper_plots/figureS1.pdf"))
    args = parser.parse_args()
    data = load_data(args.csv)
    make_plot(data, args.out)
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
