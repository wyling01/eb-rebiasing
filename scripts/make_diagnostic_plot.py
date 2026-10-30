"""Render the exchangeability and nuisance-parameter diagnostic plots."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from eb_rebiasing.data_paths import LMARENA_CSV
from eb_rebiasing.diagnostics import load_lmarena_nuisance_summary, load_lmarena_pseudo_oracle
from eb_rebiasing.plotting import plot_bias_scatter, plot_bias_vs_nuisance


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lmarena-csv", type=Path, default=LMARENA_CSV)
    parser.add_argument("--amazon-csv", type=Path,
                        default=Path("amazon/results/amazon_bias_nuisance_diagnostics.csv"))
    parser.add_argument("--gwas-csv", type=Path,
                        default=Path("gwas/results/height_gwas_bias_nuisance_diagnostics.csv"))
    parser.add_argument("--train-test-split", type=float, default=0.2)
    parser.add_argument("--split-seed", type=int, default=42)
    parser.add_argument("--pseudo-seed", type=int, default=42)
    parser.add_argument("--fitted-line", action="store_true", default=None)
    parser.add_argument("--out", type=Path,
                        default=Path("paper_plots/diagnostic_pseudo.pdf"))
    parser.add_argument("--nuisance-out", type=Path,
                        default=Path("paper_plots/diagnostic_nuisance.pdf"))
    args = parser.parse_args()

    lmarena = load_lmarena_pseudo_oracle(
        args.lmarena_csv, train_test_split=args.train_test_split,
        split_seed=args.split_seed, pseudo_seed=args.pseudo_seed,
    )
    amazon = pd.read_csv(args.amazon_csv)
    gwas = pd.read_csv(args.gwas_csv)
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.0), constrained_layout=True)
    fig.set_constrained_layout_pads(wspace=0.08)
    for ax, data, title in zip(axes, [lmarena, amazon], ["(a) LMArena", "(b) Amazon"]):
        plot_bias_scatter(data, ax=ax, fit=args.fitted_line, target_is_pseudo=True)
        ax.set_title(title, fontsize=12, fontweight="bold")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=300)
    plt.close(fig)

    lmarena_nuisance = load_lmarena_nuisance_summary(
        args.lmarena_csv, train_test_split=args.train_test_split,
        split_seed=args.split_seed,
    )
    nuisance_data = [lmarena_nuisance, amazon, gwas]
    row_titles = ["(a) LMArena", "(b) Amazon", "(c) GWAS"]
    nuisance_fig, nuisance_axes = plt.subplots(
        3, 3, figsize=(10.5, 8.4), constrained_layout=True,
    )
    nuisance_fig.set_constrained_layout_pads(hspace=0.02, wspace=0.05)
    for row_axes, data, title in zip(nuisance_axes, nuisance_data, row_titles):
        plot_bias_vs_nuisance(
            data, dataset_label=title, axes=row_axes,
            degree=2, fit=args.fitted_line,
            point_size=1 if title == "(c) GWAS" else 22,
            point_alpha=0.08 if title == "(c) GWAS" else 0.55,
        )
    nuisance_fig.canvas.draw()
    tau_ticks = [
        ([0.000, 0.005, 0.010, 0.015], ["0.000", "0.005", "0.010", "0.015"]),
        ([0.00, 0.01, 0.02, 0.03], ["0.00", "0.01", "0.02", "0.03"]),
        ([0.000, 0.001, 0.002, 0.003], ["0.000", "0.001", "0.002", "0.003"]),
    ]
    for row_axes, (ticks, labels) in zip(nuisance_axes, tau_ticks):
        row_axes[1].set_xticks(ticks, labels=labels)
    nuisance_axes[2, 2].set_xticks(
        [0.000, 0.001, 0.002, 0.003],
        labels=["0.000", "0.001", "0.002", "0.003"],
    )
    args.nuisance_out.parent.mkdir(parents=True, exist_ok=True)
    nuisance_fig.savefig(args.nuisance_out, dpi=300)
    plt.close(nuisance_fig)

    print(f"Saved {args.out} ({len(lmarena)} LMArena tasks; {len(amazon)} Amazon tasks)")
    print(f"Saved {args.nuisance_out} ({len(lmarena_nuisance)} LMArena tasks; "
          f"{len(amazon)} Amazon tasks; {len(gwas)} GWAS tasks)")


if __name__ == "__main__":
    main()
