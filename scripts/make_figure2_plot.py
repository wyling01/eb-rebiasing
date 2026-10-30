"""Render Figure 2 from the LMArena, Amazon, and GWAS plotting data.

Run from the tide_code repository root. Each panel reads from its own results
directory.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import norm

LABEL_FS = 13
TICK_FS = 11
LEGEND_FS = 13
COLTITLE_FS = 14
SAVE_DPI = 150

PARAM_COLOR = "firebrick"
NPMLE_COLOR = "steelblue"

# Per-panel tuning (histogram bins / spike-height scaling) carried over
# from the original colleague notebook so the visuals stay comparable.
# ``param_spike_height`` only fires for GWAS, where the fitted Gaussian
# prior has near-zero variance and is best rendered as a single spike.
# ``sample_n`` down-samples the bias vector before computing the
# expensive per-row marginal density (only relevant for GWAS).
PANEL_CONFIG = [
    dict(
        key="lmarena",  title="(a) LMArena",
        bias_csv="lmarena_bias_se_pt.csv",
        npmle_csv="lmarena_npmle_pt.csv",
        param_csv="lmarena_param_pt.csv",
        hist_bins=25, npmle_spike_scale=5.0,
        param_spike_height=None, sample_n=None,
        gamma_check=0.991,
    ),
    dict(
        key="amazon",   title="(b) Amazon",
        bias_csv="amazon_bias_se_pt.csv",
        npmle_csv="amazon_npmle_pt.csv",
        param_csv="amazon_param_pt.csv",
        hist_bins=20, npmle_spike_scale=20.0,
        param_spike_height=None, sample_n=None,
        gamma_check=0.826,
    ),
    dict(
        key="gwas",     title="(c) GWAS",
        bias_csv="gwas_bias_se_pt.csv",
        npmle_csv="gwas_npmle_pt.csv",
        param_csv="gwas_param_pt.csv",
        hist_bins=25, npmle_spike_scale=25.0,
        param_spike_height=25.0, sample_n=10_000,
        gamma_check=0.752,
    ),
]


# ----------------------------------------------------------------------------
# Panel renderer
# ----------------------------------------------------------------------------

def _draw_panel(ax, *, bias_csv, npmle_csv, param_csv,
                hist_bins, npmle_spike_scale, param_spike_height,
                sample_n, title, gamma_check=None, show_legend=False):
    bias_df = pd.read_csv(bias_csv)
    npmle_df = pd.read_csv(npmle_csv)
    param_df = pd.read_csv(param_csv)

    bias_full = bias_df["bias"].to_numpy(dtype=float)
    bias_sd_full = bias_df["bias_sd"].to_numpy(dtype=float)
    valid = (np.isfinite(bias_full) & np.isfinite(bias_sd_full)
             & (bias_sd_full > 0))
    bias_full = bias_full[valid]
    bias_sd_full = bias_sd_full[valid]

    # Down-sample large bias vectors so the per-row marginal-density
    # accumulation stays cheap; histogram still uses the full vector.
    if sample_n is not None and len(bias_full) > sample_n:
        rng = np.random.default_rng(123)
        idx = rng.choice(len(bias_full), size=sample_n, replace=False)
        bias = bias_full[idx]
        bias_sd = bias_sd_full[idx]
    else:
        bias = bias_full
        bias_sd = bias_sd_full

    mu_b = float(param_df["hat_mu_b"].iloc[0])
    var_b = float(param_df["hat_var_b"].iloc[0])
    sd_b = float(np.sqrt(var_b))

    atoms = npmle_df["atom"].to_numpy(dtype=float)
    weights = npmle_df["weight"].to_numpy(dtype=float)
    weights = weights / weights.sum()

    x_min = min(np.nanmin(bias), mu_b - 4 * sd_b, np.nanmin(atoms))
    x_max = max(np.nanmax(bias), mu_b + 4 * sd_b, np.nanmax(atoms))
    pad = 0.1 * (x_max - x_min)
    x_grid = np.linspace(x_min - pad, x_max + pad, 1000)

    # Parametric marginal: average over i of N(mu, var_b + bias_sd_i^2).
    param_marginal = np.mean(
        norm.pdf(x_grid[:, None], loc=mu_b,
                 scale=np.sqrt(var_b + bias_sd[None, :] ** 2)),
        axis=1,
    )
    # NPMLE marginal: average over i of the mixture sum_k w_k N(atom_k, bias_sd_i^2).
    npmle_marginal = np.zeros_like(x_grid)
    for sd_i in bias_sd:
        npmle_marginal += np.sum(
            weights[None, :] * norm.pdf(
                x_grid[:, None], loc=atoms[None, :], scale=sd_i),
            axis=1,
        )
    npmle_marginal /= len(bias_sd)

    param_prior_density = norm.pdf(x_grid, loc=mu_b, scale=sd_b)
    npmle_prior_spike_height = weights * npmle_spike_scale

    # Observed bias histogram (uses the full, un-sampled vector).
    ax.hist(bias_full, bins=hist_bins, density=True, color="gray",
            alpha=0.18, edgecolor="gray", linewidth=0.6)
    # Marginal densities of $\hat b_i$ -- dotted.
    ax.plot(x_grid, param_marginal, color=PARAM_COLOR,
            linewidth=1.6, linestyle=":")
    ax.plot(x_grid, npmle_marginal, color=NPMLE_COLOR,
            linewidth=1.6, linestyle=":")
    # Priors on $b_i$ -- solid for parametric, vlines for NPMLE atoms.
    ax.plot(x_grid, param_prior_density, color=PARAM_COLOR,
            linewidth=2.2, linestyle="-")
    if param_spike_height is not None:
        ax.vlines(mu_b, ymin=0, ymax=param_spike_height,
                  color=PARAM_COLOR, linewidth=2.2, alpha=0.9)
    ax.vlines(atoms, ymin=0, ymax=npmle_prior_spike_height,
              color=NPMLE_COLOR, linewidth=2.2, alpha=0.9)

    ax.set_xlabel(r"$\hat b_i$ or $b_i$", fontsize=LABEL_FS)
    ax.set_ylabel("Density", fontsize=LABEL_FS)
    ax.tick_params(axis="both", labelsize=TICK_FS)
    ax.set_title(title, fontsize=COLTITLE_FS, fontweight="bold", pad=4)
    ax.grid(True, alpha=0.25, which="both")

    if show_legend:
        legend_handles = [
            Line2D([0], [0], color=PARAM_COLOR, lw=2.2, label="Normal"),
            Line2D([0], [0], color=NPMLE_COLOR, lw=2.2, label="NPMLE"),
        ]
        ax.legend(handles=legend_handles, frameon=False,
                  fontsize=LEGEND_FS, handletextpad=0.4, loc="upper left",
                  bbox_to_anchor=(-0.02, 1.0))

    # Top-right: maximum absolute estimator/bias correlation
    # gamma-check, set in axes-fraction coords so it sits a small fixed
    # margin from the top/right edges regardless of data range.
    if gamma_check is not None:
        ax.text(
            0.96, 0.96, rf"$\overset{{\vee}}{{\gamma}} = {gamma_check:.3f}$",
            transform=ax.transAxes,
            ha="right", va="top",
            fontsize=LABEL_FS + 2, color="black",
        )


# ----------------------------------------------------------------------------
# Figure assembly
# ----------------------------------------------------------------------------

def make_plot(data_dirs: dict[str, Path], out_path: Path):
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.3))

    for i, (ax, cfg) in enumerate(zip(axes, PANEL_CONFIG)):
        data_dir = data_dirs[cfg["key"]]
        _draw_panel(
            ax,
            bias_csv=data_dir / cfg["bias_csv"],
            npmle_csv=data_dir / cfg["npmle_csv"],
            param_csv=data_dir / cfg["param_csv"],
            hist_bins=cfg["hist_bins"],
            npmle_spike_scale=cfg["npmle_spike_scale"],
            param_spike_height=cfg["param_spike_height"],
            sample_n=cfg["sample_n"],
            title=cfg["title"],
            gamma_check=cfg.get("gamma_check"),
            show_legend=(i == 0),
        )

    fig.tight_layout(pad=1.0)

    fig.savefig(out_path, bbox_inches="tight", dpi=SAVE_DPI)
    plt.close(fig)
    print(f"Saved plot -> {out_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lmarena-dir", type=Path,
                        default=Path("tmp/lmarena_plotting"))
    parser.add_argument("--amazon-dir", type=Path,
                        default=Path("amazon/results"))
    parser.add_argument("--gwas-dir", type=Path,
                        default=Path("gwas/results"))
    parser.add_argument("--out", type=Path, default=Path("paper_plots/figure2.pdf"))
    args = parser.parse_args()

    data_dirs = {
        "lmarena": args.lmarena_dir,
        "amazon": args.amazon_dir,
        "gwas": args.gwas_dir,
    }
    missing = [data_dirs[cfg["key"]] / cfg[key] for cfg in PANEL_CONFIG
               for key in ("bias_csv", "npmle_csv", "param_csv")
               if not (data_dirs[cfg["key"]] / cfg[key]).is_file()]
    if missing:
        parser.error("Missing input CSVs: " + ", ".join(map(str, missing)))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    make_plot(data_dirs, args.out)


if __name__ == "__main__":
    main()
