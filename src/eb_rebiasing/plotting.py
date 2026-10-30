"""Descriptive bias diagnostics, with one point per task.

Scatter inputs use columns ``theta`` and ``bias``. Nuisance inputs use
``bias_hat``, ``rho``, ``tau2`` and ``sigma2``; ``bias`` is also accepted
as the observed-bias column. Variance columns contain variances, not SDs.
Nonfinite pairs are omitted separately in each panel.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from numpy.polynomial import Polynomial


NUISANCE_COLUMNS = (("rho", r"$\rho_i$"),
                    ("tau2", r"$\tau_i^2$"),
                    ("sigma2", r"$\sigma_i^2$"))


def _bias_column(summary: pd.DataFrame) -> str:
    return "bias_hat" if "bias_hat" in summary else "bias"


def _finite_pair(summary: pd.DataFrame, x_col: str, y_col: str):
    missing = {x_col, y_col} - set(summary.columns)
    if missing:
        raise ValueError(f"Missing summary columns: {sorted(missing)}")
    x = summary[x_col].to_numpy(dtype=float)
    y = summary[y_col].to_numpy(dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 2:
        raise ValueError(f"Need at least two finite pairs for {x_col} and {y_col}")
    return x[keep], y[keep]


def _correlation(x, y) -> float:
    if np.ptp(x) == 0 or np.ptp(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def _draw_scatter(
    ax, x, y, *, degree: int, fit: bool,
    point_size: float = 22, point_alpha: float = 0.55,
):
    if not isinstance(degree, (int, np.integer)) or degree < 0:
        raise ValueError("degree must be a nonnegative integer")
    ax.scatter(x, y, s=point_size, alpha=point_alpha, color="#1f77b4",
               edgecolors="none", rasterized=True)
    ax.axhline(0, color="0.4", linestyle="--", linewidth=1, zorder=1)
    if fit:
        order = min(degree, np.unique(x).size - 1)
        padding = 0.05 * max(abs(x[0]), 1) if np.ptp(x) == 0 else 0
        grid = np.linspace(x.min()-padding, x.max()+padding, 300)
        values = (np.full_like(grid, y.mean()) if order == 0
                  else Polynomial.fit(x, y, order)(grid))
        ax.plot(grid, values, color="#d62728", linewidth=2.2)
    ax.grid(True, color="0.9", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=10)


def plot_bias_scatter(
    summary: pd.DataFrame, *, dataset_label: str = "", ax=None,
    degree: int = 2, fit: bool | None = None, target_is_pseudo: bool = True,
):
    """Plot bias against target with an optional least-squares polynomial.

    Returns the Axes. The fit degree is limited by the number of distinct
    target values; a constant target produces a horizontal mean fit.
    """
    import matplotlib.pyplot as plt

    x, y = _finite_pair(summary, "theta", "bias")
    if ax is None:
        _, ax = plt.subplots(figsize=(5.2, 3.7), constrained_layout=True)
    _draw_scatter(ax, x, y, degree=degree, fit=fit)
    prefix = "(pseudo)-primary parameter" if target_is_pseudo else "Primary parameter"
    ax.set_xlabel(prefix + r" $\theta_i$", fontsize=12)
    bias_prefix = "(pseudo)-" if target_is_pseudo else ""
    ax.set_ylabel(bias_prefix + r"prediction bias $b_i$", fontsize=12)
    title = "Exchangeability sanity check"
    ax.set_title(title + (f" for {dataset_label}" if dataset_label else ""), fontsize=14)
    return ax


def nuisance_correlations(summary: pd.DataFrame) -> pd.DataFrame:
    """Return per-panel finite-pair counts and Pearson correlations."""
    rows = []
    for column, _ in NUISANCE_COLUMNS:
        x, y = _finite_pair(summary, column, _bias_column(summary))
        rows.append({"parameter": column, "n": len(x), "correlation": _correlation(x, y)})
    return pd.DataFrame(rows)


def plot_bias_vs_nuisance(
    summary: pd.DataFrame, *, dataset_label: str = "", axes=None,
    degree: int = 2, fit: bool | None = None,
    point_size: float = 22, point_alpha: float = 0.55,
):
    """Plot observed bias against rho, tau-squared and sigma-squared.

    Here sigma2 = Var(theta_b), tau2 = Var(bias_hat), and
    rho = Corr(theta_b, bias_hat). Returns the three Axes. A polynomial fit
    of the requested degree is drawn only when ``fit`` is true.
    """
    import matplotlib.pyplot as plt

    if axes is None:
        _, axes = plt.subplots(1, 3, figsize=(11.4, 3.6), constrained_layout=True)
    axes = np.asarray(axes, dtype=object).reshape(-1)
    if len(axes) != 3 or any(ax.figure is not axes[0].figure for ax in axes):
        raise ValueError("axes must contain three Axes from the same figure")
    for ax, (column, label) in zip(axes, NUISANCE_COLUMNS):
        x, y = _finite_pair(summary, column, _bias_column(summary))
        _draw_scatter(ax, x, y, degree=degree, fit=fit,
                      point_size=point_size, point_alpha=point_alpha)
        ax.set_xlabel(label, fontsize=12)
    axes[0].set_ylabel(r"bias estimate $\hat b_i$", fontsize=12)
    if dataset_label:
        axes[0].set_title(dataset_label, loc="left", fontsize=12, fontweight="bold")
    return axes


__all__ = ["plot_bias_scatter", "plot_bias_vs_nuisance", "nuisance_correlations"]
