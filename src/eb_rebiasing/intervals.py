"""Rebiased confidence intervals.

Four methods, one per ``(prior_kind, ppi_variant)`` combination:

- ``get_eb_ppi_cis``         : PPI (lambda=1)        × Gaussian prior
- ``get_eb_pt_cis``          : PT (per-problem)      × Gaussian prior
- ``get_npmle_eb_ppi_cis``   : PPI                   × NPMLE prior
- ``get_npmle_eb_pt_cis``    : PT                    × NPMLE prior

Each function takes a ``PasDataset``-shaped object and returns
``np.ndarray`` of shape ``(M, 2)`` with ``[lower, upper]`` columns,
matching the upstream CI signature so they can be used with
``pas.experiments.run_ci_benchmark``.
"""
from __future__ import annotations

from typing import Tuple, Union

import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm

from eb_rebiasing._ppi_helpers import (
    ppi_decomposition,
    pt_lambdas,
    pt_decomposition,
)
from eb_rebiasing.prior_fit.parametric import fit_parametric_bias_prior
from eb_rebiasing.prior_fit import _HAS_NPMLE

if _HAS_NPMLE:
    from eb_rebiasing.prior_fit.npmle import (
        fit_npmle_bias_prior,
        posterior_weights_npmle,
    )


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def _zconfint(
    mean: np.ndarray,
    se: np.ndarray,
    alpha: float,
    alternative: str,
) -> np.ndarray:
    """Symmetric Wald interval. Mirrors ``pas.utils._zconfint``."""
    M = mean.shape[0]
    out = np.empty((M, 2))
    if alternative == "two-sided":
        z = norm.ppf(1 - alpha / 2)
        out[:, 0] = mean - z * se
        out[:, 1] = mean + z * se
    elif alternative == "larger":
        z = norm.ppf(1 - alpha)
        out[:, 0] = mean - z * se
        out[:, 1] = np.inf
    elif alternative == "smaller":
        z = norm.ppf(1 - alpha)
        out[:, 0] = -np.inf
        out[:, 1] = mean + z * se
    else:
        raise ValueError(
            "alternative must be 'two-sided', 'larger', or 'smaller'; "
            f"got {alternative!r}"
        )
    return out


def _select_eb_subset(
    point: np.ndarray, bias: np.ndarray, cov_mats: np.ndarray, corr_threshold: float,
    *, allow_zero_variance: bool = False,
) -> np.ndarray:
    var_prod = np.maximum(cov_mats[:, 0, 0] * cov_mats[:, 1, 1], 0.0)
    denom = np.sqrt(var_prod)
    corr = np.divide(
        cov_mats[:, 0, 1],
        denom,
        out=np.zeros_like(cov_mats[:, 0, 1]),
        where=denom > 0,
    )
    finite_cov = np.all(np.isfinite(cov_mats), axis=(1, 2))
    # The parametric fit permits zero variances because A > 0 keeps its
    # marginal likelihood well-defined. NPMLE permits Sigma_11 = 0 but
    # requires Sigma_22 > 0 because it inverts the bias-measurement variance.
    roundoff = 32 * np.finfo(float).eps * np.maximum(
        1.0, np.max(np.abs(cov_mats), axis=(1, 2)),
    )
    valid_variances = (
        (cov_mats[:, 0, 0] >= -roundoff) & (cov_mats[:, 1, 1] >= -roundoff)
        if allow_zero_variance else
        (cov_mats[:, 0, 0] > 0) & (cov_mats[:, 1, 1] > 0)
    )
    valid = (
        np.isfinite(point)
        & np.isfinite(bias)
        & finite_cov
        & valid_variances
    )
    return valid & (corr < corr_threshold)


# ---------------------------------------------------------------------------
# Parametric: symmetric Wald with shrunk variance
# ---------------------------------------------------------------------------

def _param_ci_from_components(
    point: np.ndarray,
    bias: np.ndarray,
    cov_mats: np.ndarray,
    corr_threshold: float,
    alpha: float,
    alternative: str,
) -> np.ndarray:
    """Rebiased estimator + Wald CI; shape ``(M, 2)``.

    For selected problems, we replace the point estimate by
    ``point - Sigma12 / (Sigma22 + A) * (bias - mu)`` and the variance
    by ``Sigma11 - Sigma12^2 / (Sigma22 + A)`` (this is the
    posterior-conditional variance after integrating ``b`` against
    ``N(mu, A)``).
    """
    point = np.asarray(point, dtype=float).reshape(-1)
    bias = np.asarray(bias, dtype=float).reshape(-1)
    cov_mats = np.asarray(cov_mats, dtype=float)

    use_eb = _select_eb_subset(
        point, bias, cov_mats, corr_threshold, allow_zero_variance=True,
    )
    new_point = point.copy()
    new_var = cov_mats[:, 0, 0].copy()

    if use_eb.sum() >= 2:
        mu, A = fit_parametric_bias_prior(bias[use_eb], cov_mats[use_eb])
        sigma_11 = cov_mats[use_eb, 0, 0]
        sigma_12 = cov_mats[use_eb, 0, 1]
        sigma_22 = cov_mats[use_eb, 1, 1]
        denom = sigma_22 + A
        new_point[use_eb] = point[use_eb] - sigma_12 / denom * (bias[use_eb] - mu)
        new_var[use_eb] = sigma_11 - sigma_12 ** 2 / denom

    new_var = np.maximum(new_var, 0.0)
    return _zconfint(new_point, np.sqrt(new_var), alpha, alternative)


def get_eb_ppi_cis(
    data,
    alpha: float = 0.1,
    alternative: str = "two-sided",
    corr_threshold: float = 1.0,
) -> np.ndarray:
    """Rebiased PPI (lambda=1) CI with Gaussian prior on bias."""
    point, bias, cov = ppi_decomposition(
        data.pred_unlabelled, data.pred_labelled, data.y_labelled, 1.0,
    )
    return _param_ci_from_components(point, bias, cov, corr_threshold, alpha, alternative)


def get_eb_pt_cis(
    data,
    alpha: float = 0.1,
    alternative: str = "two-sided",
    corr_threshold: float = 1.0,
    clip_lambda: Tuple[float, float] = (0.0, 1.0),
) -> np.ndarray:
    """Rebiased power-tuned PPI CI with Gaussian prior on bias."""
    lams = pt_lambdas(
        data.pred_unlabelled, data.pred_labelled, data.y_labelled, clip=clip_lambda,
    )
    _, bias, pt, cov = pt_decomposition(
        data.pred_unlabelled, data.pred_labelled, data.y_labelled, lams,
    )
    return _param_ci_from_components(pt, bias, cov, corr_threshold, alpha, alternative)


# ---------------------------------------------------------------------------
# NPMLE: equal-tail CI from a posterior-weighted normal mixture
# ---------------------------------------------------------------------------

def _weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    weights = weights / weights.sum()
    order = np.argsort(values)
    cdf = np.cumsum(weights[order])
    return float(values[order][np.searchsorted(cdf, q, side="left")])


def _mixture_normal_quantile(
    q: float, means: np.ndarray, weights: np.ndarray, sd: float,
) -> float:
    """``q``-quantile of ``sum_k w_k * N(mu_k, sd^2)``."""
    weights = weights / weights.sum()
    if not (0.0 < q < 1.0):
        raise ValueError("q must lie strictly between 0 and 1.")
    if sd <= 1e-12:
        return _weighted_quantile(means, weights, q)

    def cdf(x):
        return float(np.sum(weights * norm.cdf((x - means) / sd)))

    lo = float(np.min(means)) - 10.0 * sd
    hi = float(np.max(means)) + 10.0 * sd
    while cdf(lo) > q:
        lo -= 10.0 * sd
    while cdf(hi) < q:
        hi += 10.0 * sd
    return float(brentq(lambda x: cdf(x) - q, lo, hi))


def _npmle_ci_from_components(
    point: np.ndarray,
    bias: np.ndarray,
    cov_mats: np.ndarray,
    corr_threshold: float,
    alpha: float,
    alternative: str,
    *,
    u: Union[int, np.ndarray],
    threshold: float,
    eps: float,
    max_iter_em: int,
    verbose: bool,
) -> np.ndarray:
    """For each selected problem, the conditional distribution
    ``theta | point, bias`` is a mixture of normals indexed by the
    NPMLE atoms (with the data-update posterior weights). We invert
    the mixture CDF to get equal-tail quantiles, which form the CI.
    Problems failing the selection get the unadjusted Wald interval.
    """
    if not _HAS_NPMLE:
        raise ImportError(
            "NPMLE CIs require the [npmle] extra. "
            "Run `uv sync --extra npmle` and set MOSEKLM_LICENSE_FILE."
        )

    point = np.asarray(point, dtype=float).reshape(-1)
    bias = np.asarray(bias, dtype=float).reshape(-1)
    cov_mats = np.asarray(cov_mats, dtype=float)

    base_var = np.maximum(cov_mats[:, 0, 0], 0.0)
    out = _zconfint(point, np.sqrt(base_var), alpha, alternative)

    use_eb = _select_eb_subset(point, bias, cov_mats, corr_threshold, allow_zero_variance=True)
    use_eb &= (np.isfinite(cov_mats[:, 1, 1])& (cov_mats[:, 1, 1] > 0))
    if use_eb.sum() < 2:
        if verbose:
            print("Skipping NPMLE EB: < 2 selected problems.")
        return out

    atoms, weights = fit_npmle_bias_prior(
        bias[use_eb], cov_mats[use_eb],
        u=u, threshold=threshold, eps=eps,
        max_iter_em=max_iter_em, verbose=verbose,
    )
    post_prob = posterior_weights_npmle(
        bias[use_eb], cov_mats[use_eb], atoms, weights,
    )

    sigma_11 = cov_mats[use_eb, 0, 0]
    sigma_12 = cov_mats[use_eb, 0, 1]
    sigma_22 = cov_mats[use_eb, 1, 1]
    cond_var = np.maximum(sigma_11 - sigma_12 ** 2 / sigma_22, 0.0)

    selected = np.where(use_eb)[0]
    for local, global_idx in enumerate(selected):
        means = (
            point[global_idx]
            - (sigma_12[local] / sigma_22[local])
            * (bias[global_idx] - atoms)
        )
        ws = post_prob[local]
        sd = float(np.sqrt(cond_var[local]))

        if alternative == "two-sided":
            out[global_idx, 0] = _mixture_normal_quantile(alpha / 2.0, means, ws, sd)
            out[global_idx, 1] = _mixture_normal_quantile(1.0 - alpha / 2.0, means, ws, sd)
        elif alternative == "larger":
            out[global_idx, 0] = _mixture_normal_quantile(alpha, means, ws, sd)
            out[global_idx, 1] = np.inf
        elif alternative == "smaller":
            out[global_idx, 0] = -np.inf
            out[global_idx, 1] = _mixture_normal_quantile(1.0 - alpha, means, ws, sd)
        else:
            raise ValueError(
                "alternative must be 'two-sided', 'larger', or 'smaller'; "
                f"got {alternative!r}"
            )

    return out


def get_npmle_eb_ppi_cis(
    data,
    alpha: float = 0.1,
    alternative: str = "two-sided",
    corr_threshold: float = 1.0,
    u: Union[int, np.ndarray] = 30,
    threshold: float = 1e-6,
    eps: float = 1e-5,
    max_iter_em: int = 0,
    verbose: bool = False,
) -> np.ndarray:
    """Rebiased PPI (lambda=1) CI with NPMLE prior on bias."""
    point, bias, cov = ppi_decomposition(
        data.pred_unlabelled, data.pred_labelled, data.y_labelled, 1.0,
    )
    return _npmle_ci_from_components(
        point, bias, cov, corr_threshold, alpha, alternative,
        u=u, threshold=threshold, eps=eps,
        max_iter_em=max_iter_em, verbose=verbose,
    )


def get_npmle_eb_pt_cis(
    data,
    alpha: float = 0.1,
    alternative: str = "two-sided",
    corr_threshold: float = 1.0,
    clip_lambda: Tuple[float, float] = (0.0, 1.0),
    u: Union[int, np.ndarray] = 30,
    threshold: float = 1e-6,
    eps: float = 1e-5,
    max_iter_em: int = 0,
    verbose: bool = False,
) -> np.ndarray:
    """Rebiased power-tuned PPI CI with NPMLE prior on bias."""
    lams = pt_lambdas(
        data.pred_unlabelled, data.pred_labelled, data.y_labelled, clip=clip_lambda,
    )
    _, bias, pt, cov = pt_decomposition(
        data.pred_unlabelled, data.pred_labelled, data.y_labelled, lams,
    )
    return _npmle_ci_from_components(
        pt, bias, cov, corr_threshold, alpha, alternative,
        u=u, threshold=threshold, eps=eps,
        max_iter_em=max_iter_em, verbose=verbose,
    )
