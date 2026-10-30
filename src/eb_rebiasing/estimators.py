"""Rebiased point estimators.

For each task ``i`` we observe a (noisy) biased estimator
``theta_b_i`` and a (noisier) unbiased estimate of its bias ``b_hat_i``.
Modeling ``b_i ~ G`` with ``G`` either Gaussian (parametric) or a
discrete mixture (NPMLE), we replace the full debias correction
``theta_b_i - b_hat_i`` with a partially-Bayes correction
``theta_b_i - E_G[b_i | b_hat_i]``.

Each public estimator returns an ``(M,)`` array. The corresponding
confidence intervals live in ``eb_rebiasing.intervals``.
"""
from __future__ import annotations

from typing import Tuple

import numpy as np

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
# Internal helpers (apply the rebias correction once we have (point, bias, cov))
# ---------------------------------------------------------------------------

def _select_eb_subset(
    point: np.ndarray, bias: np.ndarray, cov_mats: np.ndarray, corr_threshold: float,
    *, allow_zero_variance: bool = False,
) -> np.ndarray:
    """Eligible tasks; Gaussian priors also admit zero measurement variance."""
    point = np.asarray(point, dtype=float)
    bias = np.asarray(bias, dtype=float)
    cov_mats = np.asarray(cov_mats, dtype=float)

    var_prod = np.maximum(cov_mats[:, 0, 0] * cov_mats[:, 1, 1], 0.0)
    denom = np.sqrt(var_prod)
    corr = np.divide(
        cov_mats[:, 0, 1],
        denom,
        out=np.zeros_like(cov_mats[:, 0, 1]),
        where=denom > 0,
    )
    finite_cov = np.all(np.isfinite(cov_mats), axis=(1, 2))
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


def _apply_param_rebias(
    point: np.ndarray, bias: np.ndarray, cov_mats: np.ndarray, corr_threshold: float
) -> np.ndarray:
    """``point - (Sigma12 / (Sigma22 + A)) * (bias - mu)`` on selected
    problems; return ``point`` unchanged elsewhere."""
    out = np.asarray(point, dtype=float).copy()
    use_eb = _select_eb_subset(
        point, bias, cov_mats, corr_threshold, allow_zero_variance=True,
    )
    if use_eb.sum() < 2:
        return out
    mu, A = fit_parametric_bias_prior(bias[use_eb], cov_mats[use_eb])
    sigma_12 = cov_mats[use_eb, 0, 1]
    sigma_22 = cov_mats[use_eb, 1, 1]
    out[use_eb] = point[use_eb] - sigma_12 / (sigma_22 + A) * (bias[use_eb] - mu)
    return out


def _apply_npmle_rebias(
    point: np.ndarray,
    bias: np.ndarray,
    cov_mats: np.ndarray,
    corr_threshold: float,
    *,
    u, threshold, eps, max_iter_em, verbose,
) -> np.ndarray:
    """``point - (Sigma12 / Sigma22) * (bias - E_G[b | bias])`` on
    selected problems, where ``E_G`` uses the NPMLE prior."""
    if not _HAS_NPMLE:
        raise ImportError(
            "NPMLE rebiasing requires the [npmle] extra. "
            "Run `uv sync --extra npmle` and set MOSEKLM_LICENSE_FILE."
        )
    out = np.asarray(point, dtype=float).copy()
    use_eb = _select_eb_subset(point, bias, cov_mats, corr_threshold,
                               allow_zero_variance=True)
    use_eb &= (np.isfinite(cov_mats[:, 1, 1])& (cov_mats[:, 1, 1] > 0))
    if use_eb.sum() < 2:
        return out

    atoms, weights = fit_npmle_bias_prior(
        bias[use_eb], cov_mats[use_eb],
        u=u, threshold=threshold, eps=eps,
        max_iter_em=max_iter_em, verbose=verbose,
    )
    post_prob = posterior_weights_npmle(
        bias[use_eb], cov_mats[use_eb], atoms, weights
    )
    post_mean = post_prob @ atoms

    sigma_12 = cov_mats[use_eb, 0, 1]
    sigma_22 = cov_mats[use_eb, 1, 1]
    out[use_eb] = point[use_eb] - (sigma_12 / sigma_22) * (bias[use_eb] - post_mean)
    return out


# ---------------------------------------------------------------------------
# Parametric rebiased point estimators
# ---------------------------------------------------------------------------

def get_param_rebiased_ppi_estimators(data, corr_threshold: float = 1.0) -> np.ndarray:
    """Rebiased PPI (lambda=1) point estimator with Gaussian bias prior."""
    point, bias, cov = ppi_decomposition(
        data.pred_unlabelled, data.pred_labelled, data.y_labelled, 1.0,
    )
    return _apply_param_rebias(point, bias, cov, corr_threshold)


def get_param_rebiased_pt_estimators(
    data,
    corr_threshold: float = 1.0,
    clip_lambda: Tuple[float, float] = (0.0, 1.0),
) -> np.ndarray:
    """Rebiased power-tuned PPI (problem-wise lambda) with Gaussian bias prior."""
    lams = pt_lambdas(
        data.pred_unlabelled, data.pred_labelled, data.y_labelled, clip=clip_lambda,
    )
    _, bias, pt, cov = pt_decomposition(
        data.pred_unlabelled, data.pred_labelled, data.y_labelled, lams,
    )
    return _apply_param_rebias(pt, bias, cov, corr_threshold)


# ---------------------------------------------------------------------------
# NPMLE rebiased point estimators
# ---------------------------------------------------------------------------

def get_npmle_rebiased_ppi_estimators(
    data,
    corr_threshold: float = 1.0,
    u=30, threshold: float = 1e-6, eps: float = 1e-5,
    max_iter_em: int = 0, verbose: bool = False,
) -> np.ndarray:
    """Rebiased PPI (lambda=1) point estimator with NPMLE bias prior."""
    point, bias, cov = ppi_decomposition(
        data.pred_unlabelled, data.pred_labelled, data.y_labelled, 1.0,
    )
    return _apply_npmle_rebias(
        point, bias, cov, corr_threshold,
        u=u, threshold=threshold, eps=eps,
        max_iter_em=max_iter_em, verbose=verbose,
    )


def get_npmle_rebiased_pt_estimators(
    data,
    corr_threshold: float = 1.0,
    clip_lambda: Tuple[float, float] = (0.0, 1.0),
    u=30, threshold: float = 1e-6, eps: float = 1e-5,
    max_iter_em: int = 0, verbose: bool = False,
) -> np.ndarray:
    """Rebiased power-tuned PPI (problem-wise lambda) with NPMLE bias prior."""
    lams = pt_lambdas(
        data.pred_unlabelled, data.pred_labelled, data.y_labelled, clip=clip_lambda,
    )
    _, bias, pt, cov = pt_decomposition(
        data.pred_unlabelled, data.pred_labelled, data.y_labelled, lams,
    )
    return _apply_npmle_rebias(
        pt, bias, cov, corr_threshold,
        u=u, threshold=threshold, eps=eps,
        max_iter_em=max_iter_em, verbose=verbose,
    )
