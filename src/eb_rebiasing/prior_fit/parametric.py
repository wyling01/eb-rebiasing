"""Parametric Gaussian prior on the bias.

Given observed bias estimates ``b_hat_i`` with measurement variances
``Var(b_hat_i) = cov_mats[i, 1, 1]`` and the model
``b_hat_i | b_i ~ N(b_i, Var(b_hat_i))``, ``b_i ~ N(mu, A)``, the
marginal MLE of ``(mu, A)`` solves a small one-dimensional optimization
in ``log A`` (with ``mu`` profiled or jointly optimized). We use BFGS
on ``(mu, log A)``.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize


def fit_parametric_bias_prior(
    bias: np.ndarray,
    cov_mats: np.ndarray,
) -> tuple[float, float]:
    """Marginal MLE for the Gaussian prior ``b ~ N(mu, A)``.

    Args:
        bias: shape ``(M,)``. Observed unbiased estimates ``b_hat_i``.
        cov_mats: shape ``(M, 2, 2)``. The ``[1, 1]`` slot per problem
            is interpreted as ``Var(b_hat_i)``.

    Returns:
        ``(mu_hat, A_hat)``. ``A_hat`` is bounded below by ``1e-12`` so
        downstream Bayes-rule computations never divide by zero.
    """
    bias = np.asarray(bias, dtype=float).reshape(-1)
    cov_mats = np.asarray(cov_mats, dtype=float)
    meas_var = cov_mats[:, 1, 1]

    init_mean = float(bias.mean())
    s2_obs = float(np.var(bias, ddof=1)) if len(bias) > 1 else 0.0
    init_var = max(s2_obs - float(meas_var.mean()), 1e-4)

    def neg_log_lik(par):
        mu_b, log_A = par
        A = np.exp(log_A)
        total_var = A + meas_var
        resid2 = (bias - mu_b) ** 2
        return 0.5 * np.sum(np.log(2 * np.pi * total_var) + resid2 / total_var)

    res = minimize(
        neg_log_lik,
        x0=np.array([init_mean, np.log(init_var)]),
        method="BFGS",
    )
    mu_hat = float(res.x[0])
    A_hat = max(float(np.exp(res.x[1])), 1e-12)
    return mu_hat, A_hat
