"""Nonparametric MLE for the bias prior.

We discretize the prior support onto an atom grid and solve the
Koenker--Mizera convex program (via the ``npeb`` package, which calls
Mosek under the hood) for the atom weights. Optionally polish the
weights with EM. Returns a discrete prior
``sum_k weights_k * delta_{atoms_k}``.

Importing this module raises ``ImportError`` if the ``[npmle]`` extra
is not installed; the package re-exports a ``_HAS_NPMLE`` flag at the
top level so callers can degrade gracefully.
"""
from __future__ import annotations

import warnings
from typing import Union

import numpy as np
from npeb import GLMixture
from scipy.special import logsumexp

import io
from contextlib import nullcontext, redirect_stdout


def _fit_on_grid(X, prec, u_grid, *, threshold, max_iter_em):
    model = GLMixture(
        prec_type="diagonal",
        homoscedastic=False,
        atoms_init=u_grid.reshape(-1, 1),
    )
    model.fit(
        X,
        prec,
        max_iter_em=max_iter_em,
        weight_thresh=threshold,
        row_condition=True,
        solver="mosek",
    )
    return model.get_params()


def fit_npmle_bias_prior(
    bias: np.ndarray,
    cov_mats: np.ndarray,
    *,
    u: Union[int, np.ndarray] = 30,
    threshold: float = 1e-6,
    eps: float = 1e-5,
    max_iter_em: int = 0,
    verbose: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Fit the NPMLE point-mass prior on the latent bias ``b``.

    Model:
        ``b_hat_i | b_i ~ N(b_i, cov_mats[i, 1, 1])``,
        ``b_i ~ sum_k weights_k * delta_{atoms_k}``.

    Args:
        bias: shape ``(M,)``. Observed bias estimates ``b_hat_i``.
        cov_mats: shape ``(M, 2, 2)``. ``cov_mats[i, 1, 1]`` is
            ``Var(b_hat_i)``.
        u: either the number of atoms (an int >= 2) or the explicit
            grid of atom locations.
        threshold: weight threshold below which atoms are pruned.
        eps: padding added to the data range when ``u`` is an int.
        max_iter_em: number of EM polishing steps inside ``npeb.fit``.
        verbose: pass-through for progress messages.

    Returns:
        ``(atoms, weights)`` of shape ``(K,)`` each, with weights
        normalized to sum to 1.
    """
    bias = np.asarray(bias, dtype=float).reshape(-1)
    cov_mats = np.asarray(cov_mats, dtype=float)

    if np.isscalar(u):
        n_grid = int(u)
        if n_grid < 2:
            raise ValueError("If u is an integer, it must be at least 2.")
        u_grid = np.linspace(bias.min() - eps, bias.max() + eps, n_grid)
    else:
        u_grid = np.asarray(u, dtype=float).reshape(-1)
        if u_grid.size < 2:
            raise ValueError("If u is array-like, it must have at least 2 atoms.")

    var_bias_hat = cov_mats[:, 1, 1]
    if np.any(~np.isfinite(var_bias_hat)) or np.any(var_bias_hat <= 0):
        raise ValueError("cov_mats[:, 1, 1] must be finite and strictly positive.")

    X = bias.reshape(-1, 1)
    prec = (1.0 / var_bias_hat).reshape(-1, 1)

    if verbose:
        print("Solving NPMLE for bias prior using npeb...")

    # Mosek's interior-point solve occasionally returns status "Unknown" on a
    # particular discretization (roughly 1 replicate in 300 in the Amazon
    # benchmark).  The convex program is unchanged by a marginal shift of the
    # atom grid -- the discretization is an approximation choice, not part of
    # the estimator -- so retry on a slightly perturbed grid rather than
    # dropping the replicate, which would preferentially discard the hardest
    # configurations and bias the reported MSE.  On the failures observed, the
    # recovered estimates agree across all of these variants to five
    # significant figures.
    if np.isscalar(u):
        retry_grids = [
            np.linspace(bias.min() - e, bias.max() + e, n)
            for e, n in (
                (eps * 100.0, n_grid),
                (eps * 100.0, max(2, n_grid - 5)),
                (eps * 100.0, n_grid + 10),
            )
        ]
    else:
        retry_grids = []

    last_error: Exception | None = None
    for attempt, grid in enumerate([u_grid, *retry_grids]):
        try:
            output_context = nullcontext() if verbose else redirect_stdout(io.StringIO())
            with output_context:
                atoms, weights = _fit_on_grid(
                    X,
                    prec,
                    grid,
                    threshold=threshold,
                    max_iter_em=max_iter_em,
                )
            if attempt > 0:
                warnings.warn(
                    f"NPMLE Mosek solve failed on the default atom grid; "
                    f"succeeded after {attempt} grid perturbation(s).",
                    RuntimeWarning,
                    stacklevel=2,
                )
            break
        except Exception as error:  # noqa: BLE001 - npeb/Mosek raise varied types
            last_error = error
    else:
        raise RuntimeError(
            "NPMLE solve failed on the default atom grid and on all "
            f"{len(retry_grids)} perturbed grids."
        ) from last_error
    atoms = np.asarray(atoms, dtype=float).reshape(-1)
    weights = np.asarray(weights, dtype=float).reshape(-1)

    keep = weights > threshold
    atoms = atoms[keep]
    weights = weights[keep]
    if weights.sum() <= 0:
        raise RuntimeError("All NPMLE weights were thresholded out. Try lowering threshold.")
    weights = np.maximum(weights, 0.0)
    weights = weights / weights.sum()
    return atoms, weights


def posterior_weights_npmle(
    bias: np.ndarray,
    cov_mats: np.ndarray,
    atoms: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    """Per-problem posterior probabilities ``P(b_i = atom_k | b_hat_i)``.

    Returns:
        Array of shape ``(M, K)``, rows sum to 1.
    """
    bias = np.asarray(bias, dtype=float).reshape(-1)
    cov_mats = np.asarray(cov_mats, dtype=float)
    atoms = np.asarray(atoms, dtype=float).reshape(-1)
    weights = np.asarray(weights, dtype=float).reshape(-1)
    weights = weights / weights.sum()
    var_bias_hat = cov_mats[:, 1, 1]

    log_lik = (
        -0.5 * np.log(2 * np.pi * var_bias_hat[:, None])
        -0.5 * (bias[:, None] - atoms[None, :]) ** 2 / var_bias_hat[:, None]
    )
    log_prior = np.log(np.maximum(weights, 1e-300))[None, :]
    log_joint = log_lik + log_prior
    log_norm = logsumexp(log_joint, axis=1, keepdims=True)
    return np.exp(log_joint - log_norm)
