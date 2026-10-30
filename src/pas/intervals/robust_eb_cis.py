"""Robust empirical-Bayes confidence intervals for unbiased PAS estimators.

The intervals implement the default ``ebci::ebci`` construction from
Armstrong, Kolesar, and Plagborg-Moller (2022).  Given the fitted rule

    theta_hat_i = mu_i + w_i (Y_i - mu_i),
    w_i = mu2 / (mu2 + se_i^2),

its robust interval is

    theta_hat_i +/- cva_alpha(se_i^2 / mu2, kappa) * w_i * se_i.

The default PAS configuration uses PT as the unbiased preliminary estimator,
shrinks toward its inverse-variance-weighted grand mean, fits the moments with
weights ``1 / se_i^2``, and uses ``fs_correction='PMT'``.
"""
from __future__ import annotations

from typing import Literal, Optional, Sequence, Union
import warnings

import numpy as np
from scipy.stats import norm

from pas.estimators.robust_eb import (
    BaseEstimator,
    FiniteSampleCorrection,
    fit_robust_eb,
    get_unbiased_estimates_and_ses,
)
from pas.intervals.robust_eb import robust_eb_critical_values


def get_robust_eb_cis(
    data,
    alpha: float = 0.1,
    alternative: str = "two-sided",
    *,
    base_estimator: BaseEstimator = "pt",
    share_var: bool = True,
    weights: Union[
        Literal["inverse_variance", "equal"],
        Sequence[float],
        np.ndarray,
    ] = "inverse_variance",
    design: Optional[np.ndarray] = None,
    shrink_to: Literal["grand_mean", "zero"] = "grand_mean",
    fs_correction: FiniteSampleCorrection = "PMT",
    kappa: Optional[float] = None,
    cv_mode: Literal["lookup", "exact"] = "exact",
) -> np.ndarray:
    """Construct robust EBCIs by shrinking an unbiased estimator only.

    ``alternative`` is restricted to ``'two-sided'`` because the source method
    derives symmetric two-sided average-coverage intervals.
    """
    if alternative != "two-sided":
        raise ValueError("robust EBCIs support only alternative='two-sided'")
    if not (0.0 < alpha < 1.0):
        raise ValueError("alpha must lie in (0, 1)")

    estimates, standard_errors = get_unbiased_estimates_and_ses(data,
                                                                base_estimator=base_estimator,
                                                                share_var=share_var,)
    estimates = np.asarray(estimates, dtype=float)
    standard_errors = np.asarray(standard_errors, dtype=float)

    if np.any(~np.isfinite(estimates)):
        raise ValueError("Estimates must be finite.")

    if np.any(~np.isfinite(standard_errors)) or np.any(standard_errors < 0):
        raise ValueError("Standard errors must be finite and non-negative.")

    positive_se = standard_errors > 0

    wald_half_length = norm.ppf(1.0 - alpha / 2.0) * standard_errors
    out = np.column_stack([estimates - wald_half_length,
                           estimates + wald_half_length,])
    if isinstance(weights, str):
        if weights == "inverse_variance":
            fit_weights = 1.0 / standard_errors[positive_se] ** 2
        elif weights == "equal":
            fit_weights = np.ones(positive_se.sum(), dtype=float)
        else:
            raise ValueError(
                "weights must be 'inverse_variance', 'equal', or an array"
            )
    else:
        weight_array = np.asarray(weights, dtype=float).reshape(-1)
        if weight_array.size != data.M:
            raise ValueError("A supplied weight array must have one value per task.")
        fit_weights = weight_array[positive_se]

    if design is None: fit_design = None
    else:
        design_array = np.asarray(design, dtype=float)
        if design_array.shape[0] != data.M:
            raise ValueError("design must have one row per task.")
        fit_design = design_array[positive_se]

    fit = fit_robust_eb(
        estimates[positive_se],
        standard_errors[positive_se],
        weights=fit_weights,
        design=fit_design,
        shrink_to=shrink_to,
        fs_correction=fs_correction,
        kappa=kappa,
    )
    if fit.mu2 <= 0.0:
        warnings.warn(
            "The robust-EB second-moment estimate is zero; returning the "
            "unshrunk Wald intervals. Use fs_correction='PMT' to avoid this "
            "finite-sample collapse.",
            RuntimeWarning,
            stacklevel=2,
        )
        return out

    normalized_bias_m2 = standard_errors[positive_se] ** 2 / fit.mu2
    critical_values = robust_eb_critical_values(
        normalized_bias_m2,
        alpha=alpha,
        kappa=fit.kappa,
        mode=cv_mode,
    )
    eb_half_length = critical_values * fit.shrinkage * fit.standard_error
    out[positive_se, 0] = fit.estimate - eb_half_length
    out[positive_se, 1] = fit.estimate + eb_half_length

    return out


__all__ = ["get_robust_eb_cis"]
