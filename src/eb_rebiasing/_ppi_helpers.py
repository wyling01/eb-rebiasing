"""Single import surface for the PPI / power-tuned PPI decomposition
helpers we re-use from the upstream `pas` package.

The rest of the `eb_rebiasing` package only reaches into upstream
internals through this module, so renames or refactors upstream are
contained to one file.

Each helper returns a ``(point_estimate, bias, cov_matrix)`` triple
(or, for the power-tuned variant, additionally the per-problem biased
estimator); the rebiased estimators / CIs in this package consume that
shape.
"""
from pas.estimators.eb_estimators import (
    _get_generic_ppi_estimators_plus_bias as ppi_decomposition,
    _get_power_tuned_lambdas as pt_lambdas,
    _get_power_tuned_ppi_estimators_plus_bias as pt_decomposition,
)

__all__ = [
    "ppi_decomposition",
    "pt_lambdas",
    "pt_decomposition",
]
