"""Prior-on-bias fitters used by the rebiased estimators / CIs."""
from eb_rebiasing.prior_fit.parametric import fit_parametric_bias_prior

try:
    from eb_rebiasing.prior_fit.npmle import (
        fit_npmle_bias_prior,
        posterior_weights_npmle,
    )
    _HAS_NPMLE = True
except ImportError:
    _HAS_NPMLE = False

__all__ = ["fit_parametric_bias_prior", "_HAS_NPMLE"]

if _HAS_NPMLE:
    __all__ += ["fit_npmle_bias_prior", "posterior_weights_npmle"]
