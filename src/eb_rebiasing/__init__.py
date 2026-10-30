"""eb-rebiasing: empirical Bayes rebiasing for prediction-powered CIs.

Importing this package automatically registers its rebiased CI methods
into ``pas.intervals.CORE_CI_METHODS``, so scripts that consume the
upstream registry pick them up without further wiring.
"""
from eb_rebiasing.intervals import (
    get_eb_ppi_cis,
    get_eb_pt_cis,
)
from eb_rebiasing.prior_fit import _HAS_NPMLE
from eb_rebiasing.registry import REBIASED_CI_METHODS, register_with_pas
from eb_rebiasing.plotting import plot_bias_scatter, plot_bias_vs_nuisance

if _HAS_NPMLE:
    from eb_rebiasing.intervals import (
        get_npmle_eb_ppi_cis,
        get_npmle_eb_pt_cis,
    )

# Register on import.
register_with_pas()


__all__ = [
    "get_eb_ppi_cis",
    "get_eb_pt_cis",
    "REBIASED_CI_METHODS",
    "register_with_pas",
    "_HAS_NPMLE",
    "plot_bias_scatter",
    "plot_bias_vs_nuisance",
]

if _HAS_NPMLE:
    __all__ += [
        "get_npmle_eb_ppi_cis",
        "get_npmle_eb_pt_cis",
    ]
