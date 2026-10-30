"""Confidence-interval registry needed by the LMArena benchmark."""

from pas.intervals.ppi_cis import get_pt_ppi_cis, get_vanilla_ppi_cis
from pas.intervals.robust_eb_cis import get_robust_eb_cis
from pas.intervals.simple_cis import get_mle_cis, get_pred_mean_cis

CORE_CI_METHODS = {
    "mle_ci": get_mle_cis,
    "pred_mean_ci": get_pred_mean_cis,
    "ppi_ci": get_vanilla_ppi_cis,
    "pt_ci": get_pt_ppi_cis,
    "robust_eb_ci": get_robust_eb_cis,
}

__all__ = [
    "get_mle_cis",
    "get_pred_mean_cis",
    "get_vanilla_ppi_cis",
    "get_pt_ppi_cis",
    "get_robust_eb_cis",
    "CORE_CI_METHODS",
]
