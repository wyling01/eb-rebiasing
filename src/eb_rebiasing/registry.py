"""Registry of CI methods contributed by this package.

Mutates ``pas.intervals.CORE_CI_METHODS`` on registration so any
script driving the upstream experiment runner sees our methods
without explicit wiring. Idempotent: re-importing this package does
not double-register.
"""
from __future__ import annotations

from pas.intervals import CORE_CI_METHODS as _PAS_CI

from eb_rebiasing.intervals import (
    get_eb_ppi_cis,
    get_eb_pt_cis,
)
from eb_rebiasing.prior_fit import _HAS_NPMLE


REBIASED_CI_METHODS = {
    "eb_ppi_ci":   get_eb_ppi_cis,
    "eb_pt_ci":    get_eb_pt_cis,
}

if _HAS_NPMLE:
    from eb_rebiasing.intervals import (
        get_npmle_eb_ppi_cis,
        get_npmle_eb_pt_cis,
    )
    REBIASED_CI_METHODS.update({
        "eb_npmle_ppi_ci":   get_npmle_eb_ppi_cis,
        "eb_npmle_pt_ci":    get_npmle_eb_pt_cis,
    })


def register_with_pas() -> None:
    """Add this package's CI methods to ``pas.intervals.CORE_CI_METHODS``.

    The rebiased CI implementations in this package are the source of
    truth, so any same-keyed method registered upstream is overwritten.
    """
    for name, fn in REBIASED_CI_METHODS.items():
        _PAS_CI[name] = fn


__all__ = ["REBIASED_CI_METHODS", "register_with_pas"]
