"""LMArena CI benchmark — run experiment and save raw per-trial results.

Saves three files in `tmp/lmarena_ci_results/` (so the outputs are not
tracked by git):
    raw.csv      -- one row per (alpha, trial), columns are
                    {method}_coverage and {method}_width for every method
    summary.csv  -- mean/SE per (alpha, method)
    meta.json    -- run config (trials, num_workers, methods, etc.)

The plot and the appendix table are produced by separate scripts that
load `raw.csv`, so plot/table iteration does not require re-running the
experiment.

Usage (from repo root):

    MOSEKLM_LICENSE_FILE=/path/to/mosek.lic \
        uv run python scripts/run_lmarena_ci.py --trials 200 --num-workers 4

The PPI/PT baselines, the two robust empirical-Bayes competitors (robust
EBCI and double shrinkage), and the LMArena dataset come from the upstream
``pas`` package. The rebiased estimators (PT-Param, PT-NPMLE) come from
``eb_rebiasing`` and auto-register on import.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=RuntimeWarning, module="numpy")
warnings.filterwarnings("ignore", category=DeprecationWarning)

import eb_rebiasing  # noqa: F401  (auto-registers methods into pas registry)
from eb_rebiasing.data_paths import LMARENA_CSV, LMARENA_RESULTS
from pas.datasets.lmarena import LMArenaDataset
from pas.experiments import run_ci_benchmark
from pas.intervals import CORE_CI_METHODS


# Shared robust-EBCI configuration. `shrink_to="grand_mean"` is the R formula
# `Y ~ 1`: each interval is centred on the base estimator shrunk toward the
# inverse-variance-weighted grand mean of the base estimates across tasks (the
# alternative, `shrink_to="zero"`, shrinks toward 0 instead).
#
# `cv_mode="exact"` solves the robust critical value per task instead of
# reading it off the packaged lookup table. The tables round m2 (and kappa)
# *up*, which would widen this competitor baseline by ~1% on average -- a
# systematic bias in our own favour. Exact mode costs ~12s per replicate,
# small next to the NPMLE method.
_ROBUST_EB_COMMON = {
    "share_var": False,
    "weights": "inverse_variance",
    "shrink_to": "grand_mean",
    "fs_correction": "PMT",
    "cv_mode": "exact",
}

# Run key -> (key in `pas.intervals.CORE_CI_METHODS`, extra kwargs).
#
# The run key is what appears in raw.csv, so the two robust-EBCI variants --
# the Armstrong et al. interval applied to the power-tuned estimator and to the
# classical one -- get distinct names while sharing one implementation.
#
# The PT-Param / PT-NPMLE rows are this paper's contribution; the rest are
# baselines from the upstream `pas` package.
METHOD_SPECS = {
    "mle_ci":                 ("mle_ci", {}),                    # Classical (CLT)
    "pred_mean_ci":           ("pred_mean_ci", {}),              # Prediction Mean
    "ppi_ci":                 ("ppi_ci", {}),                    # Vanilla PPI
    "pt_ci":                  ("pt_ci", {"share_var": False}),   # Power-Tuned PPI
    "robust_eb_pt_ci":        ("robust_eb_ci",
                               {"base_estimator": "pt", **_ROBUST_EB_COMMON}),
    "robust_eb_classical_ci": ("robust_eb_ci",
                               {"base_estimator": "classical", **_ROBUST_EB_COMMON}),
    "eb_pt_ci":               ("eb_pt_ci", {}),                  # PT-Param  (ours)
    "eb_npmle_pt_ci":         ("eb_npmle_pt_ci", {"u": 50}),            # PT-NPMLE  (ours)
    # Double shrinkage (Rosenman et al. 2023) was dropped from the headline
    # comparison in favour of the two robust-EBCI variants above. Re-enable by
    # uncommenting; the method is still registered upstream.
    # "double_shrinkage_mle_ci": ("double_shrinkage_mle_ci",
    #                             {"center": 0.5,
    #                              "finite_sample_correction": "pmt",
    #                              "cv_mode": "exact"}),
    # Other rebiased flavours / double-shrinkage hyperparameter variants kept
    # for ablations; uncomment to enable:
    # "double_shrinkage_mm1_ci": ("double_shrinkage_mm1_ci", {}),
    # "double_shrinkage_mm2_ci": ("double_shrinkage_mm2_ci", {}),
    # "double_shrinkage_ure_ci": ("double_shrinkage_ure_ci", {}),
    # "eb_ppi_ci":              ("eb_ppi_ci", {}),
    # "eb_npmle_ppi_ci":        ("eb_npmle_ppi_ci", {}),
}

RUN_METHODS = list(METHOD_SPECS)


def _register_aliases() -> None:
    """Expose each run key in ``CORE_CI_METHODS``.

    The parallel path in ``run_ci_benchmark`` ships only method *names* to its
    workers and re-resolves them against ``CORE_CI_METHODS`` there, so a run
    key that is not in the registry raises ``KeyError`` under
    ``--num-workers > 1``. Registering at import time covers the workers too:
    the spawn context re-imports this module in each child (which is also how
    the ``eb_rebiasing`` methods get registered).
    """
    for run_key, (registry_key, _) in METHOD_SPECS.items():
        if run_key == registry_key:
            continue
        if registry_key in CORE_CI_METHODS:
            CORE_CI_METHODS[run_key] = CORE_CI_METHODS[registry_key]


_register_aliases()


def _summarize(raw: pd.DataFrame, methods, reference_method="mle_ci") -> pd.DataFrame:
    records = []
    for alpha, group in raw.groupby("alpha"):
        reference_width = group[f"{reference_method}_width"].astype(float)
        if (reference_width <= 0).any():
            raise ValueError(
                f"{reference_method} contains nonpositive trial-average widths."
            )
        for method in methods:
            cov = group[f"{method}_coverage"].astype(float)
            wid = group[f"{method}_width"].astype(float)
            wid_ratio = wid / reference_width
            records.append({
                "alpha": float(alpha),
                "method": method,
                "mean_coverage": cov.mean(),
                "se_coverage": cov.sem(),
                "mean_width": wid.mean(),
                "se_width": wid.sem(),
                "mean_width_ratio": wid_ratio.mean(),
                "se_width_ratio": wid_ratio.sem(),
            })
    return pd.DataFrame(records)


def _git_revision(path):
    result = subprocess.run(
        ["git", "-C", str(path), "rev-parse", "HEAD"],
        capture_output=True, text=True,
    )
    return result.stdout.strip() if result.returncode == 0 else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=200)
    parser.add_argument("--alphas", type=float, nargs="+",
                        default=[0.05, 0.10, 0.20, 0.30])
    parser.add_argument("--train-test-split", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=42,
                        help="First split seed; trial i uses seed + i.")
    parser.add_argument("--num-workers", type=int, default=4)
    parser.add_argument("--out-dir", type=Path,
                        default=LMARENA_RESULTS)
    parser.add_argument("--csv", type=Path, default=LMARENA_CSV,
                        help="LMArena cleaned CSV; defaults to Bradley-Terry probabilities.")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    dataset_kwargs = {"file_path": args.csv.resolve(),
                      "train_test_split": args.train_test_split,
                      "split_seed": args.seed}
    dataset = LMArenaDataset(**dataset_kwargs)
    print(f"LMArena source: {dataset.file_path}")
    print(f"LMArena: M={dataset.M} problems, "
          f"n=[{dataset.ns.min()},{dataset.ns.max()}], "
          f"N=[{dataset.Ns.min()},{dataset.Ns.max()}]")

    missing = [k for k in RUN_METHODS if k not in CORE_CI_METHODS]
    if missing:
        raise RuntimeError(
            f"Missing CI methods (install npmle extra?): {missing}. "
            "Run `uv sync --extra npmle` and set MOSEKLM_LICENSE_FILE."
        )
    # Both robust-EBCI keys resolve to the same upstream callable and are told
    # apart by their kwargs. Matches
    # ebci(Y ~ 1, se=se, weights=1/se^2, fs_correction="PMT") in the reference
    # R package, with Y the chosen base estimator.
    ci_methods = {k: CORE_CI_METHODS[k] for k in RUN_METHODS}
    ci_kwargs = {
        run_key: kwargs
        for run_key, (_, kwargs) in METHOD_SPECS.items()
        if kwargs
    }

    import pas
    repo_root = Path(__file__).resolve().parent.parent
    meta = {
        "csv": str(dataset.file_path.resolve()),
        "csv_sha256": hashlib.sha256(dataset.file_path.read_bytes()).hexdigest(),
        "trials": args.trials,
        "alphas": list(args.alphas),
        "train_test_split": args.train_test_split,
        "split_seed": args.seed,
        "num_problems": dataset.M,
        "num_workers": args.num_workers,
        "methods": RUN_METHODS,
        "ci_kwargs": ci_kwargs,
        "python": sys.version,
        "repo_revision": _git_revision(repo_root),
        "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "pas_revision": _git_revision(Path(pas.__file__).parent),
        "packages": {dist.metadata["Name"]: dist.version
                     for dist in importlib.metadata.distributions()},
        "completed_alphas": [],
        "complete": False,
    }
    meta_path = args.out_dir / "meta.json"
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")

    raw_frames = []
    for alpha in args.alphas:
        print(f"\n=== alpha = {alpha} ===")
        ci_results = run_ci_benchmark(
            dataset,
            trials=args.trials,
            alpha=alpha,
            summary=False,
            ci_methods=ci_methods,
            ci_kwargs=ci_kwargs,
            num_workers=args.num_workers,
            verbose=True,
        )
        ci_results = ci_results.reset_index().rename(columns={"index": "trial"})
        ci_results.insert(0, "alpha", float(alpha))
        # Preserve completed settings if a later solve is interrupted.
        ci_results.to_csv(args.out_dir / f"raw_alpha_{alpha:g}.csv", index=False)
        raw_frames.append(ci_results)
        meta["completed_alphas"].append(float(alpha))
        meta_path.write_text(json.dumps(meta, indent=2) + "\n")

    raw = pd.concat(raw_frames, ignore_index=True)
    raw.to_csv(args.out_dir / "raw.csv", index=False)
    print(f"\nSaved raw -> {args.out_dir / 'raw.csv'} ({len(raw)} rows)")

    summary = _summarize(raw, RUN_METHODS)
    summary.to_csv(args.out_dir / "summary.csv", index=False)
    print(f"Saved summary -> {args.out_dir / 'summary.csv'}")
    print(summary.round(4).to_string(index=False))

    meta["complete"] = True
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")
    print(f"Saved meta -> {args.out_dir / 'meta.json'}")


if __name__ == "__main__":
    main()
