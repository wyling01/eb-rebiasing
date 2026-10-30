"""Generate the LM Arena prior and bias data used by the plotting notebooks."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from eb_rebiasing._ppi_helpers import pt_decomposition, pt_lambdas
from eb_rebiasing.data_paths import LMARENA_CSV
from eb_rebiasing.prior_fit.npmle import fit_npmle_bias_prior
from eb_rebiasing.prior_fit.parametric import fit_parametric_bias_prior
from pas.datasets.lmarena import LMArenaDataset


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=LMARENA_CSV)
    parser.add_argument("--train-test-split", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--npmle-grid-size", type=int, default=50)
    parser.add_argument("--out-dir", type=Path, default=Path("tmp/lmarena_plotting"))
    args = parser.parse_args()

    if not 0 < args.train_test_split < 1:
        parser.error("--train-test-split must be between 0 and 1.")
    if args.npmle_grid_size < 2:
        parser.error("--npmle-grid-size must be at least 2.")

    data = LMArenaDataset(file_path=args.csv.resolve(), train_test_split=args.train_test_split, split_seed=args.seed)
    lambdas = pt_lambdas(data.pred_unlabelled, data.pred_labelled, data.y_labelled)
    _, bias, _, cov_mats = pt_decomposition(data.pred_unlabelled, data.pred_labelled, data.y_labelled, lambdas)
    bias_sd = np.sqrt(cov_mats[:, 1, 1])
    if not np.isfinite(np.column_stack([bias, bias_sd])).all() or (bias_sd <= 0).any():
        raise ValueError("LM Arena bias estimates and standard errors must be finite, with positive standard errors.")

    atoms, weights = fit_npmle_bias_prior(bias, cov_mats, u=args.npmle_grid_size)
    hat_mu_b, hat_var_b = fit_parametric_bias_prior(bias, cov_mats)
    npmle_prior_pt = pd.DataFrame({"atom": atoms, "weight": weights})
    param_prior_pt = pd.DataFrame({"hat_mu_b": [hat_mu_b], "hat_var_b": [hat_var_b]})
    bias_se_pt = pd.DataFrame({"bias": bias, "bias_sd": bias_sd})

    args.out_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "lmarena_npmle_pt.csv": npmle_prior_pt,
        "lmarena_param_pt.csv": param_prior_pt,
        "lmarena_bias_se_pt.csv": bias_se_pt,
    }
    for filename, frame in outputs.items():
        path = args.out_dir / filename
        frame.to_csv(path, index=False)
        print(f"Saved {path}")
    print(f"LMArena: {data.M} tasks, split={args.train_test_split:g}, seed={args.seed}, NPMLE grid={args.npmle_grid_size}")


if __name__ == "__main__":
    main()
