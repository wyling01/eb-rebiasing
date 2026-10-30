"""Per-task summaries for the bias diagnostic plots."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from pas.datasets.lmarena import LMArenaDataset

from eb_rebiasing._ppi_helpers import pt_decomposition, pt_lambdas
from eb_rebiasing.data_paths import LMARENA_CSV


def load_lmarena_pseudo_oracle(
    csv_path=LMARENA_CSV, *, train_test_split: float = 0.2,
    split_seed: int = 42, pseudo_seed: int = 42,
) -> pd.DataFrame:
    """Create independent pseudo-target and pseudo-bias estimates per task.

    Each task is loaded using the requested LMArena example split, reconstructed
    to its full observation set, and randomly divided 50/50. One half estimates
    the human-preference target; the other estimates prediction minus outcome.
    """
    data = LMArenaDataset(
        file_path=Path(csv_path), train_test_split=train_test_split,
        split_seed=split_seed,
    )
    rng = np.random.default_rng(pseudo_seed)
    rows = []
    for task_id in range(data.M):
        y = np.concatenate([data.y_labelled[task_id], data.y_unlabelled[task_id]])
        pred = np.concatenate([data.pred_labelled[task_id], data.pred_unlabelled[task_id]])
        if len(y) < 2 or len(y) != len(pred) or not np.isfinite([*y, *pred]).all():
            raise ValueError(f"Invalid observations for LMArena task {task_id}")
        order = rng.permutation(len(y))
        n_bias = len(y) // 2
        bias_idx, theta_idx = order[:n_bias], order[n_bias:]
        rows.append({
            "task_id": task_id,
            "theta": float(np.mean(y[theta_idx])),
            "bias": float(np.mean(pred[bias_idx] - y[bias_idx])),
            "n_theta": len(theta_idx),
            "n_bias": len(bias_idx),
        })
    result = pd.DataFrame(rows)
    result.attrs.update({
        "file_path": str(data.file_path),
        "train_test_split": train_test_split,
        "split_seed": split_seed,
        "pseudo_seed": pseudo_seed,
    })
    return result


def nuisance_summary(data) -> pd.DataFrame:
    """Plug-in PT summary; sigma2 refers to the unlabelled prediction mean."""
    lambdas = pt_lambdas(data.pred_unlabelled, data.pred_labelled, data.y_labelled)
    theta_b, bias, pt, cov = pt_decomposition(
        data.pred_unlabelled, data.pred_labelled, data.y_labelled, lambdas,
    )
    sigma2 = np.array([np.var(p, ddof=1) / len(p) for p in data.pred_unlabelled])
    tau2 = np.maximum(cov[:, 1, 1], 0)
    # Disjoint labelled/unlabelled samples give Cov(theta_b, bias_hat) = (1-lambda)*sigma2.
    cross_cov = (1-lambdas) * sigma2
    denom = np.sqrt(sigma2 * tau2)
    rho = np.divide(cross_cov, denom, out=np.full_like(denom, np.nan), where=denom > 0)
    result = pd.DataFrame({"task_id": np.arange(data.M), "theta_b": theta_b,
                           "bias_hat": bias, "pt": pt, "lambda": lambdas,
                           "sigma2": sigma2, "tau2": tau2, "rho": np.clip(rho, -1, 1),
                           "n_labeled": data.ns, "n_unlabeled": data.Ns})
    result.attrs.update({"file_path": str(data.file_path), "split_seed": data.split_seed,
                         "train_test_split": data.train_test_split})
    return result


def load_lmarena_nuisance_summary(
    csv_path=LMARENA_CSV, *, train_test_split: float = 0.2, split_seed: int = 42,
) -> pd.DataFrame:
    return nuisance_summary(LMArenaDataset(file_path=Path(csv_path),
                                          train_test_split=train_test_split, split_seed=split_seed))


__all__ = ["load_lmarena_pseudo_oracle", "load_lmarena_nuisance_summary",
           "nuisance_summary"]
