"""Focused CI-benchmark runner extracted from the upstream PAS package."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from pas.datasets.dataset import PasDataset
from pas.utils import get_avg_ci_width, get_coverage, get_default_args


BENCHMARK_FOLDER = Path(__file__).parent.parent / "results"


def _ci_run_chunk(args):
    """Run one chunk of CI trials in a spawned worker process."""
    import copy
    from pas.intervals import CORE_CI_METHODS

    dataset, ci_method_names, ci_kwargs, alpha, base_seed, trial_indices = args
    dataset = copy.deepcopy(dataset)
    ci_methods = {key: CORE_CI_METHODS[key] for key in ci_method_names}

    rows = {}
    for i in trial_indices:
        dataset.reload_data(split_seed=i + base_seed)
        row = {}
        for ci_name, ci_func in ci_methods.items():
            kwargs = ci_kwargs.get(ci_name, {}) if ci_kwargs else {}
            ci = ci_func(dataset, alpha=alpha, **kwargs)
            row[f"{ci_name}_coverage"] = get_coverage(dataset.true_theta, ci)
            row[f"{ci_name}_width"] = get_avg_ci_width(ci)
        rows[i] = row
    return rows


def run_ci_benchmark(
    dataset: PasDataset,
    trials: int = 100,
    alpha: float = 0.1,
    summary: bool = True,
    save_results: bool = False,
    ci_methods: Optional[dict] = None,
    ci_kwargs: Optional[dict] = None,
    num_workers: int = 1,
    verbose: bool = True,
) -> pd.DataFrame:
    """Run repeated confidence-interval experiments on a PAS dataset."""
    from pas.intervals import CORE_CI_METHODS

    ci_methods = ci_methods or CORE_CI_METHODS
    columns = [
        column
        for name in ci_methods
        for column in (f"{name}_coverage", f"{name}_width")
    ]
    ci_results = pd.DataFrame(columns=columns)
    base_seed = dataset.split_seed

    if num_workers <= 1:
        for i in tqdm(range(trials), disable=not verbose):
            dataset.reload_data(split_seed=i + base_seed)
            for ci_name, ci_func in ci_methods.items():
                kwargs = ci_kwargs.get(ci_name, {}) if ci_kwargs else {}
                ci = ci_func(dataset, alpha=alpha, **kwargs)
                ci_results.loc[i, f"{ci_name}_coverage"] = get_coverage(
                    dataset.true_theta, ci
                )
                ci_results.loc[i, f"{ci_name}_width"] = get_avg_ci_width(ci)
    else:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        import multiprocessing as mp

        chunks = [
            list(map(int, chunk))
            for chunk in np.array_split(np.arange(trials), num_workers)
            if len(chunk) > 0
        ]
        context = mp.get_context("spawn")
        with ProcessPoolExecutor(max_workers=num_workers, mp_context=context) as pool:
            futures = [
                pool.submit(
                    _ci_run_chunk,
                    (
                        dataset,
                        list(ci_methods),
                        ci_kwargs,
                        alpha,
                        base_seed,
                        chunk,
                    ),
                )
                for chunk in chunks
            ]
            with tqdm(
                total=trials,
                disable=not verbose,
                desc=f"alpha={alpha} ({num_workers} workers)",
            ) as progress:
                for future in as_completed(futures):
                    chunk_rows = future.result()
                    for i, row in chunk_rows.items():
                        for column, value in row.items():
                            ci_results.loc[i, column] = value
                    progress.update(len(chunk_rows))
        ci_results = ci_results.sort_index()

    summary_text = [
        f"\nCI Results for {dataset.dataset_name} (alpha={alpha}):\n",
        f"Mean of metrics:\n {ci_results.mean().to_string(float_format='{:.8f}'.format)}",
        f"SD of metrics:\n {ci_results.std().to_string(float_format='{:.8f}'.format)}",
        f"SE of mean metrics:\n {ci_results.sem().to_string(float_format='{:.8f}'.format)}",
    ]
    if summary:
        print("\n".join(summary_text))

    if save_results:
        timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M")
        save_dir = BENCHMARK_FOLDER / f"{dataset.dataset_name}_ci_{timestamp}"
        save_dir.mkdir(parents=True, exist_ok=True)
        ci_results.to_csv(save_dir / "ci_results.csv")
        ci_args = {}
        for name, func in ci_methods.items():
            default_args = get_default_args(func)
            if ci_kwargs and name in ci_kwargs:
                default_args.update(ci_kwargs[name])
            ci_args[name] = default_args
        (save_dir / "ci_args.json").write_text(json.dumps(ci_args, indent=2) + "\n")
        (save_dir / "summary.txt").write_text("\n".join(summary_text) + "\n")

    return ci_results
