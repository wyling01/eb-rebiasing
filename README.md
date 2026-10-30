# Empirical Bayes Rebiasing

Reproduction code for paper *"Empirical Bayes Rebiasing"*.

## Code provenance

For self-contained reproducibility, this repository includes the PAS
(Prediction-Powered Adaptive Shrinkage) components required by our experiments
under `src/pas/`. The original PAS implementation is available at
[Prediction-Powered Adaptive Shrinkage](https://github.com/listar2000/prediction-powered-adaptive-shrinkage).
The local version has been adapted for the confidence-interval experiments and
empirical Bayes comparisons in this paper.

## Setup

Run the commands below from the repository root. The Python code supports
Python 3.9--3.13.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

Python dependencies are declared in `pyproject.toml`. NPMLE procedures require
a working MOSEK installation and license. If the license is not in MOSEK's
default location, set:

```bash
export MOSEKLM_LICENSE_FILE=/absolute/path/to/mosek.lic
```

The Amazon and GWAS applications additionally require R. Their dependencies
and data setup are documented separately:

- [Amazon Fine Food Reviews](amazon/README.md)
- [GWAS](gwas/README.md)

## LMArena

Run a one-trial smoke test, followed by the full experiment:

```bash
python scripts/run_lmarena_ci.py --trials 1 --alphas 0.1 --num-workers 1
python scripts/run_lmarena_ci.py --trials 200 --num-workers 4
```

Results are written to `tmp/lmarena_ci_results/`. Generate the LMArena inputs
for Figure 2 with:

```bash
python scripts/lmarena_plotting_data.py
```

## Data applications and simulations

### Amazon Fine Food Reviews

The Amazon real-data and synthetic experiments are run with R. See
[`amazon/README.md`](amazon/README.md) for dependencies, input data, parallel
execution, and output files.

```bash
cd amazon/scripts
Rscript amazon_review.R --num-workers=4
Rscript amazon_plotting_data.R
Rscript synthetic_amazon.R --num-workers=4
cd ../..
```

### GWAS

The GWAS analysis and realistic-LD simulation require external data and
additional software. See [`gwas/README.md`](gwas/README.md) for the complete
workflow.

## Paper plots

After generating the required result files, run from the repository root:

```bash
python scripts/make_figure2_plot.py
python scripts/make_figure3_plot.py
python scripts/make_figure4_plot.py
python scripts/make_figureS1_plot.py
python scripts/make_diagnostic_plot.py
```

The PDFs are saved under `paper_plots/` by default. Run any plotting script
with `--help` to see its input and output options.

## Repository layout

```text
src/eb_rebiasing/   Empirical Bayes rebiasing methods
src/pas/            Local PAS components used by the experiments
scripts/            LMArena experiment and paper-plotting scripts
amazon/             Amazon real-data and synthetic experiments
gwas/               GWAS analysis and realistic-LD simulation study
data/                Cleaned LMArena input
paper_plots/         Generated paper figures
tmp/                 Intermediate LMArena outputs
```
