# Height GWAS analysis and simulation study under realistic LD

Run all commands from this folder.

## Software

Install R, Quarto, and the R packages `data.table`, `MASS`, `ggplot2`,
`reshape2`, `knitr`, `progress`, `REBayes`, `latex2exp`, and `rmarkdown`.
NPMLE fitting also requires a working `Rmosek`/MOSEK installation and license.

## Height GWAS

### 1. Download the inputs

- [Height summary statistics](https://thessgac.com): choose "Young height"
  under "Summary Statistics for Guan et al. (2025)" and place
  `merged_nofilter.sumstats` beside `run_rb.qmd`.
- [1000G EUR PLINK reference](http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz):
  extract under `ldref/1000G/`.
- [PLINK 2](https://www.cog-genomics.org/plink/2.0/).

### 2. Run the analysis

Install `ieugwasr`. For `analysis.qmd`, generate a token on the
[OpenGWAS account page](https://api.opengwas.io/profile/), add
`OPENGWAS_JWT=<token>` to `~/.Renviron`, and restart R.

Run:

```bash
export PLINK2=/path/to/plink2
quarto render run_rb.qmd
quarto render analysis.qmd
```

Open `run_rb.html` and `analysis.html` to view the results. Discovery sets
are saved in `results/discovery_sets_bh/bh/`, and comparison outputs under
`results/`.

### 3. Export data for figures

Running `run_rb.qmd` generates the data for Figures 2 and S3 and saves it to
the CSV files listed below under `results/`.

| File | Columns |
| --- | --- |
| `gwas_bias_se_pt.csv` | `bias`, `bias_sd` |
| `gwas_npmle_pt.csv` | `atom`, `weight` |
| `gwas_param_pt.csv` | `hat_mu_b`, `hat_var_b` |
| `height_gwas_bias_nuisance.csv` | `bias_hat`, `rho`, `tau2`, `sigma2` |

## Simulation study under realistic LD

### 1. Install the genotype-simulation software

Run in R:

```r
dir.create(".Rlib-sim1000G", showWarnings = FALSE)
install.packages("remotes")
install.packages(
  "https://stat.ethz.ch/CRAN/src/contrib/00Archive/hapsim/hapsim_0.31.tar.gz",
  repos = NULL, type = "source", lib = ".Rlib-sim1000G"
)
remotes::install_github("adimitromanolakis/sim1000G@a07bd8f",
  lib = ".Rlib-sim1000G", upgrade = "never")
```

### 2. Download the reference panel

The script downloads the
[1000 Genomes Phase 3 (GRCh37) reference](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)
and [genetic maps](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip)
for chromosomes 18–22.

Install `bcftools`, then run:

```bash
bash download_reference.sh
```

Prepared VCFs and maps are saved in `data/simulation/reference/`.

### 3. Simulate parents and siblings

```bash
Rscript --vanilla simulate_genotypes.R 18 19 20 21 22
```

Files are saved in `data/simulation/` as `families.chr18.rds` through
`families.chr22.rds`.

### 4. Run the experiment

```bash
quarto render simulation.qmd
```

Open `simulation.html` to view the results. Summary tables are saved as
`results/simulation/summary.tsv` and `results/simulation/pvalue_tails.tsv`.
