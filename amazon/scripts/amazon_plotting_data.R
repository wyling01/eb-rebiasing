# ---------------------------------
# Amazon paper-plotting input data
# ---------------------------------
library(readr)
library(tibble)

source("../functions/function_rebiased.R")
source("../functions/function_load_data.R")

file <- "../data/amazon_review.h5"
out_dir <- "../results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Use one fixed split for reproducible plotting and diagnostic inputs.
sum_stat <- load_data_amazon(file, r = 0.2, seed = 3)
updated_stat_pt <- calculate_pt_and_update(sum_stat)
pt <- get_input_vectors(updated_stat_pt, "lambda_")

parametric_prior_pt <- estimate_parametric_mle(
  pt$delta,
  pt$sed
)
npmle_prior_pt <- GNmix_diffvar(pt$delta, pt$sed, u = 50)

bias_se_pt <- tibble(
  bias = as.numeric(pt$delta),
  bias_sd = as.numeric(pt$sed)
)
write_csv(bias_se_pt, file.path(out_dir, "amazon_bias_se_pt.csv"))

param_prior_pt <- tibble(
  hat_mu_b = parametric_prior_pt$hat_mu_b,
  hat_var_b = parametric_prior_pt$hat_var_b
)
write_csv(param_prior_pt, file.path(out_dir, "amazon_param_pt.csv"))

npmle_prior_output <- tibble(
  atom = as.numeric(as.matrix(npmle_prior_pt$u)),
  weight = as.numeric(npmle_prior_pt$fu)
)
write_csv(npmle_prior_output, file.path(out_dir, "amazon_npmle_pt.csv"))

diagnostic_stat <- make_amazon_diagnostic(
  file = file,
  r = 0.2,
  summary_seed = 3
)
write_csv(
  diagnostic_stat,
  file.path(out_dir, "amazon_bias_nuisance_diagnostics.csv")
)
