# --------------
# Synthetic dataset built upon Amazon review
# --------------
library(dplyr)
library(readr)
library(purrr)
library(tibble)

source("../functions/function_rebiased.R")
source("../functions/function_load_data.R")

file <-  "../data/amazon_review.h5"
sum_stat <- load_data_amazon(file, r = 0.2, seed = 123)
updated_stat_pt <- calculate_pt_and_update(sum_stat)

args <- commandArgs(trailingOnly = TRUE)
num_workers_arg <- grep("^--num-workers=", args, value = TRUE)
if (length(num_workers_arg) > 1L) {
  stop("Specify --num-workers only once.")
}
num_workers <- if (length(num_workers_arg) == 1L) {
  suppressWarnings(as.integer(sub("^--num-workers=", "", num_workers_arg)))
} else {
  4L
}
if (is.na(num_workers) || num_workers < 1L) {
  stop("--num-workers must be a positive integer.")
}
message("Using ", num_workers, " simulation worker(s).")

out_dir <- "../results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Vary b ----
## ---- Normal ----
message("\n=== Vary bias distribution: Normal ===")
result1_normal <- run_repeated_simulation_synthetic_pt(updated_stat_pt,n_repeats = 500,
                                                       M = 200, u_val = 50,
                                                       mu_theta0 = 4, sd_theta0 = 0.5,
                                                       mu_b = -0.1,sd_b = 0.01, alpha = 0.05,
                                                       seed = 1, num_workers = num_workers) 
result2_normal <- run_repeated_simulation_synthetic_pt(updated_stat_pt,n_repeats = 500,
                                                       M = 200, u_val = 50,
                                                       mu_theta0 = 4, sd_theta0 = 0.5,
                                                       mu_b = -0.1,sd_b = 0.03, alpha = 0.05,
                                                       seed = 1, num_workers = num_workers) 
result3_normal <- run_repeated_simulation_synthetic_pt(updated_stat_pt,n_repeats = 500,
                                                       M = 200, u_val = 50,
                                                       mu_theta0 = 4, sd_theta0 = 0.5,
                                                       mu_b = -0.1,sd_b = 0.05, alpha = 0.05,
                                                       seed = 1, num_workers = num_workers) 
result4_normal <- run_repeated_simulation_synthetic_pt(updated_stat_pt,n_repeats = 500,
                                                       M = 200, u_val = 50,
                                                       mu_theta0 = 4, sd_theta0 = 0.5,
                                                       mu_b = -0.1,sd_b = 0.1, alpha = 0.05,
                                                       seed = 1, num_workers = num_workers) 
normal_results_list <- list(
  `0.01` = result1_normal,
  `0.03` = result2_normal,
  `0.05` = result3_normal,
  `0.10` = result4_normal
)

normal_summary_results <- imap_dfr(
  normal_results_list,
  function(x, sigma_value) {
    df <- as.data.frame(x)
    
    if (!"Method" %in% colnames(df)) {
      df <- rownames_to_column(df, var = "Method")
    }
    
    df %>%
      as_tibble() %>%
      mutate(
        sigma = as.numeric(sigma_value),
        .before = Method
      )
  }
)

write_csv(normal_summary_results,file.path(out_dir, "synthetic_normal.csv"))
                                              
## ---- Two Point ----
message("\n=== Vary bias distribution: Two-point ===")
result1_twopoint <- run_repeated_simulation_synthetic_pt_twopoint(updated_stat_pt,n_repeats = 500,
                                                                  M = 200, u_val = 50,
                                                                  mu_theta0 = 4, sd_theta0 = 0.5,
                                                                  value = c(0,0.05), probability = c(0.5,0.5), alpha = 0.05,
                                                                  seed = 1, num_workers = num_workers) 
result2_twopoint <- run_repeated_simulation_synthetic_pt_twopoint(updated_stat_pt,n_repeats = 500,
                                                                  M = 200, u_val = 50,
                                                                  mu_theta0 = 4, sd_theta0 = 0.5,
                                                                  value = c(0,0.1), probability = c(0.5,0.5), alpha = 0.05,
                                                                  seed = 1, num_workers = num_workers) 
result3_twopoint <- run_repeated_simulation_synthetic_pt_twopoint(updated_stat_pt,n_repeats = 500,
                                                                  M = 200, u_val = 50,
                                                                  mu_theta0 = 4, sd_theta0 = 0.5,
                                                                  value = c(0,0.2), probability = c(0.5,0.5), alpha = 0.05,
                                                                  seed = 1, num_workers = num_workers) 
result4_twopoint <- run_repeated_simulation_synthetic_pt_twopoint(updated_stat_pt,n_repeats = 500,
                                                                  M = 200, u_val = 50,
                                                                  mu_theta0 = 4, sd_theta0 = 0.5,
                                                                  value = c(0,0.5), probability = c(0.5,0.5), alpha = 0.05,
                                                                  seed = 1, num_workers = num_workers)
twopoint_results_list <- list(
  `0.05` = result1_twopoint,
  `0.10` = result2_twopoint,
  `0.20` = result3_twopoint,
  `0.50` = result4_twopoint
)

twopoint_summary_results <- imap_dfr(
  twopoint_results_list,
  function(x, b0_value) {
    df <- as.data.frame(x)
    
    if (!"Method" %in% colnames(df)) {
      df <- rownames_to_column(df, var = "Method")
    }
    
    df %>%
      as_tibble() %>%
      mutate(
        b0 = as.numeric(b0_value),
        .before = Method
      )
  }
)

write_csv(twopoint_summary_results,file.path(out_dir, "synthetic_twopoint.csv"))

# ---- Vary theta ----
message("\n=== Vary target dispersion ===")
result1_theta <- run_repeated_simulation_synthetic_pt(updated_stat_pt,n_repeats = 500,
                                                       M = 200, u_val = 50,
                                                       mu_theta0 = 4, sd_theta0 = 0.01,
                                                       mu_b = -0.1,sd_b = 0.01, alpha = 0.05,
                                                       seed = 1, num_workers = num_workers) 
result2_theta <- run_repeated_simulation_synthetic_pt(updated_stat_pt,n_repeats = 500,
                                                       M = 200, u_val = 50,
                                                       mu_theta0 = 4, sd_theta0 = 0.1,
                                                       mu_b = -0.1,sd_b = 0.01, alpha = 0.05,
                                                       seed = 1, num_workers = num_workers) 
result3_theta <- result1_normal 
result4_theta <- run_repeated_simulation_synthetic_pt(updated_stat_pt,n_repeats = 500,
                                                       M = 200, u_val = 50,
                                                       mu_theta0 = 4, sd_theta0 = 1,
                                                       mu_b = -0.1,sd_b = 0.01, alpha = 0.05,
                                                       seed = 1, num_workers = num_workers) 

theta_results_list <- list(
  `0.01` = result1_theta,
  `0.10` = result2_theta,
  `0.50` = result3_theta,
  `1.00` = result4_theta
)

theta_summary_results <- imap_dfr(
  theta_results_list,
  function(x, theta_value) {
    df <- as.data.frame(x)
    
    if (!"Method" %in% colnames(df)) {
      df <- rownames_to_column(df, var = "Method")
    }
    
    df %>%
      as_tibble() %>%
      mutate(
        theta = as.numeric(theta_value),
        .before = Method
      )
  }
)

write_csv(theta_summary_results ,file.path(out_dir, "synthetic_theta.csv"))
