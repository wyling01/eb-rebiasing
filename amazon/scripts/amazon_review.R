# --------------
# Amazon review
# --------------
library(dplyr)
library(readr)
library(tibble)

source("../functions/function_rebiased.R")
source("../functions/function_load_data.R")

file <- "../data/amazon_review.h5"

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

# ---- randomly split for K=200 times ----
message("\n=== Running alpha = 0.05 ===")
result1 <- run_repeated_simulation_pt(
  file, n_repeats = 200, alpha = 0.05, seed = 123,
  num_workers = num_workers
)

message("\n=== Running alpha = 0.10 ===")
result2 <- run_repeated_simulation_pt(
  file, n_repeats = 200, alpha = 0.1, seed = 123,
  num_workers = num_workers
)

message("\n=== Running alpha = 0.20 ===")
result3 <- run_repeated_simulation_pt(
  file, n_repeats = 200, alpha = 0.2, seed = 123,
  num_workers = num_workers
)

message("\n=== Running alpha = 0.30 ===")
result4 <- run_repeated_simulation_pt(
  file, n_repeats = 200, alpha = 0.3, seed = 123,
  num_workers = num_workers
)

# ---- Save results ----
out_dir <- "../results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

results_list <- list(
  setting1 = result1,
  setting2 = result2,
  setting3 = result3,
  setting4 = result4
)

alpha_lookup <- tibble(
  Setting = paste0("setting", 1:4),
  alpha = c(0.05, 0.1, 0.2, 0.3)
)

summary_results <- bind_rows(results_list, .id = "Setting") %>%
  left_join(alpha_lookup, by = "Setting")

write_csv(summary_results, file.path(out_dir, "amazon_real.csv"))
