# Helpers for the direct-null experiment on the fixed chr18-22 panel.
source("parental_genotype_imputation.R")
source("family_gwas_lmm.R")
source("functions.R")

load_simulation_data <- function(data_dir = "data/simulation") {
  chromosomes <- 18:22
  inputs <- file.path(data_dir, sprintf("families.chr%d.rds", chromosomes))
  if (!all(file.exists(inputs))) {
    stop("Genotypes are not included. Run bash download_reference.sh, then ",
         "Rscript --vanilla simulate_genotypes.R. See README.md.")
  }
  genotype <- true_parent <- vector("list", length(chromosomes))
  variants <- vector("list", length(chromosomes))
  for (i in seq_along(chromosomes)) {
    x <- readRDS(file.path(data_dir, sprintf("families.chr%d.rds", chromosomes[i])))
    if (i == 1L) {
      family_id <- x$family_id
      pedigree <- data.frame(
        family_id = rep(family_id, each = 2L),
        individual_id = paste0(rep(family_id, each = 2L), "_S", 1:2)
      )
    }
    stopifnot(identical(x$family_id, family_id))
    child <- matrix(0L, nrow(pedigree), ncol(x$sibling1_genotype),
                    dimnames = list(pedigree$individual_id, x$variant$snp_id))
    child[seq(1L, nrow(child), by = 2L), ] <- x$sibling1_genotype
    child[seq(2L, nrow(child), by = 2L), ] <- x$sibling2_genotype
    imputed <- impute_parental_genotypes(child, x$ibd_state, pedigree)
    genotype[[i]] <- array(c(child, imputed$parental_alt_dosage),
                           dim = c(nrow(child), ncol(child), 2L))
    true_parent[[i]] <- x$father_genotype + x$mother_genotype
    variants[[i]] <- x$variant
  }
  variant <- do.call(rbind, variants)
  stopifnot(length(family_id) == 5000L, nrow(variant) == 6000L,
            !anyDuplicated(variant$snp_id))
  list(genotype = genotype, true_parent = true_parent,
       pedigree = pedigree, variant = variant,
       snp_id = lapply(variants, function(x) x$snp_id))
}

simulate_null_phenotype <- function(data, parental_variance, replicate) {
  set.seed(20260726L + replicate, kind = "Mersenne-Twister",
           normal.kind = "Inversion", sample.kind = "Rejection")
  alpha <- rnorm(nrow(data$variant))
  alpha_by_chr <- split(alpha, rep(seq_along(data$snp_id), lengths(data$snp_id)))
  parent_score <- Reduce(`+`, Map(function(g, a) drop(g %*% a),
                                  data$true_parent, alpha_by_chr))
  parent_score <- rep(parent_score, each = 2L)
  parent_score <- parent_score * sqrt(parental_variance / var(parent_score))
  residual <- rnorm(length(parent_score))
  residual <- residual * sqrt((1 - var(parent_score)) / var(residual))
  list(y = parent_score + residual, parent_score = parent_score, residual = residual)
}

simulation_gwas <- function(data, y) {
  result <- lapply(seq_along(data$genotype), function(i) {
    family_gwas_lmm(data$genotype[[i]], y, data$pedigree$family_id,
                    snp_id = data$snp_id[[i]])$results
  })
  result <- do.call(rbind, result)
  stopifnot(all(result$estimable))
  result
}

simulation_pvalues <- function(sumstats) {
  npmle <- GNmix_diffvar(sumstats$avg_NTC_Beta, sumstats$avg_NTC_SE,
                         u = 300, verbose = FALSE)
  normal <- estimate_parametric_mle(sumstats$avg_NTC_Beta, sumstats$avg_NTC_SE)
  p <- list(
    npmle = vapply(seq_len(nrow(sumstats)), function(j) {
      cal_p_value_npmle_diffvar(as.matrix(npmle$u), npmle$fu,
        sumstats$direct_Beta[j], sumstats$avg_NTC_Beta[j],
        sumstats$direct_SE[j], sumstats$avg_NTC_SE[j], sumstats$r_direct_avg_NTC[j])
    }, numeric(1)),
    parametric = vapply(seq_len(nrow(sumstats)), function(j) {
      cal_p_value_param_diffvar(normal$hat_mu_b, normal$hat_var_b,
        sumstats$direct_Beta[j], sumstats$avg_NTC_Beta[j],
        sumstats$direct_SE[j], sumstats$avg_NTC_SE[j], sumstats$r_direct_avg_NTC[j])
    }, numeric(1)),
    direct_z = cal_p_value_z_diffvar(sumstats$direct_Beta, sumstats$direct_SE),
    population_z = cal_p_value_z_diffvar(sumstats$population_Beta, sumstats$population_SE)
  )
  stopifnot(npmle$status == "OPTIMAL", normal$convergence == 0L,
            all(is.finite(unlist(p))), all(unlist(p) >= 0 & unlist(p) <= 1))
  list(p_values = p, npmle = npmle, normal = normal)
}

summarize_null <- function(p_values, q = 0.05) {
  do.call(rbind, lapply(names(p_values), function(method) {
    p <- p_values[[method]]
    counts <- apply(p, 1L, function(x) sum(p.adjust(x, "BH") <= q))
    n <- nrow(p)
    ci <- binom.test(sum(counts > 0), n)$conf.int
    prop <- rowMeans(p < 0.05)
    half_width <- if (n > 1) qt(0.975, n - 1) * sd(prop) / sqrt(n) else NA_real_
    data.frame(method = method, replicates = n, bh_level = q,
      rejection_replicates = sum(counts > 0), fdr = mean(counts > 0),
      fdr_ci_lower = ci[1], fdr_ci_upper = ci[2],
      mean_bh_rejections = mean(counts), maximum_bh_rejections = max(counts),
      p_lt_0_05_proportion = mean(prop),
      p_lt_0_05_ci_lower = mean(prop) - half_width,
      p_lt_0_05_ci_upper = mean(prop) + half_width)
  }))
}

summarize_tails <- function(p_values, thresholds = c(0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001)) {
  do.call(rbind, lapply(names(p_values), function(method) {
    do.call(rbind, lapply(thresholds, function(threshold) {
      prop <- rowMeans(p_values[[method]] < threshold)
      n <- length(prop)
      half_width <- if (n > 1) qt(0.975, n - 1) * sd(prop) / sqrt(n) else NA_real_
      data.frame(method = method, threshold = threshold,
        proportion = mean(prop), ci_lower = mean(prop) - half_width,
        ci_upper = mean(prop) + half_width)
    }))
  }))
}

run_null_study <- function(data, parental_variance, n_replicates, output_dir) {
  if (dir.exists(output_dir)) stop("Choose a new output directory: ", output_dir)
  dir.create(output_dir, recursive = TRUE)
  methods <- c("npmle", "parametric", "direct_z", "population_z")
  p_values <- setNames(lapply(methods, function(x) {
    matrix(NA_real_, n_replicates, nrow(data$variant),
           dimnames = list(seq_len(n_replicates), data$variant$snp_id))
  }), methods)
  prior <- phenotype_qc <- vector("list", n_replicates)
  for (replicate in seq_len(n_replicates)) {
    phenotype <- simulate_null_phenotype(data, parental_variance, replicate)
    sumstats <- simulation_gwas(data, phenotype$y)
    fit <- simulation_pvalues(sumstats)
    for (method in methods) p_values[[method]][replicate, ] <- fit$p_values[[method]]
    prior[[replicate]] <- data.frame(replicate,
      normal_mean = fit$normal$hat_mu_b, normal_variance = fit$normal$hat_var_b,
      normal_convergence = fit$normal$convergence, npmle_status = fit$npmle$status)
    phenotype_qc[[replicate]] <- data.frame(replicate,
      parental_variance = var(phenotype$parent_score),
      residual_variance = var(phenotype$residual), phenotype_variance = var(phenotype$y))
    if (replicate == 1L) {
      data.table::fwrite(sumstats, file.path(output_dir, "sumstats.rep1.tsv"), sep = "\t")
      saveRDS(fit[c("npmle", "normal")], file.path(output_dir, "priors.rep1.rds"))
    }
    if (replicate == 1L || replicate %% 10L == 0L)
      message("Parental variance ", parental_variance, ": ", replicate, "/", n_replicates)
  }
  saveRDS(p_values, file.path(output_dir, "p_values.rds"))
  data.table::fwrite(summarize_null(p_values), file.path(output_dir, "summary.tsv"), sep = "\t")
  data.table::fwrite(summarize_tails(p_values), file.path(output_dir, "pvalue_tails.tsv"), sep = "\t")
  data.table::fwrite(do.call(rbind, prior), file.path(output_dir, "prior_summary.tsv"), sep = "\t")
  data.table::fwrite(do.call(rbind, phenotype_qc), file.path(output_dir, "phenotype_qc.tsv"), sep = "\t")
  saveRDS(list(parental_variance = parental_variance, replicates = seq_len(n_replicates),
               seed = 20260726L, chromosomes = 18:22, bh_level = 0.05),
          file.path(output_dir, "settings.rds"))
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  invisible(p_values)
}
