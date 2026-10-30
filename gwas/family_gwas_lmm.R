# Family-GWAS with one sibling random-intercept variance component.
#
# Model for individual i in family f at SNP j:
#   y_i = w_i' gamma + delta_j * g_ij +
#         alpha_j * g_parent_fj + u_f + e_i
#
# where u_f ~ N(0, sigma2_sib), e_i ~ N(0, sigma2_e), and
# Cov(y_i, y_k) = sigma2_sib when i and k are siblings.

.solve_chol <- function(chol_a, b) {
  backsolve(
    chol_a,
    forwardsolve(t(chol_a), b)
  )
}

.apply_sib_correlation_inverse <- function(
    x,
    sibling_correlation,
    family_index,
    family_size) {
  # Invert R = (1 - rho) I + rho K by family block, where
  # K[i, k] = 1 when individuals i and k have the same family_id.
  was_vector <- is.null(dim(x))
  if (was_vector) {
    x <- matrix(x, ncol = 1L)
  }

  family_sum <- rowsum(
    x,
    group = family_index,
    reorder = FALSE
  )
  family_weight <- sibling_correlation / (
    1 - sibling_correlation +
      family_size * sibling_correlation
  )
  adjustment <- sweep(
    family_sum[family_index, , drop = FALSE],
    MARGIN = 1L,
    STATS = family_weight[family_index],
    FUN = "*"
  )
  answer <- (
    x - adjustment
  ) / (1 - sibling_correlation)

  singleton <- family_size[family_index] == 1L
  answer[singleton, ] <- x[singleton, , drop = FALSE]

  if (was_vector) {
    return(drop(answer))
  }
  answer
}

.make_fixed_effect_matrix <- function(n, covariates) {
  if (is.null(covariates)) {
    fixed <- matrix(1, nrow = n, ncol = 1L)
    colnames(fixed) <- "intercept"
    return(fixed)
  }

  covariates <- as.matrix(covariates)
  if (!is.numeric(covariates) || nrow(covariates) != n) {
    stop("covariates must be a numeric matrix with n rows.")
  }
  if (anyNA(covariates) ||
      any(!is.finite(covariates))) {
    stop("covariates cannot contain missing values.")
  }

  fixed <- cbind(intercept = 1, covariates)
  if (qr(fixed)$rank != ncol(fixed)) {
    stop(
      "The intercept and covariates are not full rank. ",
      "Do not include an intercept column in covariates."
    )
  }
  fixed
}

.fit_sibling_null_lmm <- function(
    phenotype,
    family_id,
    covariates = NULL,
    rho_upper = 1 - 1e-8) {
  n <- length(phenotype)

  if (!is.numeric(phenotype) ||
      anyNA(phenotype) ||
      any(!is.finite(phenotype))) {
    stop("phenotype must be a complete numeric vector.")
  }
  if (length(family_id) != n || anyNA(family_id)) {
    stop("family_id must be a complete vector of length n.")
  }

  family_id <- as.character(family_id)
  family_index <- match(family_id, unique(family_id))
  family_size <- tabulate(family_index)
  if (!any(family_size >= 2L)) {
    stop(
      "At least one family must contain two or more siblings ",
      "to estimate the sibling variance component."
    )
  }

  fixed <- .make_fixed_effect_matrix(n, covariates)
  n_fixed <- ncol(fixed)
  residual_df <- n - n_fixed
  if (residual_df <= 2L) {
    stop("There are too few residual degrees of freedom.")
  }

  profile_reml <- function(rho, return_fit = FALSE) {
    rinv_y <- .apply_sib_correlation_inverse(
      phenotype,
      sibling_correlation = rho,
      family_index = family_index,
      family_size = family_size
    )
    rinv_fixed <- .apply_sib_correlation_inverse(
      fixed,
      sibling_correlation = rho,
      family_index = family_index,
      family_size = family_size
    )

    fixed_information <- crossprod(fixed, rinv_fixed)
    fixed_chol <- tryCatch(
      chol(fixed_information),
      error = function(error) NULL
    )
    if (is.null(fixed_chol)) {
      return(Inf)
    }

    fixed_score <- crossprod(fixed, rinv_y)
    fixed_estimate <- .solve_chol(fixed_chol, fixed_score)
    p_y <- rinv_y - drop(rinv_fixed %*% fixed_estimate)
    residual_quadratic <- sum(phenotype * p_y)

    if (!is.finite(residual_quadratic) ||
        residual_quadratic <= 0) {
      return(Inf)
    }

    total_variance <- residual_quadratic / residual_df
    logdet_correlation <-
      sum(
        (family_size - 1L) * log1p(-rho) +
          log1p((family_size - 1L) * rho)
      )
    logdet_fixed_information <-
      2 * sum(log(diag(fixed_chol)))

    objective <-
      logdet_correlation +
      logdet_fixed_information +
      residual_df * log(total_variance)

    if (!return_fit) {
      return(objective)
    }

    list(
      objective = objective,
      sibling_correlation = rho,
      total_variance = total_variance,
      sibling_variance = rho * total_variance,
      residual_variance = (1 - rho) * total_variance,
      family_index = family_index,
      family_size = family_size,
      fixed = fixed,
      rinv_fixed = rinv_fixed,
      fixed_chol = fixed_chol,
      p_y = p_y
    )
  }

  optimized <- stats::optimize(
    f = profile_reml,
    interval = c(0, rho_upper),
    tol = 1e-8
  )

  candidates <- c(0, optimized$minimum, rho_upper)
  candidate_objective <- vapply(
    candidates,
    profile_reml,
    numeric(1)
  )
  rho_hat <- candidates[[which.min(candidate_objective)]]
  profile_reml(rho_hat, return_fit = TRUE)
}

#' Run a sibling-random-intercept family GWAS.
#'
#' @param genotype Numeric array with dimensions n x p x 2. The third
#'   dimension is ordered as proband genotype, combined parental genotype.
#' @param phenotype Numeric phenotype vector of length n.
#' @param family_id Length-n pedigree family or sibship identifier. Individuals
#'   with the same value share the sibling random effect. Singleton values are
#'   allowed, but at least one family must contain two or more siblings.
#' @param covariates Optional numeric n-row matrix. Do not include an intercept.
#' @param snp_id Optional vector of p SNP identifiers.
#'
#' @return A list with:
#'   - results: one row per SNP with direct, average non-transmitted, and
#'     population effects; their standard errors; sampling correlations;
#'     and an estimability flag;
#'   - variance_components: sibling, residual, and total variance;
#'   - sibling_correlation: fitted residual phenotype correlation of siblings.
family_gwas_lmm <- function(
    genotype,
    phenotype,
    family_id,
    covariates = NULL,
    snp_id = NULL) {
  genotype_dim <- dim(genotype)
  if (length(genotype_dim) != 3L ||
      genotype_dim[[3L]] != 2L) {
    stop("genotype must be an n by p by 2 numeric array.")
  }
  if (!is.numeric(genotype) ||
      anyNA(genotype) ||
      any(!is.finite(genotype))) {
    stop("genotype must be a complete numeric array.")
  }
  if (any(genotype[, , 1L] < 0 | genotype[, , 1L] > 2) ||
      any(genotype[, , 2L] < 0 | genotype[, , 2L] > 4)) {
    stop(
      "Proband dosages must be in [0, 2] and combined parental ",
      "dosages must be in [0, 4]."
    )
  }

  n <- genotype_dim[[1L]]
  n_snps <- genotype_dim[[2L]]
  if (length(phenotype) != n ||
      length(family_id) != n) {
    stop("phenotype and family_id must both have length n.")
  }

  if (is.null(snp_id)) {
    snp_id <- seq_len(n_snps)
  }
  if (length(snp_id) != n_snps) {
    stop("snp_id must have length p.")
  }

  null_model <- .fit_sibling_null_lmm(
    phenotype = phenotype,
    family_id = family_id,
    covariates = covariates
  )

  delta <- rep(NA_real_, n_snps)
  alpha <- rep(NA_real_, n_snps)
  delta_se <- rep(NA_real_, n_snps)
  alpha_se <- rep(NA_real_, n_snps)
  delta_alpha_cor <- rep(NA_real_, n_snps)
  population_beta <- rep(NA_real_, n_snps)
  population_se <- rep(NA_real_, n_snps)
  direct_population_cor <- rep(NA_real_, n_snps)

  for (snp in seq_len(n_snps)) {
    snp_design <- cbind(
      proband = genotype[, snp, 1L],
      parent_sum = genotype[, snp, 2L]
    )

    rinv_snp <- .apply_sib_correlation_inverse(
      snp_design,
      sibling_correlation =
        null_model$sibling_correlation,
      family_index = null_model$family_index,
      family_size = null_model$family_size
    )
    fixed_snp_information <- crossprod(
      null_model$fixed,
      rinv_snp
    )
    p_snp <- rinv_snp -
      null_model$rinv_fixed %*%
      .solve_chol(
        null_model$fixed_chol,
        fixed_snp_information
      )

    snp_information <- crossprod(snp_design, p_snp)
    snp_information <- 0.5 * (
      snp_information + t(snp_information)
    )
    if (any(!is.finite(snp_information)) ||
        any(diag(snp_information) <= 0) ||
        rcond(snp_information) < 1e-10) {
      next
    }
    snp_chol <- tryCatch(
      chol(snp_information),
      error = function(error) NULL
    )
    if (is.null(snp_chol)) {
      next
    }

    snp_score <- crossprod(
      snp_design,
      null_model$p_y
    )
    estimate <- drop(
      .solve_chol(snp_chol, snp_score)
    )
    estimate_covariance <-
      null_model$total_variance *
      chol2inv(snp_chol)
    estimate_se <- sqrt(diag(estimate_covariance))
    population_variance <-
      estimate_covariance[1L, 1L] +
      estimate_covariance[2L, 2L] +
      2 * estimate_covariance[1L, 2L]
    if (!is.finite(population_variance) ||
        population_variance <= 0) {
      next
    }
    current_population_se <- sqrt(population_variance)

    delta[[snp]] <- estimate[[1L]]
    alpha[[snp]] <- estimate[[2L]]
    delta_se[[snp]] <- estimate_se[[1L]]
    alpha_se[[snp]] <- estimate_se[[2L]]
    delta_alpha_cor[[snp]] <-
      estimate_covariance[1L, 2L] /
        (estimate_se[[1L]] * estimate_se[[2L]])
    population_beta[[snp]] <- sum(estimate)
    population_se[[snp]] <- current_population_se
    direct_population_cor[[snp]] <-
      (
        estimate_covariance[1L, 1L] +
          estimate_covariance[1L, 2L]
      ) /
        (estimate_se[[1L]] * current_population_se)
  }

  list(
    results = data.frame(
      snp_id = snp_id,
      direct_Beta = delta,
      direct_SE = delta_se,
      avg_NTC_Beta = alpha,
      avg_NTC_SE = alpha_se,
      population_Beta = population_beta,
      population_SE = population_se,
      r_direct_avg_NTC = delta_alpha_cor,
      r_direct_population = direct_population_cor,
      estimable = is.finite(delta) & is.finite(alpha),
      stringsAsFactors = FALSE
    ),
    variance_components = c(
      sibling = null_model$sibling_variance,
      residual = null_model$residual_variance,
      total = null_model$total_variance
    ),
    sibling_correlation =
      null_model$sibling_correlation
  )
}
