.compute_summary <- function(Y, Z, r) {
  Tn <- length(Y)
  
  theta_true <- mean(Y)
  theta_f_true <- mean(Z)

  n_lab <- max(1L, min(Tn - 1L, floor(r * Tn)))
  lab_i  <- sample.int(Tn, n_lab)
  unlab_i <- setdiff(seq_len(Tn), lab_i)
  n_unl <- length(unlab_i)
  
  y_lab <- Y[lab_i]
  z_lab <- Z[lab_i]
  z_unl <- Z[unlab_i]

  Z_tilde <- mean(z_unl)                          
  Z_bar   <- mean(z_lab)                          
  Y_bar   <- mean(y_lab)                          

  b_hat     <- Z_bar - Y_bar
  theta_hat <- Z_tilde - b_hat

  sZ2_unl <- if (n_unl >= 2) var(z_unl) else NA_real_ 
  sY2_lab <- if (n_lab >= 2) var(y_lab) else NA_real_ 
  sZ2_lab <- if (n_lab >= 2) var(z_lab) else NA_real_ 
  cov_lab <- if (n_lab >= 2) cov(y_lab, z_lab) else NA_real_ 
  sR2     <- if (n_lab >= 2) var(y_lab - z_lab) else NA_real_
  
  var_b      <- if (!is.na(sR2)) sR2 / n_lab else NA_real_
  var_Ztilde <- if (!is.na(sZ2_unl)) sZ2_unl / n_unl else NA_real_

  var_theta  <- if (is.na(var_b) | is.na(var_Ztilde)) NA_real_ else (ifelse(is.na(var_b), 0, var_b) + ifelse(is.na(var_Ztilde), 0, var_Ztilde))
  var_Ybar   <- if (!is.na(sY2_lab)) sY2_lab / n_lab else NA_real_
  
  corr_theta_b <- if (!is.na(var_theta) && !is.na(var_b) && var_theta > 0 && var_b >= 0) {
    -sqrt(var_b / var_theta)
  } else NA_real_
  
  var_f_global <- var(Z) 
  
  data.frame(
    theta_true,
    theta_f_true,
    n_total      = Tn,
    n_labeled    = n_lab,
    n_unlabeled  = n_unl,

    Z_tilde      = Z_tilde,
    Z_bar        = Z_bar,
    Y_bar        = Y_bar,
    b_hat        = b_hat,
    theta_hat    = theta_hat,

    var_theta_hat= var_theta,
    var_b_hat    = var_b,
    corr_theta_b = corr_theta_b,
    var_Z_tilde  = var_Ztilde,
    var_Y_bar    = var_Ybar,

    var_f_global = var_f_global,  
    var_z_unlab    = sZ2_unl,
    var_z_lab    = sZ2_lab,      
    var_y_lab    = sY2_lab,       
    cov_yf_lab   = cov_lab,       
    stringsAsFactors = FALSE
  )
}

calculate_pt_and_update <- function(sum_stat) {
  
  library(dplyr)
  
  updated_stat <- sum_stat %>%
    mutate(
      temp_cov_Ybar_Zbar = cov_yf_lab / n_labeled,
      
      temp_lambda_denom =
        var_f_global * (1 / n_labeled + 1 / n_unlabeled),
      
      lambda_value = pmin(
        pmax(
          ifelse(
            !is.na(temp_lambda_denom) & temp_lambda_denom > 0,
            temp_cov_Ybar_Zbar / temp_lambda_denom,
            NA_real_
          ),
          0
        ),
        1
      ),

      lambda_theta_hat =
        Y_bar + lambda_value * (Z_tilde - Z_bar),

      lambda_b_hat =
        Z_tilde - lambda_theta_hat,

      temp_var_Ybar_local =
        var_y_lab / n_labeled,
      
      temp_var_Zbar_local =
        var_z_lab / n_labeled,
      
      temp_var_Ztilde_local =
        var_z_unlab / n_unlabeled,
      
      temp_cov_Ybar_Zbar_local =
        cov_yf_lab / n_labeled,

      lambda_var_theta_hat =
        temp_var_Ybar_local +
        lambda_value^2 * temp_var_Zbar_local +
        lambda_value^2 * temp_var_Ztilde_local -
        2 * lambda_value * temp_cov_Ybar_Zbar_local,

      lambda_var_b_hat =
        temp_var_Ybar_local +
        lambda_value^2 * temp_var_Zbar_local +
        (1 - lambda_value)^2 * temp_var_Ztilde_local -
        2 * lambda_value * temp_cov_Ybar_Zbar_local,

      temp_cov_theta_b =
        -temp_var_Ybar_local +
        2 * lambda_value * temp_cov_Ybar_Zbar_local -
        lambda_value^2 * temp_var_Zbar_local +
        lambda_value * (1 - lambda_value) * temp_var_Ztilde_local,
      
      lambda_corr_theta_b = ifelse(
        lambda_var_theta_hat > 0 & lambda_var_b_hat > 0,
        temp_cov_theta_b / sqrt(lambda_var_theta_hat * lambda_var_b_hat),
        NA_real_
      )
    ) %>%
    dplyr::select(
      -temp_cov_Ybar_Zbar,
      -temp_lambda_denom,
      -temp_var_Ybar_local,
      -temp_var_Zbar_local,
      -temp_var_Ztilde_local,
      -temp_cov_Ybar_Zbar_local,
      -temp_cov_theta_b
    )
  
  return(updated_stat)
}

get_input_vectors <- function(df, prefix = "") {
  list(
    beta  = df[[paste0(prefix, "theta_hat")]],
    delta = df[[paste0(prefix, "b_hat")]],
    seb   = sqrt(df[[paste0(prefix, "var_theta_hat")]]),
    sed   = sqrt(df[[paste0(prefix, "var_b_hat")]]),
    rho   = df[[paste0(prefix, "corr_theta_b")]]
  )
}

load_data_amazon <- function(file, r, seed = NULL) {
  stopifnot(file.exists(file), r > 0, r < 1)
  if (!is.null(seed)) set.seed(seed)
  ls_tbl <- rhdf5::h5ls(file, recursive = FALSE)
  keys   <- ls_tbl$name[ls_tbl$otype == "H5I_DATASET"]
  if (length(keys) == 0L) stop("No root-level datasets found.")
  
  rows <- lapply(keys, function(k) {
    m <- rhdf5::h5read(file, paste0("/", k))  # exactly your call
    if (!is.matrix(m) || ncol(m) < 2)
      stop(sprintf("Dataset '%s' is not a 2-column matrix.", k))
    Y <- as.numeric(m[, 1])
    Z <- as.numeric(m[, 2])
    cbind(data.frame(group = k, stringsAsFactors = FALSE),
          .compute_summary(Y, Z, r = r))
  })
  
  do.call(rbind, rows)
}



#----Real data analysis----
get_asymmetric_results_pt <- function(df, u_val = 50,alpha=0.05) {
  std <- get_input_vectors(df)
  uni <- get_input_vectors(df, "lambda_")
  
  uni_npmle_prior <- GNmix_diffvar(uni$delta, uni$sed, u = u_val,verbose = FALSE)
  uni_param_prior <- estimate_parametric_mle(uni$delta, uni$sed)
  
  res_list <- list(
    classic        = mapply(cal_ci_normal, df$Y_bar, sqrt(df$var_Y_bar), MoreArgs = list(alpha = alpha)),
    classic_bias        = mapply(cal_ci_normal, df$Z_tilde, sqrt(df$var_Z_tilde), MoreArgs = list(alpha = alpha)),
    ppi            = mapply(cal_ci_normal, std$beta, std$seb, MoreArgs = list(alpha = alpha)),
    pt        = mapply(cal_ci_normal, uni$beta, uni$seb, MoreArgs = list(alpha = alpha)),
    pt_param      = mapply(cal_ci_parametric, uni$beta, uni$delta, uni$seb, uni$sed, uni$rho, 
                           MoreArgs = list(a = uni_param_prior$hat_mu_b, b = uni_param_prior$hat_var_b, alpha = alpha)),
    pt_npmle      = mapply(cal_ci_npmle, uni$beta, uni$delta, uni$seb, uni$sed, uni$rho, 
                           MoreArgs = list(v = as.matrix(uni_npmle_prior$u), p = uni_npmle_prior$fu, alpha = alpha))
  )

  ci_mats <- lapply(res_list, function(x) {
    m <- t(x)
    colnames(m) <- c("lower", "upper")
    return(m)
  })

  tmp <- future_lapply(1:2, function(k) {
    if (k == 1) {
      cal_ci_ebci(df$Y_bar,sqrt(df$var_Y_bar),alpha = alpha)
    } else {
      cal_ci_ebci(uni$beta,uni$seb,alpha = alpha)
    } 
  })
  ci_mats$classic_ebci <- tmp[[1]]
  ci_mats$pt_ebci      <- tmp[[2]]
  
  
  return(list(ci_matrices = ci_mats, priors = list(uni_npmle = uni_npmle_prior, uni_param = uni_param_prior)))
}

run_repeated_simulation_pt <- function(file, n_repeats = 200, r_split = 0.2, u_val = 50,
                                       alpha = 0.05, seed = NULL, num_workers = 1) {
  stopifnot(length(num_workers) == 1L, is.finite(num_workers), num_workers >= 1)
  num_workers <- as.integer(num_workers)

  run_one <- function(i) {
    sum_stat <- load_data_amazon(
      file,
      r = r_split,
      seed = if (is.null(seed)) NULL else seed + i - 1L
    )
    updated_stat <- calculate_pt_and_update(sum_stat)
    
    updated_stat_clean_2 <- updated_stat[
      complete.cases(updated_stat) & 
        updated_stat$var_b_hat > 1e-10, 
    ]
    
    asym_res <- get_asymmetric_results_pt(updated_stat_clean_2, u_val = u_val,alpha=alpha)
    theta_true <- updated_stat_clean_2$theta_true
    
    list(
      coverage = unlist(lapply(asym_res$ci_matrices, cov_rate, theta = theta_true)),
      length = unlist(lapply(asym_res$ci_matrices, function(m) {
        mean(m[, "upper"] - m[, "lower"], na.rm = TRUE)
      })),
      length_ratio = calculate_avg_of_ratios(asym_res)
    )
  }

  if (num_workers == 1L) {
    pb <- txtProgressBar(min = 0, max = n_repeats, style = 3)
    on.exit(close(pb), add = TRUE)
    simulation_results <- lapply(seq_len(n_repeats), function(i) {
      result_i <- run_one(i)
      setTxtProgressBar(pb, i)
      result_i
    })
  } else {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = num_workers)

    old_progressr_enable <- getOption("progressr.enable")
    old_progressr_handlers <- progressr::handlers()
    on.exit({
      options(progressr.enable = old_progressr_enable)
      progressr::handlers(old_progressr_handlers)
    }, add = TRUE)
    options(progressr.enable = TRUE)
    progressr::handlers(progressr::handler_txtprogressbar(clear = FALSE))

    simulation_results <- progressr::with_progress({
      progress <- progressr::progressor(steps = n_repeats)
      progress(amount = 0)
      future.apply::future_lapply(
        seq_len(n_repeats),
        function(i) {
          result_i <- run_one(i)
          progress()
          result_i
        },
        future.seed = NULL,
        future.chunk.size = 1
      )
    })
  }

  trial_cov_asym <- lapply(simulation_results, `[[`, "coverage")
  trial_len_asym <- lapply(simulation_results, `[[`, "length")
  trial_rat_asym <- lapply(simulation_results, `[[`, "length_ratio")
  
  c_mat <- do.call(rbind, trial_cov_asym)
  l_mat <- do.call(rbind, trial_len_asym)
  r_mat <- do.call(rbind, trial_rat_asym)
  
  mc_se <- function(x) {
    n_eff <- sum(!is.na(x))
    if (n_eff <= 1) return(NA_real_)
    sd(x, na.rm = TRUE) / sqrt(n_eff)
  }
  
  df_asymmetric <- data.frame(
    Method = colnames(c_mat),
    SE_Coverage = apply(c_mat, 2, mc_se),
    Avg_Coverage = colMeans(c_mat, na.rm = TRUE),
    SE_Length = apply(l_mat, 2, mc_se),
    Avg_Length = colMeans(l_mat, na.rm = TRUE),
    SE_Length_Ratio = apply(r_mat, 2, mc_se),
    Avg_Length_Ratio = colMeans(r_mat, na.rm = TRUE)
  )
  
  return(df_asymmetric)
}


#----Synthetic data analysis----
simulate_ppi_data_pt <- function(M = 200, updated_stat_pt,
                                 # Latent distributions:
                                 # true_mu_f = theta_i^f
                                 # true_b    = theta_i^f - theta_i
                                 mu_theta0 = 4, sd_theta0 = 0.01,
                                 mu_b = -0.1, sd_b = 0.1,
                                 seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  n_real <- nrow(updated_stat_pt)
  sampled_indices <- sample(seq_len(n_real), size = M, replace = TRUE)
  sim_settings <- updated_stat_pt[sampled_indices, ]
  
  sig_y_lab <- sqrt(sim_settings$var_y_lab)
  sig_z_lab <- sqrt(sim_settings$var_z_lab)
  cov_lab   <- sim_settings$cov_yf_lab
  
  var_mean_unl <- sim_settings$var_Z_tilde
  
  n_lab <- sim_settings$n_labeled
  n_unl <- sim_settings$n_unlabeled
  
  # Problem-specific PT lambda_i
  lambda_vec <- sim_settings$lambda_value
  
  true_theta <- rnorm(M, mean = mu_theta0, sd = sd_theta0)
  true_b     <- rnorm(M, mean = mu_b, sd = sd_b)
  true_mu_f  <- true_theta + true_b
  
  Y_bar_vec   <- numeric(M)
  Z_bar_vec   <- numeric(M)
  Z_tilde_vec <- numeric(M)
  
  for (i in seq_len(M)) {
    
    # Labeled data: (Y_bar, Z_bar)
    Sigma_lab <- matrix(
      c(sig_y_lab[i]^2, cov_lab[i],
        cov_lab[i],     sig_z_lab[i]^2),
      nrow = 2
    ) / n_lab[i]
    
    mu_vec <- c(true_theta[i], true_mu_f[i])
    obs_lab <- MASS::mvrnorm(1, mu = mu_vec, Sigma = Sigma_lab)
    
    Y_bar_vec[i] <- obs_lab[1]
    Z_bar_vec[i] <- obs_lab[2]

    Z_tilde_vec[i] <- rnorm(1,mean = true_mu_f[i],sd = sqrt(var_mean_unl[i]))
  }
  
  b_hat_vec <- Z_bar_vec - Y_bar_vec
  theta_hat_vec <- Z_tilde_vec - b_hat_vec
  lambda_theta_hat_vec <-Y_bar_vec + lambda_vec * (Z_tilde_vec - Z_bar_vec)
  lambda_b_hat_vec <- Z_tilde_vec - lambda_theta_hat_vec
  
  VY <- sim_settings$var_y_lab / n_lab
  VZ <- sim_settings$var_z_lab / n_lab
  VU <- sim_settings$var_z_unlab / n_unl
  C  <- sim_settings$cov_yf_lab / n_lab
  
  calc_lambda_var_theta_hat <-
    VY +
    lambda_vec^2 * VZ +
    lambda_vec^2 * VU -
    2 * lambda_vec * C
  
  calc_lambda_var_b_hat <-
    VY +
    lambda_vec^2 * VZ +
    (1 - lambda_vec)^2 * VU -
    2 * lambda_vec * C
  
  calc_lambda_cov_theta_b <-
    -VY +
    2 * lambda_vec * C -
    lambda_vec^2 * VZ +
    lambda_vec * (1 - lambda_vec) * VU
  
  calc_lambda_corr_theta_b <- ifelse(
    is.finite(calc_lambda_var_theta_hat) &
      is.finite(calc_lambda_var_b_hat) &
      calc_lambda_var_theta_hat > 0 &
      calc_lambda_var_b_hat > 0,
    calc_lambda_cov_theta_b /
      sqrt(calc_lambda_var_theta_hat * calc_lambda_var_b_hat),
    NA_real_
  )
  
  df_sim <- data.frame(
    id = seq_len(M),
    
    n_labeled = n_lab,
    n_unlabeled = n_unl,
    
    Y_bar = Y_bar_vec,
    Z_bar = Z_bar_vec,
    Z_tilde = Z_tilde_vec,
    
    var_Y_bar = sim_settings$var_y_lab / n_lab,
    var_Z_tilde = var_mean_unl,
    
    b_hat = b_hat_vec,
    theta_hat = theta_hat_vec,
    var_b_hat = sim_settings$var_b_hat,
    var_theta_hat = sim_settings$var_theta_hat,
    corr_theta_b = sim_settings$corr_theta_b,
    
    lambda_value = lambda_vec,
    lambda_b_hat = lambda_b_hat_vec,
    lambda_theta_hat = lambda_theta_hat_vec,
    lambda_var_b_hat = calc_lambda_var_b_hat,
    lambda_var_theta_hat = calc_lambda_var_theta_hat,
    lambda_corr_theta_b = calc_lambda_corr_theta_b,
    
    true_theta = true_theta,
    true_theta_f = true_mu_f,
    true_b = true_b,
    stringsAsFactors = FALSE
  )
  
  return(df_sim)
}


get_asymmetric_results_pt_synthetic <- function(df, u_val = 50,
                                                mu_b = -0.1, sd_b = 0.1,
                                                alpha = 0.05) {
  std <- get_input_vectors(df)
  uni <- get_input_vectors(df, "lambda_")
  
  uni_npmle_prior <- GNmix_diffvar(uni$delta, uni$sed,u = u_val,verbose = FALSE)
  uni_param_prior <- estimate_parametric_mle(uni$delta,uni$sed)
  
  
  res_list <- list(
    classic = mapply(cal_ci_normal,df$Y_bar,sqrt(df$var_Y_bar),MoreArgs = list(alpha = alpha)),
    classic_bias = mapply(cal_ci_normal,df$Z_tilde,sqrt(df$var_Z_tilde),MoreArgs = list(alpha = alpha)),
    ppi = mapply(cal_ci_normal,std$beta,std$seb,MoreArgs = list(alpha = alpha)),
    pt = mapply(cal_ci_normal,uni$beta,uni$seb,MoreArgs = list(alpha = alpha)),
    pt_oracle = mapply(cal_ci_parametric,uni$beta,uni$delta,uni$seb,uni$sed,uni$rho,
                       MoreArgs = list(a = mu_b,b = sd_b^2,alpha = alpha)),
    pt_param = mapply(cal_ci_parametric,uni$beta,uni$delta,uni$seb,uni$sed,uni$rho,
                      MoreArgs = list(a = uni_param_prior$hat_mu_b,b = uni_param_prior$hat_var_b,alpha = alpha)),
    pt_npmle = mapply(cal_ci_npmle,uni$beta,uni$delta,uni$seb,uni$sed,uni$rho,
                      MoreArgs = list(v = as.matrix(uni_npmle_prior$u),p = uni_npmle_prior$fu,alpha = alpha)
    )
  )
  
  ci_mats <- lapply(res_list, function(x) {
    m <- t(x)
    colnames(m) <- c("lower", "upper")
    return(m)
  })

  tmp <- future_lapply(1:2, function(k) {
    if (k == 1) {
      cal_ci_ebci(df$Y_bar,sqrt(df$var_Y_bar),alpha = alpha)
    } else {
      cal_ci_ebci(uni$beta,uni$seb,alpha = alpha)
    }
  })
  
  ci_mats$classic_ebci <- tmp[[1]]
  ci_mats$pt_ebci      <- tmp[[2]]
  
  return(list(
    ci_matrices = ci_mats,
    priors = list(uni_npmle = uni_npmle_prior, uni_param = uni_param_prior,
                  oracle = list(mu = mu_b,var = sd_b^2))
  )
  )
}

run_repeated_simulation_synthetic_pt <- function(updated_stat_pt,n_repeats = 200,
                                                 M = 200, u_val = 50,
                                                 mu_theta0 = 4, sd_theta0 = 0.01,
                                                 mu_b = -0.1,sd_b = 0.1, alpha = 0.05,
                                                 seed = NULL, num_workers = 1) {
  stopifnot(length(num_workers) == 1L, is.finite(num_workers), num_workers >= 1)
  num_workers <- as.integer(num_workers)

  run_one <- function(i) {
    df_sim <- simulate_ppi_data_pt(M = M,updated_stat_pt = updated_stat_pt,
                                   mu_theta0 = mu_theta0,sd_theta0 = sd_theta0,
                                   mu_b = mu_b,sd_b = sd_b,
                                   seed = if (is.null(seed)) NULL else seed + i - 1L)
    
    df_sim_clean <- df_sim[
      complete.cases(df_sim) &
        df_sim$var_b_hat > 1e-10 &
        df_sim$var_theta_hat > 1e-10 &
        df_sim$lambda_var_b_hat > 1e-10 &
        df_sim$lambda_var_theta_hat > 1e-10 &
        df_sim$var_Y_bar > 1e-10 &
        df_sim$var_Z_tilde > 1e-10,
    ]
    
    asym_res <- get_asymmetric_results_pt_synthetic(df_sim_clean,u_val = u_val,
                                                    mu_b = mu_b,sd_b = sd_b,alpha = alpha)
    
    theta_true <- df_sim_clean$true_theta

    list(
      coverage = unlist(lapply(asym_res$ci_matrices, cov_rate, theta = theta_true)),
      length = unlist(lapply(asym_res$ci_matrices, function(m) {
        mean(m[, "upper"] - m[, "lower"], na.rm = TRUE)
      })),
      length_ratio = calculate_avg_of_ratios(asym_res)
    )
  }

  if (num_workers == 1L) {
    pb <- txtProgressBar(min = 0, max = n_repeats, style = 3)
    on.exit(close(pb), add = TRUE)
    simulation_results <- lapply(seq_len(n_repeats), function(i) {
      result_i <- run_one(i)
      setTxtProgressBar(pb, i)
      result_i
    })
  } else {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = num_workers)

    old_progressr_enable <- getOption("progressr.enable")
    old_progressr_handlers <- progressr::handlers()
    on.exit({
      options(progressr.enable = old_progressr_enable)
      progressr::handlers(old_progressr_handlers)
    }, add = TRUE)
    options(progressr.enable = TRUE)
    progressr::handlers(progressr::handler_txtprogressbar(clear = FALSE))

    simulation_results <- progressr::with_progress({
      progress <- progressr::progressor(steps = n_repeats)
      progress(amount = 0)
      future.apply::future_lapply(
        seq_len(n_repeats),
        function(i) {
          result_i <- run_one(i)
          progress()
          result_i
        },
        future.seed = NULL,
        future.chunk.size = 1
      )
    })
  }

  trial_cov_asym <- lapply(simulation_results, `[[`, "coverage")
  trial_len_asym <- lapply(simulation_results, `[[`, "length")
  trial_rat_asym <- lapply(simulation_results, `[[`, "length_ratio")

  c_mat <- do.call(rbind, trial_cov_asym)
  l_mat <- do.call(rbind, trial_len_asym)
  r_mat <- do.call(rbind, trial_rat_asym)
  
  mc_se <- function(x) {
    n_eff <- sum(!is.na(x))
    if (n_eff <= 1) return(NA_real_)
    sd(x, na.rm = TRUE) / sqrt(n_eff)
  }
  
  df_asymmetric <- data.frame(
    SE_Coverage = apply(c_mat, 2, mc_se),
    Avg_Coverage = colMeans(c_mat, na.rm = TRUE),
    SE_Length = apply(l_mat, 2, mc_se),
    Avg_Length = colMeans(l_mat, na.rm = TRUE),
    SE_Length_Ratio = apply(r_mat, 2, mc_se),
    Avg_Length_Ratio = colMeans(r_mat, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  return(df_asymmetric)
}



simulate_ppi_data_pt_twopoint <- function(M = 200, updated_stat_pt,
                                          mu_theta0 = 4, sd_theta0 = 0.01,
                                          value = c(0,0.05), probability = c(0.5,0.5),
                                          seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  n_real <- nrow(updated_stat_pt)
  sampled_indices <- sample(seq_len(n_real), size = M, replace = TRUE)
  sim_settings <- updated_stat_pt[sampled_indices, ]

  sig_y_lab <- sqrt(sim_settings$var_y_lab)
  sig_z_lab <- sqrt(sim_settings$var_z_lab)
  cov_lab   <- sim_settings$cov_yf_lab

  var_mean_unl <- sim_settings$var_Z_tilde

  n_lab <- sim_settings$n_labeled
  n_unl <- sim_settings$n_unlabeled

  lambda_vec <- sim_settings$lambda_value
  
  true_theta <- rnorm(M, mean = mu_theta0, sd = sd_theta0)
  true_b     <- sample(value, size = M, replace = TRUE, prob =  probability)
  true_mu_f  <- true_theta + true_b
  
  Y_bar_vec   <- numeric(M)
  Z_bar_vec   <- numeric(M)
  Z_tilde_vec <- numeric(M)
  
  for (i in seq_len(M)) {

    Sigma_lab <- matrix(
      c(sig_y_lab[i]^2, cov_lab[i],
        cov_lab[i],     sig_z_lab[i]^2),
      nrow = 2
    ) / n_lab[i]
    
    mu_vec <- c(true_theta[i], true_mu_f[i])
    
    obs_lab <- MASS::mvrnorm(1, mu = mu_vec, Sigma = Sigma_lab)
    
    Y_bar_vec[i] <- obs_lab[1]
    Z_bar_vec[i] <- obs_lab[2]

    Z_tilde_vec[i] <- rnorm(1,mean = true_mu_f[i],sd = sqrt(var_mean_unl[i]))
  }
  
  b_hat_vec <- Z_bar_vec - Y_bar_vec
  theta_hat_vec <- Z_tilde_vec - b_hat_vec
  lambda_theta_hat_vec <-Y_bar_vec + lambda_vec * (Z_tilde_vec - Z_bar_vec)
  lambda_b_hat_vec <-Z_tilde_vec - lambda_theta_hat_vec
  
  VY <- sim_settings$var_y_lab / n_lab
  VZ <- sim_settings$var_z_lab / n_lab
  VU <- sim_settings$var_z_unlab / n_unl
  C  <- sim_settings$cov_yf_lab / n_lab
  
  calc_lambda_var_theta_hat <-VY +lambda_vec^2 * VZ +lambda_vec^2 * VU -2 * lambda_vec * C
  calc_lambda_var_b_hat <-VY +lambda_vec^2 * VZ +(1 - lambda_vec)^2 * VU -2 * lambda_vec * C
  calc_lambda_cov_theta_b <--VY +2 * lambda_vec * C -lambda_vec^2 * VZ +lambda_vec * (1 - lambda_vec) * VU
  calc_lambda_corr_theta_b <- ifelse(
    is.finite(calc_lambda_var_theta_hat) &
      is.finite(calc_lambda_var_b_hat) &
      calc_lambda_var_theta_hat > 0 &
      calc_lambda_var_b_hat > 0,
    calc_lambda_cov_theta_b /
      sqrt(calc_lambda_var_theta_hat * calc_lambda_var_b_hat),
    NA_real_
  )
  
  df_sim <- data.frame(
    id = seq_len(M),
    
    n_labeled = n_lab,
    n_unlabeled = n_unl,
    
    Y_bar = Y_bar_vec,
    Z_bar = Z_bar_vec,
    Z_tilde = Z_tilde_vec,
    
    var_Y_bar = sim_settings$var_y_lab / n_lab,
    var_Z_tilde = var_mean_unl,

    b_hat = b_hat_vec,
    theta_hat = theta_hat_vec,
    var_b_hat = sim_settings$var_b_hat,
    var_theta_hat = sim_settings$var_theta_hat,
    corr_theta_b = sim_settings$corr_theta_b,

    lambda_value = lambda_vec,
    lambda_b_hat = lambda_b_hat_vec,
    lambda_theta_hat = lambda_theta_hat_vec,
    lambda_var_b_hat = calc_lambda_var_b_hat,
    lambda_var_theta_hat = calc_lambda_var_theta_hat,
    lambda_corr_theta_b = calc_lambda_corr_theta_b,

    true_theta = true_theta,
    true_theta_f = true_mu_f,
    true_b = true_b,
    stringsAsFactors = FALSE
  )
  
  return(df_sim)
}

get_asymmetric_results_pt_synthetic_twopoint <- function(df, u_val = 50,
                                                         value = c(0,0.05), probability = c(0.5,0.5),
                                                         alpha = 0.05) {
  std <- get_input_vectors(df)
  uni <- get_input_vectors(df, "lambda_")
  
  uni_npmle_prior <- GNmix_diffvar(uni$delta, uni$sed,u = u_val,verbose = FALSE)
  uni_param_prior <- estimate_parametric_mle(uni$delta,uni$sed)
  
  
  res_list <- list(
    classic = mapply(cal_ci_normal,df$Y_bar,sqrt(df$var_Y_bar),MoreArgs = list(alpha = alpha)),
    classic_bias = mapply(cal_ci_normal,df$Z_tilde,sqrt(df$var_Z_tilde),MoreArgs = list(alpha = alpha)),
    ppi = mapply(cal_ci_normal,std$beta,std$seb,MoreArgs = list(alpha = alpha)),
    pt = mapply(cal_ci_normal,uni$beta,uni$seb,MoreArgs = list(alpha = alpha)),
    pt_oracle = mapply(cal_ci_npmle,uni$beta,uni$delta,uni$seb,uni$sed,uni$rho,
                       MoreArgs = list(v = value,p = probability,alpha = alpha)),
    pt_param = mapply(cal_ci_parametric,uni$beta,uni$delta,uni$seb,uni$sed,uni$rho,
                      MoreArgs = list(a = uni_param_prior$hat_mu_b,b = uni_param_prior$hat_var_b,alpha = alpha)),
    pt_npmle = mapply(cal_ci_npmle,uni$beta,uni$delta,uni$seb,uni$sed,uni$rho,
                      MoreArgs = list(v = as.matrix(uni_npmle_prior$u),p = uni_npmle_prior$fu,alpha = alpha)
    )
  )
  
  ci_mats <- lapply(res_list, function(x) {
    m <- t(x)
    colnames(m) <- c("lower", "upper")
    return(m)
  })

  tmp <- future_lapply(1:2, function(k) {
    if (k == 1) {
      cal_ci_ebci(df$Y_bar,sqrt(df$var_Y_bar),alpha = alpha)
    } else {
      cal_ci_ebci(uni$beta,uni$seb,alpha = alpha)
    }
  })
  
  ci_mats$classic_ebci <- tmp[[1]]
  ci_mats$pt_ebci      <- tmp[[2]]
  
  return(list(
    ci_matrices = ci_mats,
    priors = list(uni_npmle = uni_npmle_prior, uni_param = uni_param_prior,
                  oracle = list(value = value,probability = probability))
  )
  )
}

run_repeated_simulation_synthetic_pt_twopoint <- function(updated_stat_pt,n_repeats = 200,
                                                          M = 200, u_val = 50,
                                                          mu_theta0 = 4, sd_theta0 = 0.01,
                                                          value = c(0,0.05), probability = c(0.5,0.5),
                                                          alpha = 0.05, seed = NULL,
                                                          num_workers = 1) {
  stopifnot(length(num_workers) == 1L, is.finite(num_workers), num_workers >= 1)
  num_workers <- as.integer(num_workers)

  run_one <- function(i) {
    df_sim <- simulate_ppi_data_pt_twopoint(M = M,updated_stat_pt = updated_stat_pt,
                                            mu_theta0 = mu_theta0,sd_theta0 = sd_theta0,
                                            value = value,probability = probability,
                                            seed = if (is.null(seed)) NULL else seed + i - 1L)
    
    df_sim_clean <- df_sim[
      complete.cases(df_sim) &
        df_sim$var_b_hat > 1e-10 &
        df_sim$var_theta_hat > 1e-10 &
        df_sim$lambda_var_b_hat > 1e-10 &
        df_sim$lambda_var_theta_hat > 1e-10 &
        df_sim$var_Y_bar > 1e-10 &
        df_sim$var_Z_tilde > 1e-10,
    ]
    
    asym_res <- get_asymmetric_results_pt_synthetic_twopoint(df_sim_clean,u_val = u_val,
                                                             value = value,probability = probability,alpha = alpha)
    
    theta_true <- df_sim_clean$true_theta

    list(
      coverage = unlist(lapply(asym_res$ci_matrices, cov_rate, theta = theta_true)),
      length = unlist(lapply(asym_res$ci_matrices, function(m) {
        mean(m[, "upper"] - m[, "lower"], na.rm = TRUE)
      })),
      length_ratio = calculate_avg_of_ratios(asym_res)
    )
  }

  if (num_workers == 1L) {
    pb <- txtProgressBar(min = 0, max = n_repeats, style = 3)
    on.exit(close(pb), add = TRUE)
    simulation_results <- lapply(seq_len(n_repeats), function(i) {
      result_i <- run_one(i)
      setTxtProgressBar(pb, i)
      result_i
    })
  } else {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = num_workers)

    old_progressr_enable <- getOption("progressr.enable")
    old_progressr_handlers <- progressr::handlers()
    on.exit({
      options(progressr.enable = old_progressr_enable)
      progressr::handlers(old_progressr_handlers)
    }, add = TRUE)
    options(progressr.enable = TRUE)
    progressr::handlers(progressr::handler_txtprogressbar(clear = FALSE))

    simulation_results <- progressr::with_progress({
      progress <- progressr::progressor(steps = n_repeats)
      progress(amount = 0)
      future.apply::future_lapply(
        seq_len(n_repeats),
        function(i) {
          result_i <- run_one(i)
          progress()
          result_i
        },
        future.seed = NULL,
        future.chunk.size = 1
      )
    })
  }

  trial_cov_asym <- lapply(simulation_results, `[[`, "coverage")
  trial_len_asym <- lapply(simulation_results, `[[`, "length")
  trial_rat_asym <- lapply(simulation_results, `[[`, "length_ratio")

  c_mat <- do.call(rbind, trial_cov_asym)
  l_mat <- do.call(rbind, trial_len_asym)
  r_mat <- do.call(rbind, trial_rat_asym)
  
  mc_se <- function(x) {
    n_eff <- sum(!is.na(x))
    if (n_eff <= 1) return(NA_real_)
    sd(x, na.rm = TRUE) / sqrt(n_eff)
  }
  
  df_asymmetric <- data.frame(
    SE_Coverage = apply(c_mat, 2, mc_se),
    Avg_Coverage = colMeans(c_mat, na.rm = TRUE),
    SE_Length = apply(l_mat, 2, mc_se),
    Avg_Length = colMeans(l_mat, na.rm = TRUE),
    SE_Length_Ratio = apply(r_mat, 2, mc_se),
    Avg_Length_Ratio = colMeans(r_mat, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  return(df_asymmetric)
}

# ---- Diagnostic data ----
make_amazon_diagnostic <- function(file, r = 0.2, summary_seed = 5, pseudo_seed = 123) {
  stopifnot(file.exists(file), r > 0, r < 1)
  sum_stat <- load_data_amazon(file, r = r, seed = summary_seed)
  if (!is.null(pseudo_seed)) set.seed(pseudo_seed)
  ls_tbl <- rhdf5::h5ls(file, recursive = FALSE)
  keys <- ls_tbl$name[ls_tbl$otype == "H5I_DATASET"]
  if (length(keys) == 0L) stop("No root-level datasets found.")
  pseudo_rows <- lapply(keys, function(k) {
    m <- rhdf5::h5read(file, paste0("/", k))
    if (!is.matrix(m) || ncol(m) < 2L) stop(sprintf("Dataset '%s' is not a two-column matrix.", k))
    Y <- as.numeric(m[, 1])
    Z <- as.numeric(m[, 2])
    n <- length(Y)
    if (n < 2L) stop(sprintf("Dataset '%s' has fewer than two observations.", k))
    idx <- sample.int(n)
    n_bias <- floor(n / 2)
    idx_bias <- idx[seq_len(n_bias)]
    idx_theta <- idx[(n_bias + 1L):n]
    data.frame(
      group = k,
      theta = mean(Y[idx_theta]),
      bias = mean(Z[idx_bias] - Y[idx_bias]),
      stringsAsFactors = FALSE
    )
  })
  pseudo_stat <- do.call(rbind, pseudo_rows)
  rownames(pseudo_stat) <- NULL
  matched <- match(pseudo_stat$group, sum_stat$group)
  if (anyNA(matched)) stop("Task IDs do not match between pseudo and summary splits.")
  diagnostic <- data.frame(
    theta = pseudo_stat$theta,
    bias = pseudo_stat$bias,
    bias_hat = sum_stat$b_hat[matched],
    rho = sum_stat$corr_theta_b[matched],
    tau2 = sum_stat$var_b_hat[matched],
    sigma2 = sum_stat$var_theta_hat[matched]
  )
  rownames(diagnostic) <- NULL
  diagnostic
}
