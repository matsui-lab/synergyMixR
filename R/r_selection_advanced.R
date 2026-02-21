#' Ensure MPCA model has W field (reconstruct from P and D if needed)
#'
#' @param mpca_fit A fitted MPCA model
#' @return The same model with W field guaranteed to exist
#' @keywords internal
ensure_mpca_has_W <- function(mpca_fit) {
  if (is.null(mpca_fit$W) && !is.null(mpca_fit$P) && !is.null(mpca_fit$D)) {
    K <- length(mpca_fit$P)
    mpca_fit$W <- vector("list", K)
    for (k in seq_len(K)) {
      mpca_fit$W[[k]] <- mpca_fit$P[[k]] %*% mpca_fit$D[[k]]
    }
  }
  return(mpca_fit)
}


#' Compute Smoothness Penalty for Factor Scores
#'
#' Computes the temporal smoothness penalty for factor scores using second-order
#' differences. This penalty measures how "wiggly" the factor activation profiles are,
#' with higher values indicating more high-frequency noise.
#'
#' @param factor_scores A list of N matrices, each (T_i x r), containing factor scores
#'   for each subject.
#'
#' @return A numeric scalar representing the total smoothness penalty across all
#'   subjects and factors.
#'
#' @details
#' The smoothness penalty is computed as:
#' \deqn{
#'   \text{smoothness} = \sum_{i=1}^{N} \sum_{j=1}^{r} \sum_{t=3}^{T_i} (z_{i,j,t} - 2z_{i,j,t-1} + z_{i,j,t-2})^2
#' }
#' where \eqn{z_{i,j,t}} is the j-th factor score for subject i at time t.
#'
#' This is equivalent to summing the squared second-order differences across all
#' time series. When r is too large, models tend to fit high-frequency noise,
#' producing wiggly factors with high smoothness penalties.
#'
#' @examples
#' \dontrun{
#' # Compute factor scores first
#' mfa_fit <- compute_factor_scores_mfa(mfa_fit, list_of_data)
#' # Then compute smoothness
#' smooth_val <- compute_smoothness_penalty(mfa_fit$factor_scores)
#' }
#'
#' @export
compute_smoothness_penalty <- function(factor_scores) {
  if (!is.list(factor_scores)) {
    stop("factor_scores must be a list of matrices")
  }
  
  N <- length(factor_scores)
  if (N == 0) {
    return(0)
  }
  
  total_smoothness <- 0
  
  for (i in seq_len(N)) {
    Z_i <- factor_scores[[i]]
    if (is.null(Z_i) || !is.matrix(Z_i)) {
      next
    }
    
    T_i <- nrow(Z_i)
    r <- ncol(Z_i)
    
    # Need at least 3 time points for second-order differences
    if (T_i < 3) {
      next
    }
    
    # Compute second-order differences for each factor
    for (j in seq_len(r)) {
      z_j <- Z_i[, j]
      # diff(z_j, differences = 2) computes second-order differences
      second_diff <- diff(z_j, differences = 2)
      total_smoothness <- total_smoothness + sum(second_diff^2)
    }
  }
  
  return(total_smoothness)
}


#' Compute VSS (Very Simple Structure) Criterion for Factor Loadings
#'
#' Computes the Very Simple Structure (VSS) criterion, which evaluates how well
#' factor loadings exhibit simple structure (each variable loads strongly on one
#' factor, weakly on others).
#'
#' @param Lambda A list of K loading matrices, each (M x r), representing factor
#'   loadings for each cluster.
#' @param threshold Threshold for determining "strong" vs "weak" loadings.
#'   Default is 0.5. Loadings with absolute value above this threshold are
#'   considered strong.
#' @param method Method for computing VSS. Options: "threshold" (count variables with 
#'   exactly one strong loading above threshold), "ratio" (compute concentration ratio 
#'   max(lambda^2)/sum(lambda^2) for each variable). Default is "threshold".
#' @param rotate Rotation method to apply before computing VSS. Options: "none" (no rotation),
#'   "varimax" (orthogonal rotation). Default is "none". VSS on unrotated FA loadings is
#'   rotation-dependent and often uninformative. Rotation is recommended for meaningful VSS.
#'
#' @return A numeric scalar representing the VSS criterion. Higher values indicate
#'   better simple structure. Returns NA if rotation fails or loadings are degenerate.
#'
#' @details
#' The VSS criterion evaluates simple structure in factor loadings. Two methods are available:
#' 
#' **Threshold method (default):** For each variable (muscle), count how many factors have 
#' |loading| > threshold. Simple structure is achieved when this count is exactly 1. 
#' The VSS score is the proportion of variables across all clusters that exhibit simple structure.
#' 
#' **Ratio method:** For each variable (muscle), compute the concentration ratio 
#' max(lambda^2)/sum(lambda^2), which measures how concentrated the loadings are on a single factor.
#' The VSS score is the average concentration ratio across all variables. This method is more 
#' robust when loadings are small in absolute magnitude, as it doesn't depend on an arbitrary threshold.
#'
#' This implementation is inspired by the psych::VSS() function, adapted for mixture models 
#' where each cluster has its own loading matrix.
#'
#' @references
#' Revelle, W., & Rocklin, T. (1979). Very Simple Structure: An Alternative
#' Procedure For Estimating The Optimal Number Of Interpretable Factors.
#' Multivariate Behavioral Research, 14(4), 403-414.
#'
#' @examples
#' \dontrun{
#' # For an MFA model
#' vss_val <- compute_vss_criterion(mfa_fit$Lambda, threshold = 0.5)
#' }
#'
#' @export
compute_vss_criterion <- function(Lambda, threshold = 0.5, method = c("threshold", "ratio"), rotate = c("none", "varimax")) {
  method <- match.arg(method)
  rotate <- match.arg(rotate)
  
  if (!is.list(Lambda)) {
    stop("Lambda must be a list of loading matrices")
  }
  
  K <- length(Lambda)
  if (K == 0) {
    return(0)
  }
  
  # Apply rotation if requested
  if (rotate != "none") {
    Lambda_rotated <- vector("list", K)
    for (k in seq_len(K)) {
      Lambda_k <- Lambda[[k]]
      if (is.null(Lambda_k) || !is.matrix(Lambda_k)) {
        Lambda_rotated[[k]] <- Lambda_k
        next
      }
      
      r <- ncol(Lambda_k)
      if (r < 2) {
        # Cannot rotate with fewer than 2 factors
        Lambda_rotated[[k]] <- Lambda_k
        next
      }
      
      # Apply varimax rotation
      tryCatch({
        vout <- stats::varimax(Lambda_k)
        Lambda_rotated[[k]] <- as.matrix(unclass(vout$loadings))
      }, error = function(e) {
        warning(sprintf("Varimax rotation failed for cluster %d: %s. Using unrotated loadings.", k, e$message))
        Lambda_rotated[[k]] <- Lambda_k
      })
    }
    Lambda <- Lambda_rotated
  }
  
  if (method == "threshold") {
    # Original threshold-based method
    total_vars <- 0
    simple_structure_count <- 0
    
    for (k in seq_len(K)) {
      Lambda_k <- Lambda[[k]]
      if (is.null(Lambda_k) || !is.matrix(Lambda_k)) {
        next
      }
      
      M <- nrow(Lambda_k)
      r <- ncol(Lambda_k)
      
      if (r < 2) {
        # With only 1 factor, all variables trivially have simple structure
        simple_structure_count <- simple_structure_count + M
        total_vars <- total_vars + M
        next
      }
      
      # For each variable (muscle), count how many factors it loads strongly on
      for (m in seq_len(M)) {
        loadings_m <- abs(Lambda_k[m, ])
        strong_loadings <- sum(loadings_m > threshold)
        
        # Simple structure: exactly one strong loading
        if (strong_loadings == 1) {
          simple_structure_count <- simple_structure_count + 1
        }
        
        total_vars <- total_vars + 1
      }
    }
    
    if (total_vars == 0) {
      return(0)
    }
    
    # Return proportion of variables with simple structure
    vss_score <- simple_structure_count / total_vars
    return(vss_score)
    
  } else if (method == "ratio") {
    # Ratio-based method: for each muscle, compute max(lambda_mj^2) / sum(lambda_mj^2)
    # Higher ratio means more concentrated loading (better simple structure)
    total_vars <- 0
    sum_ratios <- 0
    
    for (k in seq_len(K)) {
      Lambda_k <- Lambda[[k]]
      if (is.null(Lambda_k) || !is.matrix(Lambda_k)) {
        next
      }
      
      M <- nrow(Lambda_k)
      r <- ncol(Lambda_k)
      
      if (r < 1) {
        next
      }
      
      # For each variable (muscle), compute concentration ratio
      for (m in seq_len(M)) {
        loadings_m_sq <- Lambda_k[m, ]^2  # Squared loadings (sign-invariant)
        sum_sq <- sum(loadings_m_sq)
        
        if (sum_sq > 0) {
          ratio_m <- max(loadings_m_sq) / sum_sq
          sum_ratios <- sum_ratios + ratio_m
          total_vars <- total_vars + 1
        }
      }
    }
    
    if (total_vars == 0) {
      return(NA_real_)
    }
    
    # Return average concentration ratio across all variables
    vss_score <- sum_ratios / total_vars
    return(vss_score)
  }
}


#' Select Optimal r by Smoothness Penalty for MFA
#'
#' Selects the optimal number of factors (r) for a given number of clusters (K)
#' by minimizing the smoothness penalty of factor activation profiles, subject
#' to constraints on VAF or BIC.
#'
#' @param list_of_data A list of N matrices, each (T_i x M).
#' @param K Number of clusters (fixed).
#' @param rvec Vector of candidate r values to evaluate.
#' @param lambda_smooth Smoothness penalty weight. Default is 1e-3.
#' @param vaf_threshold Minimum VAF threshold. Models with VAF below this are excluded.
#'   Default is 0.9.
#' @param max_iter Maximum EM iterations.
#' @param nIterFA Sub-iterations for FA update.
#' @param tol Convergence tolerance.
#' @param n_init Number of random initializations.
#' @param use_kmeans_init Whether to use k-means initialization.
#' @param subject_rdim_for_kmeans Dimension for k-means features.
#' @param mc_cores Number of cores for parallelization.
#' @param n_threads Number of OpenMP threads.
#'
#' @return A list with:
#'   \item{best_r}{The selected r value.}
#'   \item{summary}{A data.frame with columns r, BIC, VAF, smoothness, penalized_score.}
#'   \item{best_model}{The fitted model for the selected r.}
#'
#' @details
#' This function fits MFA models for each r in rvec with fixed K, computes
#' factor scores, and evaluates the smoothness penalty. The optimal r is selected
#' by minimizing:
#' \deqn{
#'   \text{score} = \text{BIC} + \lambda_{\text{smooth}} \times \text{smoothness}
#' }
#' among models satisfying VAF >= vaf_threshold.
#'
#' @examples
#' \dontrun{
#' result <- select_r_by_smoothness_mfa(
#'   list_of_data, K = 3, rvec = 1:5,
#'   lambda_smooth = 1e-3, vaf_threshold = 0.9
#' )
#' print(result$summary)
#' best_model <- result$best_model
#' }
#'
#' @export
select_r_by_smoothness_mfa <- function(
    list_of_data,
    K,
    rvec = 1:6,
    lambda_smooth = 1e-3,
    vaf_threshold = 0.9,
    max_iter = 50,
    nIterFA = 5,
    tol = 1e-3,
    n_init = 1,
    use_kmeans_init = TRUE,
    subject_rdim_for_kmeans = 2,
    mc_cores = 1,
    n_threads = 1
) {
  N <- length(list_of_data)
  N_total_rows <- sum(sapply(list_of_data, nrow))
  M <- ncol(list_of_data[[1]])
  
  results_list <- list()
  
  for (r in rvec) {
    # Fit MFA model
    fit_mfa <- mfa_em_fit(
      list_of_data = list_of_data,
      K = K,
      r = r,
      max_iter = max_iter,
      nIterFA = nIterFA,
      tol = tol,
      n_init = n_init,
      use_kmeans_init = use_kmeans_init,
      subject_rdim_for_kmeans = subject_rdim_for_kmeans,
      mc_cores = mc_cores,
      n_threads = n_threads
    )
    
    # Compute BIC
    loglik_val <- fit_mfa$logLik
    bic_val <- compute_BIC_mfa(loglik_val, K, r, M, N_total_rows)
    
    # Compute VAF
    vaf_val <- compute_global_vaf_mfa(list_of_data, fit_mfa)
    
    # Compute factor scores
    fit_mfa <- compute_factor_scores_mfa(fit_mfa, list_of_data)
    
    # Compute smoothness penalty
    smoothness_val <- compute_smoothness_penalty(fit_mfa$factor_scores)
    
    # Penalized score
    penalized_score <- bic_val + lambda_smooth * smoothness_val
    
    results_list[[length(results_list) + 1]] <- list(
      r = r,
      model = fit_mfa,
      logLik = loglik_val,
      BIC = bic_val,
      VAF = vaf_val,
      smoothness = smoothness_val,
      penalized_score = penalized_score
    )
  }
  
  # Build summary data frame
  df_summary <- do.call(rbind, lapply(results_list, function(res) {
    data.frame(
      r = res$r,
      logLik = res$logLik,
      BIC = res$BIC,
      VAF = res$VAF,
      smoothness = res$smoothness,
      penalized_score = res$penalized_score
    )
  }))
  
  # Filter by VAF threshold
  df_valid <- df_summary[df_summary$VAF >= vaf_threshold, ]
  
  if (nrow(df_valid) == 0) {
    warning("No models satisfy VAF threshold. Returning model with highest VAF.")
    best_idx <- which.max(df_summary$VAF)
  } else {
    # Select r with minimum penalized score among valid models
    best_idx_in_valid <- which.min(df_valid$penalized_score)
    best_r <- df_valid$r[best_idx_in_valid]
    best_idx <- which(df_summary$r == best_r)[1]
  }
  
  best_model_info <- results_list[[best_idx]]
  
  cat("=== Best r by smoothness penalty (MFA) ===\n")
  cat(sprintf("Selected r = %d (K = %d)\n", best_model_info$r, K))
  cat(sprintf("BIC = %.2f, VAF = %.4f, Smoothness = %.2f, Penalized Score = %.2f\n",
              best_model_info$BIC, best_model_info$VAF,
              best_model_info$smoothness, best_model_info$penalized_score))
  
  list(
    best_r = best_model_info$r,
    summary = df_summary,
    best_model = best_model_info$model
  )
}


#' Select Optimal r by Smoothness Penalty for MPCA
#'
#' Selects the optimal number of factors (r) for a given number of clusters (K)
#' by minimizing the smoothness penalty of PC activation profiles, subject
#' to constraints on VAF or BIC.
#'
#' @param list_of_data A list of N matrices, each (T_i x M).
#' @param K Number of clusters (fixed).
#' @param rvec Vector of candidate r values to evaluate.
#' @param lambda_smooth Smoothness penalty weight. Default is 1e-3.
#' @param vaf_threshold Minimum VAF threshold. Models with VAF below this are excluded.
#'   Default is 0.9.
#' @param max_iter Maximum EM iterations.
#' @param nIterPCA Sub-iterations for PCA update.
#' @param tol Convergence tolerance.
#' @param method "EM" or "closed_form".
#' @param n_init Number of random initializations.
#' @param use_kmeans_init Whether to use k-means initialization.
#' @param subject_rdim_for_kmeans Dimension for k-means features.
#' @param mc_cores Number of cores for parallelization.
#' @param n_threads Number of OpenMP threads.
#'
#' @return A list with:
#'   \item{best_r}{The selected r value.}
#'   \item{summary}{A data.frame with columns r, BIC, VAF, smoothness, penalized_score.}
#'   \item{best_model}{The fitted model for the selected r.}
#'
#' @details
#' This function fits MPCA models for each r in rvec with fixed K, computes
#' PC scores, and evaluates the smoothness penalty. The optimal r is selected
#' by minimizing:
#' \deqn{
#'   \text{score} = \text{BIC} + \lambda_{\text{smooth}} \times \text{smoothness}
#' }
#' among models satisfying VAF >= vaf_threshold.
#'
#' @examples
#' \dontrun{
#' result <- select_r_by_smoothness_mpca(
#'   list_of_data, K = 3, rvec = 1:5,
#'   lambda_smooth = 1e-3, vaf_threshold = 0.9
#' )
#' print(result$summary)
#' best_model <- result$best_model
#' }
#'
#' @export
select_r_by_smoothness_mpca <- function(
    list_of_data,
    K,
    rvec = 1:6,
    lambda_smooth = 1e-3,
    vaf_threshold = 0.9,
    max_iter = 50,
    nIterPCA = 5,
    tol = 1e-3,
    method = "EM",
    n_init = 1,
    use_kmeans_init = TRUE,
    subject_rdim_for_kmeans = 2,
    mc_cores = 1,
    n_threads = 1
) {
  N <- length(list_of_data)
  N_total_rows <- sum(sapply(list_of_data, nrow))
  M <- ncol(list_of_data[[1]])
  
  results_list <- list()
  
  for (r in rvec) {
    # Fit MPCA model
    fit_mpca <- mixture_pca_em_fit(
      list_of_data = list_of_data,
      K = K,
      r = r,
      max_iter = max_iter,
      nIterPCA = nIterPCA,
      tol = tol,
      method = method,
      n_init = n_init,
      use_kmeans_init = use_kmeans_init,
      subject_rdim_for_kmeans = subject_rdim_for_kmeans,
      mc_cores = mc_cores,
      n_threads = n_threads
    )
    
    # Ensure W field exists (reconstruct from P and D if needed)
    fit_mpca <- ensure_mpca_has_W(fit_mpca)
    
    # Compute BIC
    loglik_val <- fit_mpca$logLik
    bic_val <- compute_BIC_mpca(loglik_val, K, r, M, N_total_rows)
    
    # Compute VAF
    vaf_val <- compute_global_vaf_mpca(list_of_data, fit_mpca)
    
    # Compute factor scores
    fit_mpca <- compute_factor_scores_mpca(fit_mpca, list_of_data)
    
    # Compute smoothness penalty
    smoothness_val <- compute_smoothness_penalty(fit_mpca$factor_scores)
    
    # Penalized score
    penalized_score <- bic_val + lambda_smooth * smoothness_val
    
    results_list[[length(results_list) + 1]] <- list(
      r = r,
      model = fit_mpca,
      logLik = loglik_val,
      BIC = bic_val,
      VAF = vaf_val,
      smoothness = smoothness_val,
      penalized_score = penalized_score
    )
  }
  
  # Build summary data frame
  df_summary <- do.call(rbind, lapply(results_list, function(res) {
    data.frame(
      r = res$r,
      logLik = res$logLik,
      BIC = res$BIC,
      VAF = res$VAF,
      smoothness = res$smoothness,
      penalized_score = res$penalized_score
    )
  }))
  
  # Filter by VAF threshold
  df_valid <- df_summary[df_summary$VAF >= vaf_threshold, ]
  
  if (nrow(df_valid) == 0) {
    warning("No models satisfy VAF threshold. Returning model with highest VAF.")
    best_idx <- which.max(df_summary$VAF)
  } else {
    # Select r with minimum penalized score among valid models
    best_idx_in_valid <- which.min(df_valid$penalized_score)
    best_r <- df_valid$r[best_idx_in_valid]
    best_idx <- which(df_summary$r == best_r)[1]
  }
  
  best_model_info <- results_list[[best_idx]]
  
  cat("=== Best r by smoothness penalty (MPCA) ===\n")
  cat(sprintf("Selected r = %d (K = %d)\n", best_model_info$r, K))
  cat(sprintf("BIC = %.2f, VAF = %.4f, Smoothness = %.2f, Penalized Score = %.2f\n",
              best_model_info$BIC, best_model_info$VAF,
              best_model_info$smoothness, best_model_info$penalized_score))
  
  list(
    best_r = best_model_info$r,
    summary = df_summary,
    best_model = best_model_info$model
  )
}


#' Select Optimal r by VSS Criterion for MFA
#'
#' Selects the optimal number of factors (r) for a given number of clusters (K)
#' by maximizing the VSS (Very Simple Structure) criterion, subject to constraints
#' on VAF or BIC.
#'
#' @param list_of_data A list of N matrices, each (T_i x M).
#' @param K Number of clusters (fixed).
#' @param rvec Vector of candidate r values to evaluate.
#' @param vss_threshold Threshold for determining strong vs weak loadings. Default is 0.5.
#' @param vaf_threshold Minimum VAF threshold. Models with VAF below this are excluded.
#'   Default is 0.9.
#' @param max_iter Maximum EM iterations.
#' @param nIterFA Sub-iterations for FA update.
#' @param tol Convergence tolerance.
#' @param n_init Number of random initializations.
#' @param use_kmeans_init Whether to use k-means initialization.
#' @param subject_rdim_for_kmeans Dimension for k-means features.
#' @param mc_cores Number of cores for parallelization.
#' @param n_threads Number of OpenMP threads.
#'
#' @return A list with:
#'   \item{best_r}{The selected r value.}
#'   \item{summary}{A data.frame with columns r, BIC, VAF, VSS.}
#'   \item{best_model}{The fitted model for the selected r.}
#'
#' @details
#' This function fits MFA models for each r in rvec with fixed K, and evaluates
#' the VSS criterion. The optimal r is selected by maximizing VSS among models
#' satisfying VAF >= vaf_threshold.
#'
#' Higher VSS values indicate better simple structure, where each muscle loads
#' strongly on exactly one synergy.
#'
#' @examples
#' \dontrun{
#' result <- select_r_by_vss_mfa(
#'   list_of_data, K = 3, rvec = 1:5,
#'   vss_threshold = 0.5, vaf_threshold = 0.9
#' )
#' print(result$summary)
#' best_model <- result$best_model
#' }
#'
#' @export
select_r_by_vss_mfa <- function(
    list_of_data,
    K,
    rvec = 1:6,
    vss_threshold = 0.5,
    vaf_threshold = 0.9,
    max_iter = 50,
    nIterFA = 5,
    tol = 1e-3,
    n_init = 1,
    use_kmeans_init = TRUE,
    subject_rdim_for_kmeans = 2,
    mc_cores = 1,
    n_threads = 1
) {
  N <- length(list_of_data)
  N_total_rows <- sum(sapply(list_of_data, nrow))
  M <- ncol(list_of_data[[1]])
  
  results_list <- list()
  
  for (r in rvec) {
    # Fit MFA model
    fit_mfa <- mfa_em_fit(
      list_of_data = list_of_data,
      K = K,
      r = r,
      max_iter = max_iter,
      nIterFA = nIterFA,
      tol = tol,
      n_init = n_init,
      use_kmeans_init = use_kmeans_init,
      subject_rdim_for_kmeans = subject_rdim_for_kmeans,
      mc_cores = mc_cores,
      n_threads = n_threads
    )
    
    # Compute BIC
    loglik_val <- fit_mfa$logLik
    bic_val <- compute_BIC_mfa(loglik_val, K, r, M, N_total_rows)
    
    # Compute VAF
    vaf_val <- compute_global_vaf_mfa(list_of_data, fit_mfa)
    
    # Compute VSS criterion
    vss_val <- compute_vss_criterion(fit_mfa$Lambda, threshold = vss_threshold)
    
    results_list[[length(results_list) + 1]] <- list(
      r = r,
      model = fit_mfa,
      logLik = loglik_val,
      BIC = bic_val,
      VAF = vaf_val,
      VSS = vss_val
    )
  }
  
  # Build summary data frame
  df_summary <- do.call(rbind, lapply(results_list, function(res) {
    data.frame(
      r = res$r,
      logLik = res$logLik,
      BIC = res$BIC,
      VAF = res$VAF,
      VSS = res$VSS
    )
  }))
  
  # Filter by VAF threshold
  df_valid <- df_summary[df_summary$VAF >= vaf_threshold, ]
  
  if (nrow(df_valid) == 0) {
    warning("No models satisfy VAF threshold. Returning model with highest VAF.")
    best_idx <- which.max(df_summary$VAF)
  } else {
    # Select r with maximum VSS among valid models
    best_idx_in_valid <- which.max(df_valid$VSS)
    best_r <- df_valid$r[best_idx_in_valid]
    best_idx <- which(df_summary$r == best_r)[1]
  }
  
  best_model_info <- results_list[[best_idx]]
  
  cat("=== Best r by VSS criterion (MFA) ===\n")
  cat(sprintf("Selected r = %d (K = %d)\n", best_model_info$r, K))
  cat(sprintf("BIC = %.2f, VAF = %.4f, VSS = %.4f\n",
              best_model_info$BIC, best_model_info$VAF, best_model_info$VSS))
  
  list(
    best_r = best_model_info$r,
    summary = df_summary,
    best_model = best_model_info$model
  )
}


#' Select Optimal r by VSS Criterion for MPCA
#'
#' Selects the optimal number of factors (r) for a given number of clusters (K)
#' by maximizing the VSS (Very Simple Structure) criterion, subject to constraints
#' on VAF or BIC.
#'
#' @param list_of_data A list of N matrices, each (T_i x M).
#' @param K Number of clusters (fixed).
#' @param rvec Vector of candidate r values to evaluate.
#' @param vss_threshold Threshold for determining strong vs weak loadings. Default is 0.5.
#' @param vaf_threshold Minimum VAF threshold. Models with VAF below this are excluded.
#'   Default is 0.9.
#' @param max_iter Maximum EM iterations.
#' @param nIterPCA Sub-iterations for PCA update.
#' @param tol Convergence tolerance.
#' @param method "EM" or "closed_form".
#' @param n_init Number of random initializations.
#' @param use_kmeans_init Whether to use k-means initialization.
#' @param subject_rdim_for_kmeans Dimension for k-means features.
#' @param mc_cores Number of cores for parallelization.
#' @param n_threads Number of OpenMP threads.
#'
#' @return A list with:
#'   \item{best_r}{The selected r value.}
#'   \item{summary}{A data.frame with columns r, BIC, VAF, VSS.}
#'   \item{best_model}{The fitted model for the selected r.}
#'
#' @details
#' This function fits MPCA models for each r in rvec with fixed K, and evaluates
#' the VSS criterion. The optimal r is selected by maximizing VSS among models
#' satisfying VAF >= vaf_threshold.
#'
#' Higher VSS values indicate better simple structure, where each muscle loads
#' strongly on exactly one synergy.
#'
#' @examples
#' \dontrun{
#' result <- select_r_by_vss_mpca(
#'   list_of_data, K = 3, rvec = 1:5,
#'   vss_threshold = 0.5, vaf_threshold = 0.9
#' )
#' print(result$summary)
#' best_model <- result$best_model
#' }
#'
#' @export
select_r_by_vss_mpca <- function(
    list_of_data,
    K,
    rvec = 1:6,
    vss_threshold = 0.5,
    vaf_threshold = 0.9,
    max_iter = 50,
    nIterPCA = 5,
    tol = 1e-3,
    method = "EM",
    n_init = 1,
    use_kmeans_init = TRUE,
    subject_rdim_for_kmeans = 2,
    mc_cores = 1,
    n_threads = 1
) {
  N <- length(list_of_data)
  N_total_rows <- sum(sapply(list_of_data, nrow))
  M <- ncol(list_of_data[[1]])
  
  results_list <- list()
  
  for (r in rvec) {
    # Fit MPCA model
    fit_mpca <- mixture_pca_em_fit(
      list_of_data = list_of_data,
      K = K,
      r = r,
      max_iter = max_iter,
      nIterPCA = nIterPCA,
      tol = tol,
      method = method,
      n_init = n_init,
      use_kmeans_init = use_kmeans_init,
      subject_rdim_for_kmeans = subject_rdim_for_kmeans,
      mc_cores = mc_cores,
      n_threads = n_threads
    )
    
    # Ensure W field exists (reconstruct from P and D if needed)
    fit_mpca <- ensure_mpca_has_W(fit_mpca)
    
    # Compute BIC
    loglik_val <- fit_mpca$logLik
    bic_val <- compute_BIC_mpca(loglik_val, K, r, M, N_total_rows)
    
    # Compute VAF
    vaf_val <- compute_global_vaf_mpca(list_of_data, fit_mpca)
    
    # Compute VSS criterion
    vss_val <- compute_vss_criterion(fit_mpca$W, threshold = vss_threshold)
    
    results_list[[length(results_list) + 1]] <- list(
      r = r,
      model = fit_mpca,
      logLik = loglik_val,
      BIC = bic_val,
      VAF = vaf_val,
      VSS = vss_val
    )
  }
  
  # Build summary data frame
  df_summary <- do.call(rbind, lapply(results_list, function(res) {
    data.frame(
      r = res$r,
      logLik = res$logLik,
      BIC = res$BIC,
      VAF = res$VAF,
      VSS = res$VSS
    )
  }))
  
  # Filter by VAF threshold
  df_valid <- df_summary[df_summary$VAF >= vaf_threshold, ]
  
  if (nrow(df_valid) == 0) {
    warning("No models satisfy VAF threshold. Returning model with highest VAF.")
    best_idx <- which.max(df_summary$VAF)
  } else {
    # Select r with maximum VSS among valid models
    best_idx_in_valid <- which.max(df_valid$VSS)
    best_r <- df_valid$r[best_idx_in_valid]
    best_idx <- which(df_summary$r == best_r)[1]
  }
  
  best_model_info <- results_list[[best_idx]]
  
  cat("=== Best r by VSS criterion (MPCA) ===\n")
  cat(sprintf("Selected r = %d (K = %d)\n", best_model_info$r, K))
  cat(sprintf("BIC = %.2f, VAF = %.4f, VSS = %.4f\n",
              best_model_info$BIC, best_model_info$VAF, best_model_info$VSS))
  
  list(
    best_r = best_model_info$r,
    summary = df_summary,
    best_model = best_model_info$model
  )
}


#' Select Optimal r Using Multi-Criteria Approach for MFA
#'
#' Selects the optimal number of factors (r) for a given number of clusters (K)
#' by combining multiple criteria: BIC, VAF, smoothness penalty, and VSS.
#'
#' @param list_of_data A list of N matrices, each (T_i x M).
#' @param K Number of clusters (fixed).
#' @param rvec Vector of candidate r values to evaluate. Default is 1:6.
#' @param vaf_threshold Minimum VAF threshold. Models with VAF below this are excluded.
#'   Default is 0.9.
#' @param lambda_smooth Smoothness penalty weight. Default is 1e-3.
#' @param vss_threshold Threshold for VSS criterion. Default is 0.5.
#' @param weights Named vector of weights for each criterion. Default is
#'   c(BIC = 1, VAF = 1, smoothness = 1, VSS = 1).
#' @param max_iter Maximum EM iterations.
#' @param nIterFA Sub-iterations for FA update.
#' @param tol Convergence tolerance.
#' @param n_init Number of random initializations.
#' @param use_kmeans_init Whether to use k-means initialization.
#' @param subject_rdim_for_kmeans Dimension for k-means features.
#' @param mc_cores Number of cores for parallelization.
#' @param n_threads Number of OpenMP threads.
#'
#' @return A list with:
#'   \item{best_r}{The selected r value.}
#'   \item{summary}{A data.frame with columns r, BIC, VAF, smoothness, VSS, composite_score.}
#'   \item{best_model}{The fitted model for the selected r.}
#'
#' @details
#' This function combines multiple criteria to select r:
#' \itemize{
#'   \item BIC: Lower is better (penalizes complexity)
#'   \item VAF: Higher is better (rewards fit quality)
#'   \item Smoothness: Lower is better (penalizes wiggly factors)
#'   \item VSS: Higher is better (rewards simple structure)
#' }
#'
#' Each criterion is normalized to [0, 1] range, then combined using weights:
#' \deqn{
#'   \text{score} = w_{\text{BIC}} \times \text{BIC}_{\text{norm}} +
#'                  w_{\text{VAF}} \times (1 - \text{VAF}_{\text{norm}}) +
#'                  w_{\text{smooth}} \times \text{smooth}_{\text{norm}} +
#'                  w_{\text{VSS}} \times (1 - \text{VSS}_{\text{norm}})
#' }
#'
#' The r with the minimum composite score is selected among models satisfying
#' VAF >= vaf_threshold.
#'
#' @examples
#' \dontrun{
#' result <- select_r_multicriteria_mfa(
#'   list_of_data, K = 3, rvec = 1:5,
#'   vaf_threshold = 0.9,
#'   lambda_smooth = 1e-3,
#'   weights = c(BIC = 1, VAF = 1, smoothness = 1, VSS = 1)
#' )
#' print(result$summary)
#' best_model <- result$best_model
#' }
#'
#' @export
select_r_multicriteria_mfa <- function(
    list_of_data,
    K,
    rvec = 1:6,
    vaf_threshold = 0.9,
    lambda_smooth = 1e-3,
    vss_threshold = 0.5,
    weights = c(BIC = 1, VAF = 1, smoothness = 1, VSS = 1),
    max_iter = 50,
    nIterFA = 5,
    tol = 1e-3,
    n_init = 1,
    use_kmeans_init = TRUE,
    subject_rdim_for_kmeans = 2,
    mc_cores = 1,
    n_threads = 1
) {
  N <- length(list_of_data)
  N_total_rows <- sum(sapply(list_of_data, nrow))
  M <- ncol(list_of_data[[1]])
  
  results_list <- list()
  
  for (r in rvec) {
    # Fit MFA model
    fit_mfa <- mfa_em_fit(
      list_of_data = list_of_data,
      K = K,
      r = r,
      max_iter = max_iter,
      nIterFA = nIterFA,
      tol = tol,
      n_init = n_init,
      use_kmeans_init = use_kmeans_init,
      subject_rdim_for_kmeans = subject_rdim_for_kmeans,
      mc_cores = mc_cores,
      n_threads = n_threads
    )
    
    # Compute BIC
    loglik_val <- fit_mfa$logLik
    bic_val <- compute_BIC_mfa(loglik_val, K, r, M, N_total_rows)
    
    # Compute VAF
    vaf_val <- compute_global_vaf_mfa(list_of_data, fit_mfa)
    
    # Compute factor scores
    fit_mfa <- compute_factor_scores_mfa(fit_mfa, list_of_data)
    
    # Compute smoothness penalty
    smoothness_val <- compute_smoothness_penalty(fit_mfa$factor_scores)
    
    # Compute VSS criterion
    vss_val <- compute_vss_criterion(fit_mfa$Lambda, threshold = vss_threshold)
    
    results_list[[length(results_list) + 1]] <- list(
      r = r,
      model = fit_mfa,
      logLik = loglik_val,
      BIC = bic_val,
      VAF = vaf_val,
      smoothness = smoothness_val,
      VSS = vss_val
    )
  }
  
  # Build summary data frame
  df_summary <- do.call(rbind, lapply(results_list, function(res) {
    data.frame(
      r = res$r,
      logLik = res$logLik,
      BIC = res$BIC,
      VAF = res$VAF,
      smoothness = res$smoothness,
      VSS = res$VSS
    )
  }))
  
  # Normalize each criterion to [0, 1]
  normalize <- function(x) {
    if (length(unique(x)) == 1) return(rep(0.5, length(x)))
    (x - min(x)) / (max(x) - min(x))
  }
  
  bic_norm <- normalize(df_summary$BIC)
  vaf_norm <- normalize(df_summary$VAF)
  smooth_norm <- normalize(df_summary$smoothness)
  vss_norm <- normalize(df_summary$VSS)
  
  # Compute composite score (lower is better)
  # BIC: lower is better -> use as is
  # VAF: higher is better -> use (1 - normalized)
  # Smoothness: lower is better -> use as is
  # VSS: higher is better -> use (1 - normalized)
  composite_score <- 
    weights["BIC"] * bic_norm +
    weights["VAF"] * (1 - vaf_norm) +
    weights["smoothness"] * smooth_norm +
    weights["VSS"] * (1 - vss_norm)
  
  df_summary$composite_score <- composite_score
  
  # Filter by VAF threshold
  df_valid <- df_summary[df_summary$VAF >= vaf_threshold, ]
  
  if (nrow(df_valid) == 0) {
    warning("No models satisfy VAF threshold. Returning model with highest VAF.")
    best_idx <- which.max(df_summary$VAF)
  } else {
    # Select r with minimum composite score among valid models
    best_idx_in_valid <- which.min(df_valid$composite_score)
    best_r <- df_valid$r[best_idx_in_valid]
    best_idx <- which(df_summary$r == best_r)[1]
  }
  
  best_model_info <- results_list[[best_idx]]
  
  cat("=== Best r by multi-criteria selection (MFA) ===\n")
  cat(sprintf("Selected r = %d (K = %d)\n", best_model_info$r, K))
  cat(sprintf("BIC = %.2f, VAF = %.4f, Smoothness = %.2f, VSS = %.4f\n",
              best_model_info$BIC, best_model_info$VAF,
              best_model_info$smoothness, best_model_info$VSS))
  cat(sprintf("Composite Score = %.4f\n", df_summary$composite_score[best_idx]))
  
  list(
    best_r = best_model_info$r,
    summary = df_summary,
    best_model = best_model_info$model
  )
}


#' Select Optimal r Using Multi-Criteria Approach for MPCA
#'
#' Selects the optimal number of factors (r) for a given number of clusters (K)
#' by combining multiple criteria: BIC, VAF, smoothness penalty, and VSS.
#'
#' @param list_of_data A list of N matrices, each (T_i x M).
#' @param K Number of clusters (fixed).
#' @param rvec Vector of candidate r values to evaluate. Default is 1:6.
#' @param vaf_threshold Minimum VAF threshold. Models with VAF below this are excluded.
#'   Default is 0.9.
#' @param lambda_smooth Smoothness penalty weight. Default is 1e-3.
#' @param vss_threshold Threshold for VSS criterion. Default is 0.5.
#' @param weights Named vector of weights for each criterion. Default is
#'   c(BIC = 1, VAF = 1, smoothness = 1, VSS = 1).
#' @param max_iter Maximum EM iterations.
#' @param nIterPCA Sub-iterations for PCA update.
#' @param tol Convergence tolerance.
#' @param method "EM" or "closed_form".
#' @param n_init Number of random initializations.
#' @param use_kmeans_init Whether to use k-means initialization.
#' @param subject_rdim_for_kmeans Dimension for k-means features.
#' @param mc_cores Number of cores for parallelization.
#' @param n_threads Number of OpenMP threads.
#'
#' @return A list with:
#'   \item{best_r}{The selected r value.}
#'   \item{summary}{A data.frame with columns r, BIC, VAF, smoothness, VSS, composite_score.}
#'   \item{best_model}{The fitted model for the selected r.}
#'
#' @details
#' This function combines multiple criteria to select r:
#' \itemize{
#'   \item BIC: Lower is better (penalizes complexity)
#'   \item VAF: Higher is better (rewards fit quality)
#'   \item Smoothness: Lower is better (penalizes wiggly factors)
#'   \item VSS: Higher is better (rewards simple structure)
#' }
#'
#' Each criterion is normalized to [0, 1] range, then combined using weights:
#' \deqn{
#'   \text{score} = w_{\text{BIC}} \times \text{BIC}_{\text{norm}} +
#'                  w_{\text{VAF}} \times (1 - \text{VAF}_{\text{norm}}) +
#'                  w_{\text{smooth}} \times \text{smooth}_{\text{norm}} +
#'                  w_{\text{VSS}} \times (1 - \text{VSS}_{\text{norm}})
#' }
#'
#' The r with the minimum composite score is selected among models satisfying
#' VAF >= vaf_threshold.
#'
#' @examples
#' \dontrun{
#' result <- select_r_multicriteria_mpca(
#'   list_of_data, K = 3, rvec = 1:5,
#'   vaf_threshold = 0.9,
#'   lambda_smooth = 1e-3,
#'   weights = c(BIC = 1, VAF = 1, smoothness = 1, VSS = 1)
#' )
#' print(result$summary)
#' best_model <- result$best_model
#' }
#'
#' @export
select_r_multicriteria_mpca <- function(
    list_of_data,
    K,
    rvec = 1:6,
    vaf_threshold = 0.9,
    lambda_smooth = 1e-3,
    vss_threshold = 0.5,
    weights = c(BIC = 1, VAF = 1, smoothness = 1, VSS = 1),
    max_iter = 50,
    nIterPCA = 5,
    tol = 1e-3,
    method = "EM",
    n_init = 1,
    use_kmeans_init = TRUE,
    subject_rdim_for_kmeans = 2,
    mc_cores = 1,
    n_threads = 1
) {
  N <- length(list_of_data)
  N_total_rows <- sum(sapply(list_of_data, nrow))
  M <- ncol(list_of_data[[1]])
  
  results_list <- list()
  
  for (r in rvec) {
    # Fit MPCA model
    fit_mpca <- mixture_pca_em_fit(
      list_of_data = list_of_data,
      K = K,
      r = r,
      max_iter = max_iter,
      nIterPCA = nIterPCA,
      tol = tol,
      method = method,
      n_init = n_init,
      use_kmeans_init = use_kmeans_init,
      subject_rdim_for_kmeans = subject_rdim_for_kmeans,
      mc_cores = mc_cores,
      n_threads = n_threads
    )
    
    # Ensure W field exists (reconstruct from P and D if needed)
    fit_mpca <- ensure_mpca_has_W(fit_mpca)
    
    # Compute BIC
    loglik_val <- fit_mpca$logLik
    bic_val <- compute_BIC_mpca(loglik_val, K, r, M, N_total_rows)
    
    # Compute VAF
    vaf_val <- compute_global_vaf_mpca(list_of_data, fit_mpca)
    
    # Compute factor scores
    fit_mpca <- compute_factor_scores_mpca(fit_mpca, list_of_data)
    
    # Compute smoothness penalty
    smoothness_val <- compute_smoothness_penalty(fit_mpca$factor_scores)
    
    # Compute VSS criterion
    vss_val <- compute_vss_criterion(fit_mpca$W, threshold = vss_threshold)
    
    results_list[[length(results_list) + 1]] <- list(
      r = r,
      model = fit_mpca,
      logLik = loglik_val,
      BIC = bic_val,
      VAF = vaf_val,
      smoothness = smoothness_val,
      VSS = vss_val
    )
  }
  
  # Build summary data frame
  df_summary <- do.call(rbind, lapply(results_list, function(res) {
    data.frame(
      r = res$r,
      logLik = res$logLik,
      BIC = res$BIC,
      VAF = res$VAF,
      smoothness = res$smoothness,
      VSS = res$VSS
    )
  }))
  
  # Normalize each criterion to [0, 1]
  normalize <- function(x) {
    if (length(unique(x)) == 1) return(rep(0.5, length(x)))
    (x - min(x)) / (max(x) - min(x))
  }
  
  bic_norm <- normalize(df_summary$BIC)
  vaf_norm <- normalize(df_summary$VAF)
  smooth_norm <- normalize(df_summary$smoothness)
  vss_norm <- normalize(df_summary$VSS)
  
  # Compute composite score (lower is better)
  # BIC: lower is better -> use as is
  # VAF: higher is better -> use (1 - normalized)
  # Smoothness: lower is better -> use as is
  # VSS: higher is better -> use (1 - normalized)
  composite_score <- 
    weights["BIC"] * bic_norm +
    weights["VAF"] * (1 - vaf_norm) +
    weights["smoothness"] * smooth_norm +
    weights["VSS"] * (1 - vss_norm)
  
  df_summary$composite_score <- composite_score
  
  # Filter by VAF threshold
  df_valid <- df_summary[df_summary$VAF >= vaf_threshold, ]
  
  if (nrow(df_valid) == 0) {
    warning("No models satisfy VAF threshold. Returning model with highest VAF.")
    best_idx <- which.max(df_summary$VAF)
  } else {
    # Select r with minimum composite score among valid models
    best_idx_in_valid <- which.min(df_valid$composite_score)
    best_r <- df_valid$r[best_idx_in_valid]
    best_idx <- which(df_summary$r == best_r)[1]
  }
  
  best_model_info <- results_list[[best_idx]]
  
  cat("=== Best r by multi-criteria selection (MPCA) ===\n")
  cat(sprintf("Selected r = %d (K = %d)\n", best_model_info$r, K))
  cat(sprintf("BIC = %.2f, VAF = %.4f, Smoothness = %.2f, VSS = %.4f\n",
              best_model_info$BIC, best_model_info$VAF,
              best_model_info$smoothness, best_model_info$VSS))
  cat(sprintf("Composite Score = %.4f\n", df_summary$composite_score[best_idx]))
  
  list(
    best_r = best_model_info$r,
    summary = df_summary,
    best_model = best_model_info$model
  )
}
