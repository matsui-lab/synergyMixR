#' Search for the Best (K, r) for MFA via Grid Search (Using mfa_em_fit)
#'
#' This function tries all combinations of K in Kvec and r in rvec,
#' fits a Mixture Factor Analysis (MFA) model using mfa_em_fit
#' for each combination, computes BIC, and returns the best model by BIC.
#'
#' @param list_of_data A list of data matrices, each (T_i x M).
#' @param Kvec An integer vector of candidate K values.
#' @param rvec An integer vector of candidate r values.
#' @param max_iter Maximum EM iterations.
#' @param nIterFA Sub-iterations within the factor-analyzer update in C++.
#' @param tol Convergence tolerance.
#' @param n_init Number of random initializations to try inside mfa_em_fit.
#' @param use_kmeans_init Logical; whether to also try a k-means-based init.
#' @param subject_rdim_for_kmeans Dimension for PCA-based features if use_kmeans_init=TRUE.
#' @param mc_cores_grid Number of CPU cores for grid-level parallelization (default: 1).
#' @param mc_cores Number of CPU cores for initialization-level parallelization within each (K, r) (default: 1).
#' @param n_threads Number of OpenMP threads for EM-level parallelization (default: 1 to avoid nested parallelism).
#'
#' @return A list with:
#'   \item{summary}{A data.frame with columns K, r, logLik, BIC.}
#'   \item{best_model_info}{A list with (K, r, model, logLik, BIC) for the best BIC.}
#'   \item{all_models}{A list of all fits, each (K, r, model, logLik, BIC).}
#'
#' @details
#' - Calls mfa_em_fit(..., n_init, use_kmeans_init, ...) for each (K, r).
#' - After fitting, we compute log-likelihood by summing over all data
#'   (you can use the built-in logLik from the fit, or re-compute).
#' - Then compute BIC = -2*logLik + p*log(N_total_rows).
#'   Here p = K*(M*r + 2*M) + (K - 1) by default (diagonal noise).
#'
#' @export
select_optimal_K_r_mfa <- function(
    list_of_data,
    Kvec = 1:5,
    rvec = 1:5,
    max_iter = 50,
    nIterFA  = 5,
    tol      = 1e-3,
    n_init   = 1,
    use_kmeans_init = TRUE,
    subject_rdim_for_kmeans = 2,
    mc_cores_grid = 1,
    mc_cores = 1,
    n_threads = 1
){
  # Input validation
  stopifnot(
    "list_of_data must be a list" = is.list(list_of_data),
    "No data: list_of_data cannot be empty" = length(list_of_data) > 0,
    "All elements of list_of_data must be matrices" = all(sapply(list_of_data, is.matrix)),
    "Kvec must be a numeric vector" = is.numeric(Kvec),
    "All values in Kvec must be positive integers" = all(Kvec > 0 & Kvec == floor(Kvec)),
    "rvec must be a numeric vector" = is.numeric(rvec),
    "All values in rvec must be positive integers" = all(rvec > 0 & rvec == floor(rvec)),
    "max_iter must be positive" = is.numeric(max_iter) && max_iter > 0,
    "nIterFA must be positive" = is.numeric(nIterFA) && nIterFA > 0,
    "tol must be positive" = is.numeric(tol) && tol > 0,
    "n_init must be a positive integer" = is.numeric(n_init) && n_init > 0 && n_init == floor(n_init),
    "use_kmeans_init must be logical" = is.logical(use_kmeans_init),
    "mc_cores_grid must be a positive integer" = is.numeric(mc_cores_grid) && mc_cores_grid > 0 && mc_cores_grid == floor(mc_cores_grid),
    "mc_cores must be a positive integer" = is.numeric(mc_cores) && mc_cores > 0 && mc_cores == floor(mc_cores),
    "n_threads must be a non-negative integer" = is.numeric(n_threads) && n_threads >= 0 && n_threads == floor(n_threads)
  )
  
  # Windows parallel fallback for grid-level parallelization
  if (.Platform$OS.type == "windows" && mc_cores_grid > 1) {
    warning("Windows detected: forcing sequential grid execution (mc_cores_grid = 1)")
    mc_cores_grid <- 1
  }
  
  # Parallel safety logic: prevent nested over-parallelization
  # If grid-level parallelization is enabled, force n_threads to 1
  if (mc_cores_grid > 1 && n_threads > 1) {
    message("Warning: Nested parallelism detected: forcing n_threads = 1 for stability.")
    n_threads <- 1
  }
  
  # Core allocation logic: avoid nested parallelization
  # If grid-level parallelization is enabled, disable internal parallelization
  mc_cores_internal <- if (mc_cores_grid > 1) 1 else mc_cores
  
  N <- length(list_of_data)

  # Total time samples (for BIC)
  # Here we use "N_total_rows = sum of T_i"
  N_total_rows <- sum(sapply(list_of_data, nrow))

  M <- ncol(list_of_data[[1]])  # Number of channels

  # Create parameter grid for all (K, r) combinations
  param_grid <- expand.grid(K = Kvec, r = rvec)
  total_steps <- nrow(param_grid)

  # Check if progressr is available for progress tracking
  has_progressr <- requireNamespace("progressr", quietly = TRUE)
  
  # Execute grid search with optional parallelization and progress tracking
  if (has_progressr) {
    # Use progressr for progress tracking (works with both sequential and parallel)
    # Note: Do not call handlers(global = TRUE) here as it conflicts with testthat
    # and other contexts where handlers are already on the stack.
    # Progress will work if the user has configured handlers externally.
    results_list <- progressr::with_progress({
      p <- progressr::progressor(steps = total_steps)
      
      # Define function to fit one (K, r) combination with progress
      fit_one_combination <- function(idx) {
        K <- param_grid$K[idx]
        r <- param_grid$r[idx]
        
        # Report progress
        p(message = sprintf("Evaluating K=%d, r=%d", K, r))
        
        # --- 1) Fit via mfa_em_fit (multiple inits, k-means init possible) ---
        fit_mfa <- mfa_em_fit(
          list_of_data = list_of_data,
          K = K,
          r = r,
          max_iter = max_iter,
          nIterFA  = nIterFA,
          tol      = tol,
          n_init   = n_init,
          use_kmeans_init = use_kmeans_init,
          subject_rdim_for_kmeans = subject_rdim_for_kmeans,
          mc_cores = mc_cores_internal,
          n_threads = n_threads
        )

        # --- 2) Use fit_mfa$logLik directly ---
        loglik_val <- fit_mfa$logLik

        # --- 3) compute BIC using corrected formula with factor covariances ---
        bic_val <- compute_BIC_mfa(loglik_val, K, r, M, N_total_rows)

        list(
          K=K, r=r,
          model=fit_mfa,
          logLik=loglik_val,
          BIC=bic_val
        )
      }
      
      # Execute with appropriate parallelization
      if (mc_cores_grid == 1) {
        # Sequential execution
        lapply(seq_len(total_steps), fit_one_combination)
      } else {
        # Parallel execution
        if(!requireNamespace("parallel", quietly=TRUE)){
          stop("Package 'parallel' is required for parallel grid search. Please install it.")
        }
        parallel::mclapply(
          seq_len(total_steps),
          fit_one_combination,
          mc.cores = mc_cores_grid
        )
      }
    })
  } else {
    # Fallback: use txtProgressBar for sequential, simple message for parallel
    
    # Define function to fit one (K, r) combination (no progress tracking)
    fit_one_combination <- function(idx) {
      K <- param_grid$K[idx]
      r <- param_grid$r[idx]
      
      # --- 1) Fit via mfa_em_fit (multiple inits, k-means init possible) ---
      fit_mfa <- mfa_em_fit(
        list_of_data = list_of_data,
        K = K,
        r = r,
        max_iter = max_iter,
        nIterFA  = nIterFA,
        tol      = tol,
        n_init   = n_init,
        use_kmeans_init = use_kmeans_init,
        subject_rdim_for_kmeans = subject_rdim_for_kmeans,
        mc_cores = mc_cores_internal,
        n_threads = n_threads
      )

      # --- 2) Use fit_mfa$logLik directly ---
      loglik_val <- fit_mfa$logLik

      # --- 3) compute BIC using corrected formula with factor covariances ---
      bic_val <- compute_BIC_mfa(loglik_val, K, r, M, N_total_rows)

      list(
        K=K, r=r,
        model=fit_mfa,
        logLik=loglik_val,
        BIC=bic_val
      )
    }
    
    if (mc_cores_grid == 1) {
      # Sequential execution with txtProgressBar
      pb <- txtProgressBar(min=0, max=total_steps, style=3)
      results_list <- lapply(seq_len(total_steps), function(idx) {
        setTxtProgressBar(pb, idx)
        fit_one_combination(idx)
      })
      close(pb)
    } else {
      # Parallel execution with simple message
      if(!requireNamespace("parallel", quietly=TRUE)){
        stop("Package 'parallel' is required for parallel grid search. Please install it.")
      }
      message(sprintf("Running grid search in parallel with %d cores...", mc_cores_grid))
      results_list <- parallel::mclapply(
        seq_len(total_steps),
        fit_one_combination,
        mc.cores = mc_cores_grid
      )
    }
  }

  # Build summary data frame
  df_summary <- do.call(rbind, lapply(results_list, function(res) {
    data.frame(K=res$K, r=res$r, logLik=res$logLik, BIC=res$BIC)
  }))

  # Sort by BIC
  df_summary <- df_summary[order(df_summary$BIC), ]
  best_row <- df_summary[1, ]
  message("=== Best model by BIC (MFA) ===")
  message(paste(utils::capture.output(print(best_row)), collapse = "\n"))

  best_K <- best_row$K
  best_r <- best_row$r

  # Find the best model
  best_model_index <- which(
    sapply(results_list, function(x) x$K) == best_K &
      sapply(results_list, function(x) x$r) == best_r
  )
  best_model_info <- results_list[[ best_model_index[1] ]]  # Use the first match

  list(
    summary         = df_summary,
    best_model_info = best_model_info,
    all_models      = results_list
  )
}
#' Search for the Best (K, r) for Mixture PCA via Grid Search
#'
#' This function tries all combinations of K in Kvec and r in rvec,
#' fits a Mixture PCA (PPCA) model using mixture_pca_em_fit
#' for each combination, computes BIC, and returns the best model by BIC.
#'
#' @param list_of_data A list of data matrices, each (T_i x M).
#' @param Kvec An integer vector of candidate K values.
#' @param rvec An integer vector of candidate r values.
#' @param max_iter Maximum EM iterations.
#' @param nIterPCA Sub-iterations within the PCA/PPCA update in C++.
#' @param tol Convergence tolerance.
#' @param method "EM" or "closed_form" (passed to mixture_pca_em_fit).
#' @param n_init Number of random initializations to try inside mixture_pca_em_fit.
#' @param use_kmeans_init Logical; whether to also try a k-means-based init.
#' @param subject_rdim_for_kmeans Dimension for PCA-based features if use_kmeans_init=TRUE.
#' @param mc_cores_grid Number of CPU cores for grid-level parallelization (default: 1).
#' @param mc_cores Number of CPU cores for initialization-level parallelization within each (K, r) (default: 1).
#' @param n_threads Number of OpenMP threads for EM-level parallelization (default: 1 to avoid nested parallelism).
#'
#' @return A list with:
#'   \item{summary}{A data.frame with columns K, r, logLik, BIC.}
#'   \item{best_model_info}{A list with (K, r, model, logLik, BIC) for the best BIC.}
#'   \item{all_models}{A list of all fits, each (K, r, model, logLik, BIC).}
#'
#' @details
#' - Calls mixture_pca_em_fit(..., n_init, use_kmeans_init, ...) for each (K, r).
#' - Uses fit$logLik as the final log-likelihood.
#' - BIC formula can be adjusted if needed. By default, we might do
#'   p = K*(M*r + r + M) + (K-1), for example, if we consider each cluster having:
#'   - W is M x r (or equivalently P x D)
#'   - mu is length M
#'   - sigma2 is 1 scalar
#'   - mixing proportions (K-1)
#'
#' @export
select_optimal_K_r_mpca <- function(
    list_of_data,
    Kvec = 1:5,
    rvec = 1:5,
    max_iter = 50,
    nIterPCA = 5,
    tol      = 1e-3,
    method   = "EM",
    n_init   = 1,
    use_kmeans_init = TRUE,
    subject_rdim_for_kmeans = 2,
    mc_cores_grid = 1,
    mc_cores = 1,
    n_threads = 1
){
  # Input validation
  stopifnot(
    "list_of_data must be a list" = is.list(list_of_data),
    "No data: list_of_data cannot be empty" = length(list_of_data) > 0,
    "All elements of list_of_data must be matrices" = all(sapply(list_of_data, is.matrix)),
    "Kvec must be a numeric vector" = is.numeric(Kvec),
    "All values in Kvec must be positive integers" = all(Kvec > 0 & Kvec == floor(Kvec)),
    "rvec must be a numeric vector" = is.numeric(rvec),
    "All values in rvec must be positive integers" = all(rvec > 0 & rvec == floor(rvec)),
    "max_iter must be positive" = is.numeric(max_iter) && max_iter > 0,
    "nIterPCA must be positive" = is.numeric(nIterPCA) && nIterPCA > 0,
    "tol must be positive" = is.numeric(tol) && tol > 0,
    "method must be 'EM' or 'closed_form'" = method %in% c("EM", "closed_form"),
    "n_init must be a positive integer" = is.numeric(n_init) && n_init > 0 && n_init == floor(n_init),
    "use_kmeans_init must be logical" = is.logical(use_kmeans_init),
    "mc_cores_grid must be a positive integer" = is.numeric(mc_cores_grid) && mc_cores_grid > 0 && mc_cores_grid == floor(mc_cores_grid),
    "mc_cores must be a positive integer" = is.numeric(mc_cores) && mc_cores > 0 && mc_cores == floor(mc_cores),
    "n_threads must be a non-negative integer" = is.numeric(n_threads) && n_threads >= 0 && n_threads == floor(n_threads)
  )
  
  # Windows parallel fallback for grid-level parallelization
  if (.Platform$OS.type == "windows" && mc_cores_grid > 1) {
    warning("Windows detected: forcing sequential grid execution (mc_cores_grid = 1)")
    mc_cores_grid <- 1
  }
  
  # Parallel safety logic: prevent nested over-parallelization
  # If grid-level parallelization is enabled, force n_threads to 1
  if (mc_cores_grid > 1 && n_threads > 1) {
    message("Warning: Nested parallelism detected: forcing n_threads = 1 for stability.")
    n_threads <- 1
  }
  
  # Core allocation logic: avoid nested parallelization
  # If grid-level parallelization is enabled, disable internal parallelization
  mc_cores_internal <- if (mc_cores_grid > 1) 1 else mc_cores
  
  N <- length(list_of_data)

  # 総タイムサンプル数
  N_total_rows <- sum(sapply(list_of_data, nrow))
  M <- ncol(list_of_data[[1]])

  # Create parameter grid for all (K, r) combinations
  param_grid <- expand.grid(K = Kvec, r = rvec)
  total_steps <- nrow(param_grid)

  # Check if progressr is available for progress tracking
  has_progressr <- requireNamespace("progressr", quietly = TRUE)
  
  # Execute grid search with optional parallelization and progress tracking
  if (has_progressr) {
    # Use progressr for progress tracking (works with both sequential and parallel)
    # Note: Do not call handlers(global = TRUE) here as it conflicts with testthat
    # and other contexts where handlers are already on the stack.
    # Progress will work if the user has configured handlers externally.
    results_list <- progressr::with_progress({
      p <- progressr::progressor(steps = total_steps)
      
      # Define function to fit one (K, r) combination with progress
      fit_one_combination <- function(idx) {
        K <- param_grid$K[idx]
        r <- param_grid$r[idx]
        
        # Report progress
        p(message = sprintf("Evaluating K=%d, r=%d", K, r))
        
        # --- 1) Fit via mixture_pca_em_fit ---
        fit_pca <- mixture_pca_em_fit(
          list_of_data = list_of_data,
          K = K,
          r = r,
          max_iter = max_iter,
          nIterPCA = nIterPCA,
          tol      = tol,
          method   = method,
          n_init   = n_init,
          use_kmeans_init = use_kmeans_init,
          subject_rdim_for_kmeans = subject_rdim_for_kmeans,
          mc_cores = mc_cores_internal,
          n_threads = n_threads
        )

        # --- 2) logLik ---
        loglik_val <- fit_pca$logLik

        # --- 3) compute BIC using corrected formula with factor covariances ---
        bic_val <- compute_BIC_mpca(loglik_val, K, r, M, N_total_rows)

        list(
          K=K, r=r,
          model=fit_pca,
          logLik=loglik_val,
          BIC=bic_val
        )
      }
      
      # Execute with appropriate parallelization
      if (mc_cores_grid == 1) {
        # Sequential execution
        lapply(seq_len(total_steps), fit_one_combination)
      } else {
        # Parallel execution
        if(!requireNamespace("parallel", quietly=TRUE)){
          stop("Package 'parallel' is required for parallel grid search. Please install it.")
        }
        parallel::mclapply(
          seq_len(total_steps),
          fit_one_combination,
          mc.cores = mc_cores_grid
        )
      }
    })
  } else {
    # Fallback: use txtProgressBar for sequential, simple message for parallel
    
    # Define function to fit one (K, r) combination (no progress tracking)
    fit_one_combination <- function(idx) {
      K <- param_grid$K[idx]
      r <- param_grid$r[idx]
      
      # --- 1) Fit via mixture_pca_em_fit ---
      fit_pca <- mixture_pca_em_fit(
        list_of_data = list_of_data,
        K = K,
        r = r,
        max_iter = max_iter,
        nIterPCA = nIterPCA,
        tol      = tol,
        method   = method,
        n_init   = n_init,
        use_kmeans_init = use_kmeans_init,
        subject_rdim_for_kmeans = subject_rdim_for_kmeans,
        mc_cores = mc_cores_internal,
        n_threads = n_threads
      )

      # --- 2) logLik ---
      loglik_val <- fit_pca$logLik

      # --- 3) compute BIC using corrected formula with factor covariances ---
      bic_val <- compute_BIC_mpca(loglik_val, K, r, M, N_total_rows)

      list(
        K=K, r=r,
        model=fit_pca,
        logLik=loglik_val,
        BIC=bic_val
      )
    }
    
    if (mc_cores_grid == 1) {
      # Sequential execution with txtProgressBar
      pb <- txtProgressBar(min=0, max=total_steps, style=3)
      results_list <- lapply(seq_len(total_steps), function(idx) {
        setTxtProgressBar(pb, idx)
        fit_one_combination(idx)
      })
      close(pb)
    } else {
      # Parallel execution with simple message
      if(!requireNamespace("parallel", quietly=TRUE)){
        stop("Package 'parallel' is required for parallel grid search. Please install it.")
      }
      message(sprintf("Running grid search in parallel with %d cores...", mc_cores_grid))
      results_list <- parallel::mclapply(
        seq_len(total_steps),
        fit_one_combination,
        mc.cores = mc_cores_grid
      )
    }
  }

  # Build summary data frame
  df_summary <- do.call(rbind, lapply(results_list, function(res) {
    data.frame(K=res$K, r=res$r, logLik=res$logLik, BIC=res$BIC)
  }))

  # Sort by BIC
  df_summary <- df_summary[order(df_summary$BIC), ]
  best_row <- df_summary[1, ]
  message("=== Best model by BIC (Mixture PCA) ===")
  message(paste(utils::capture.output(print(best_row)), collapse = "\n"))

  best_K <- best_row$K
  best_r <- best_row$r

  # 最良モデル
  best_model_index <- which(
    sapply(results_list, function(x) x$K) == best_K &
      sapply(results_list, function(x) x$r) == best_r
  )
  best_model_info <- results_list[[ best_model_index[1] ]]

  list(
    summary         = df_summary,
    best_model_info = best_model_info,
    all_models      = results_list
  )
}

#' Two-Stage Model Selection for MFA: Select K First, Then r
#'
#' This function implements a more stable two-stage approach to model selection:
#' Stage 1: Fix r (e.g., r=2) and search for optimal K using BIC or ICL
#' Stage 2: Fix K at the optimal value and search for optimal r using BIC
#'
#' This approach is more stable than simultaneous (K, r) grid search because
#' K and r have strong interactions in mixture models.
#'
#' @param list_of_data A list of data matrices, each (T_i x M).
#' @param Kvec An integer vector of candidate K values for Stage 1.
#' @param rvec An integer vector of candidate r values for Stage 2.
#' @param r_fixed Fixed r value to use in Stage 1 (default: 2).
#' @param criterion Model selection criterion for K selection in Stage 1.
#'   Either "BIC" (default) or "ICL". ICL adds an entropy penalty that
#'   discourages uncertain cluster assignments, making it more conservative
#'   for K selection. Stage 2 (r selection) always uses BIC.
#' @param max_iter Maximum EM iterations.
#' @param icl_max_iter Maximum EM iterations for ICL entropy calculation (default: 10).
#'   Only used when criterion="ICL". This limits iterations to maintain soft
#'   cluster assignments (non-zero entropy), making ICL effective. If NULL,
#'   defaults to 10. Uses random initialization (not K-means) internally.
#' @param icl_temperature Temperature parameter for ICL entropy scaling (default: 2.0).
#'   Higher values favor more clusters by scaling down the entropy penalty.
#' @param nIterFA Sub-iterations within the factor-analyzer update in C++.
#' @param tol Convergence tolerance.
#' @param n_init Number of random initializations to try.
#' @param use_kmeans_init Logical; whether to also try a k-means-based init.
#' @param subject_rdim_for_kmeans Dimension for PCA-based features if use_kmeans_init=TRUE.
#' @param mc_cores Number of CPU cores for initialization-level parallelization (default: 1).
#' @param n_threads Number of OpenMP threads for EM-level parallelization (default: 1).
#'
#' @return A list with:
#'   \item{stage1_summary}{Data frame with K search results (r fixed). Contains
#'     BIC and optionally ICL if criterion="ICL".}
#'   \item{stage2_summary}{Data frame with r search results (K fixed).}
#'   \item{best_K}{Optimal K from Stage 1.}
#'   \item{best_r}{Optimal r from Stage 2.}
#'   \item{best_model_info}{Final model with optimal (K, r).}
#'   \item{best_model}{The fitted model object.}
#'   \item{criterion}{The criterion used for K selection.}
#'
#' @details
#' This two-stage approach is recommended in the literature for mixture models
#' because it avoids the instability of simultaneous (K, r) optimization.
#' The default r_fixed=2 is a reasonable starting point for most applications.
#'
#' When criterion="ICL", the Integrated Complete-data Likelihood is used for
#' K selection. ICL adds an entropy penalty to BIC that discourages models
#' with uncertain cluster assignments. This makes ICL more conservative than
#' BIC and can help avoid overestimating K.
#'
#' @references
#' Biernacki, C., Celeux, G., & Govaert, G. (2000). Assessing a mixture model
#' for clustering with the integrated completed likelihood. IEEE Transactions
#' on Pattern Analysis and Machine Intelligence, 22(7), 719-725.
#'
#' @export
select_optimal_K_r_mfa_twostage <- function(
    list_of_data,
    Kvec = 2:5,
    rvec = 1:5,
    r_fixed = 2,
    criterion = c("ICL", "BIC"),
    max_iter = 50,
    icl_max_iter = NULL,
    icl_temperature = 2.0,
    nIterFA  = 5,
    tol      = 1e-3,
    n_init   = 1,
    use_kmeans_init = TRUE,
    subject_rdim_for_kmeans = 2,
    mc_cores = 1,
    n_threads = 1
){
  criterion <- match.arg(criterion)


  # ICL with limited iterations to maintain soft assignments
  if (criterion == "ICL" && is.null(icl_max_iter)) {
    icl_max_iter <- 10  # Default: limit to 10 iterations for ICL
  }

  N <- length(list_of_data)
  M <- ncol(list_of_data[[1]])
  N_total_rows <- sum(sapply(list_of_data, nrow))

  message("=== Two-Stage Model Selection for MFA ===")
  message(sprintf("Criterion for K selection: %s", criterion))
  message(sprintf("Stage 1: Searching for optimal K with r fixed at %d", r_fixed))

  # ===== Stage 1: Fix r, search for K =====
  stage1_results <- list()
  for(i in seq_along(Kvec)) {
    K <- Kvec[i]
    message(sprintf("  Fitting K=%d, r=%d...", K, r_fixed))

    fit_mfa <- mfa_em_fit(
      list_of_data = list_of_data,
      K = K,
      r = r_fixed,
      max_iter = max_iter,
      nIterFA  = nIterFA,
      tol      = tol,
      n_init   = n_init,
      use_kmeans_init = use_kmeans_init,
      subject_rdim_for_kmeans = subject_rdim_for_kmeans,
      mc_cores = mc_cores,
      n_threads = n_threads
    )

    loglik_val <- fit_mfa$logLik
    bic_val <- compute_BIC_mfa(loglik_val, K, r_fixed, M, N_total_rows)

    # Compute ICL if requested
    if (criterion == "ICL") {
      # Run limited iterations with RANDOM init to get soft assignments for ICL
      # K-means init leads to fast hard convergence, so use random init instead
      fit_icl <- mfa_em_fit(
        list_of_data = list_of_data,
        K = K,
        r = r_fixed,
        max_iter = icl_max_iter,
        nIterFA  = nIterFA,
        tol      = 1e-10,  # Very small to prevent early stopping
        n_init   = 3,      # Multiple random inits
        use_kmeans_init = FALSE,  # Random init for soft assignments
        subject_rdim_for_kmeans = subject_rdim_for_kmeans,
        mc_cores = mc_cores,
        n_threads = n_threads
      )

      # Apply temperature scaling if specified
      resp_for_icl <- fit_icl$resp
      if (icl_temperature != 1.0) {
        resp_for_icl <- apply_temperature(resp_for_icl, icl_temperature)
      }

      icl_val <- compute_ICL_mfa(fit_icl$logLik, K, r_fixed, M, N_total_rows, resp_for_icl)
      entropy_val <- compute_entropy(resp_for_icl)
    } else {
      icl_val <- NA
      entropy_val <- NA
    }

    stage1_results[[i]] <- list(
      K = K,
      r = r_fixed,
      logLik = loglik_val,
      BIC = bic_val,
      ICL = icl_val,
      entropy = entropy_val,
      model = fit_mfa
    )
  }

  # Find best K based on criterion
  stage1_df <- do.call(rbind, lapply(stage1_results, function(res) {
    data.frame(K=res$K, r=res$r, logLik=res$logLik, BIC=res$BIC,
               ICL=res$ICL, entropy=res$entropy)
  }))

  if (criterion == "ICL") {
    stage1_df <- stage1_df[order(stage1_df$ICL), ]
    best_K <- stage1_df$K[1]
    message(sprintf("\nStage 1 Result: Optimal K = %d (ICL = %.2f, BIC = %.2f)",
                best_K, stage1_df$ICL[1], stage1_df$BIC[1]))
  } else {
    stage1_df <- stage1_df[order(stage1_df$BIC), ]
    best_K <- stage1_df$K[1]
    message(sprintf("\nStage 1 Result: Optimal K = %d (BIC = %.2f)", best_K, stage1_df$BIC[1]))
  }

  # ===== Stage 2: Fix K at best_K, search for r =====
  # Stage 2 always uses BIC (r selection is less problematic)
  message(sprintf("\nStage 2: Searching for optimal r with K fixed at %d (using BIC)", best_K))

  stage2_results <- list()
  for(j in seq_along(rvec)) {
    r <- rvec[j]
    message(sprintf("  Fitting K=%d, r=%d...", best_K, r))

    fit_mfa <- mfa_em_fit(
      list_of_data = list_of_data,
      K = best_K,
      r = r,
      max_iter = max_iter,
      nIterFA  = nIterFA,
      tol      = tol,
      n_init   = n_init,
      use_kmeans_init = use_kmeans_init,
      subject_rdim_for_kmeans = subject_rdim_for_kmeans,
      mc_cores = mc_cores,
      n_threads = n_threads
    )

    loglik_val <- fit_mfa$logLik
    bic_val <- compute_BIC_mfa(loglik_val, best_K, r, M, N_total_rows)

    stage2_results[[j]] <- list(
      K = best_K,
      r = r,
      logLik = loglik_val,
      BIC = bic_val,
      model = fit_mfa
    )
  }

  # Find best r
  stage2_df <- do.call(rbind, lapply(stage2_results, function(res) {
    data.frame(K=res$K, r=res$r, logLik=res$logLik, BIC=res$BIC)
  }))
  stage2_df <- stage2_df[order(stage2_df$BIC), ]
  best_r <- stage2_df$r[1]

  message(sprintf("\nStage 2 Result: Optimal r = %d (BIC = %.2f)", best_r, stage2_df$BIC[1]))
  message(sprintf("\n=== Final Model: K=%d, r=%d ===", best_K, best_r))

  # Get the best model
  best_model_idx <- which(sapply(stage2_results, function(x) x$r) == best_r)
  best_model_info <- stage2_results[[best_model_idx[1]]]

  list(
    stage1_summary = stage1_df,
    stage2_summary = stage2_df,
    best_K = best_K,
    best_r = best_r,
    best_model_info = best_model_info,
    best_model = best_model_info$model,
    criterion = criterion
  )
}

#' Two-Stage Model Selection for MPCA: Select K First, Then r
#'
#' This function implements a more stable two-stage approach to model selection:
#' Stage 1: Fix r (e.g., r=2) and search for optimal K using BIC or ICL
#' Stage 2: Fix K at the optimal value and search for optimal r using BIC
#'
#' This approach is more stable than simultaneous (K, r) grid search because
#' K and r have strong interactions in mixture models.
#'
#' @param list_of_data A list of data matrices, each (T_i x M).
#' @param Kvec An integer vector of candidate K values for Stage 1.
#' @param rvec An integer vector of candidate r values for Stage 2.
#' @param r_fixed Fixed r value to use in Stage 1 (default: 2).
#' @param criterion Model selection criterion for K selection in Stage 1.
#'   Either "BIC" (default) or "ICL". ICL adds an entropy penalty that
#'   discourages uncertain cluster assignments, making it more conservative
#'   for K selection. Stage 2 (r selection) always uses BIC.
#' @param max_iter Maximum EM iterations.
#' @param icl_max_iter Maximum EM iterations for ICL entropy calculation (default: 10).
#'   Only used when criterion="ICL". This limits iterations to maintain soft
#'   cluster assignments (non-zero entropy), making ICL effective. If NULL,
#'   defaults to 10. Uses random initialization (not K-means) internally.
#' @param icl_temperature Temperature parameter for ICL entropy scaling (default: 2.0).
#'   Higher values favor more clusters by scaling down the entropy penalty.
#' @param nIterPCA Sub-iterations within the PCA/PPCA update in C++.
#' @param tol Convergence tolerance.
#' @param method "EM" or "closed_form" (passed to mixture_pca_em_fit).
#' @param n_init Number of random initializations to try.
#' @param use_kmeans_init Logical; whether to also try a k-means-based init.
#' @param subject_rdim_for_kmeans Dimension for PCA-based features if use_kmeans_init=TRUE.
#' @param mc_cores Number of CPU cores for initialization-level parallelization (default: 1).
#' @param n_threads Number of OpenMP threads for EM-level parallelization (default: 1).
#'
#' @return A list with:
#'   \item{stage1_summary}{Data frame with K search results (r fixed). Contains
#'     BIC and optionally ICL if criterion="ICL".}
#'   \item{stage2_summary}{Data frame with r search results (K fixed).}
#'   \item{best_K}{Optimal K from Stage 1.}
#'   \item{best_r}{Optimal r from Stage 2.}
#'   \item{best_model_info}{Final model with optimal (K, r).}
#'   \item{best_model}{The fitted model object.}
#'   \item{criterion}{The criterion used for K selection.}
#'
#' @details
#' This two-stage approach is recommended in the literature for mixture models
#' because it avoids the instability of simultaneous (K, r) optimization.
#' The default r_fixed=2 is a reasonable starting point for most applications.
#'
#' When criterion="ICL", the Integrated Complete-data Likelihood is used for
#' K selection. ICL adds an entropy penalty to BIC that discourages models
#' with uncertain cluster assignments. This makes ICL more conservative than
#' BIC and can help avoid overestimating K.
#'
#' @references
#' Biernacki, C., Celeux, G., & Govaert, G. (2000). Assessing a mixture model
#' for clustering with the integrated completed likelihood. IEEE Transactions
#' on Pattern Analysis and Machine Intelligence, 22(7), 719-725.
#'
#' @export
select_optimal_K_r_mpca_twostage <- function(
    list_of_data,
    Kvec = 2:5,
    rvec = 1:5,
    r_fixed = 2,
    criterion = c("ICL", "BIC"),
    max_iter = 50,
    icl_max_iter = NULL,
    icl_temperature = 2.0,
    nIterPCA = 5,
    tol      = 1e-3,
    method   = "EM",
    n_init   = 1,
    use_kmeans_init = TRUE,
    subject_rdim_for_kmeans = 2,
    mc_cores = 1,
    n_threads = 1
){
  criterion <- match.arg(criterion)

  # ICL with limited iterations to maintain soft assignments
  if (criterion == "ICL" && is.null(icl_max_iter)) {
    icl_max_iter <- 10  # Default: limit to 10 iterations for ICL
  }

  N <- length(list_of_data)
  M <- ncol(list_of_data[[1]])
  N_total_rows <- sum(sapply(list_of_data, nrow))

  message("=== Two-Stage Model Selection for MPCA ===")
  message(sprintf("Criterion for K selection: %s", criterion))
  message(sprintf("Stage 1: Searching for optimal K with r fixed at %d", r_fixed))

  # ===== Stage 1: Fix r, search for K =====
  stage1_results <- list()
  for(i in seq_along(Kvec)) {
    K <- Kvec[i]
    message(sprintf("  Fitting K=%d, r=%d...", K, r_fixed))

    fit_pca <- mixture_pca_em_fit(
      list_of_data = list_of_data,
      K = K,
      r = r_fixed,
      max_iter = max_iter,
      nIterPCA = nIterPCA,
      tol      = tol,
      method   = method,
      n_init   = n_init,
      use_kmeans_init = use_kmeans_init,
      subject_rdim_for_kmeans = subject_rdim_for_kmeans,
      mc_cores = mc_cores,
      n_threads = n_threads
    )

    loglik_val <- fit_pca$logLik
    bic_val <- compute_BIC_mpca(loglik_val, K, r_fixed, M, N_total_rows)

    # Compute ICL if requested
    if (criterion == "ICL") {
      # Run limited iterations with RANDOM init to get soft assignments for ICL
      # K-means init leads to fast hard convergence, so use random init instead
      fit_icl <- mixture_pca_em_fit(
        list_of_data = list_of_data,
        K = K,
        r = r_fixed,
        max_iter = icl_max_iter,
        nIterPCA = nIterPCA,
        tol      = 1e-10,  # Very small to prevent early stopping
        method   = method,
        n_init   = 3,      # Multiple random inits
        use_kmeans_init = FALSE,  # Random init for soft assignments
        subject_rdim_for_kmeans = subject_rdim_for_kmeans,
        mc_cores = mc_cores,
        n_threads = n_threads
      )

      # Apply temperature scaling if specified
      resp_for_icl <- fit_icl$resp
      if (icl_temperature != 1.0) {
        resp_for_icl <- apply_temperature(resp_for_icl, icl_temperature)
      }

      icl_val <- compute_ICL_mpca(fit_icl$logLik, K, r_fixed, M, N_total_rows, resp_for_icl)
      entropy_val <- compute_entropy(resp_for_icl)
    } else {
      icl_val <- NA
      entropy_val <- NA
    }

    stage1_results[[i]] <- list(
      K = K,
      r = r_fixed,
      logLik = loglik_val,
      BIC = bic_val,
      ICL = icl_val,
      entropy = entropy_val,
      model = fit_pca
    )
  }

  # Find best K based on criterion
  stage1_df <- do.call(rbind, lapply(stage1_results, function(res) {
    data.frame(K=res$K, r=res$r, logLik=res$logLik, BIC=res$BIC,
               ICL=res$ICL, entropy=res$entropy)
  }))

  if (criterion == "ICL") {
    stage1_df <- stage1_df[order(stage1_df$ICL), ]
    best_K <- stage1_df$K[1]
    message(sprintf("\nStage 1 Result: Optimal K = %d (ICL = %.2f, BIC = %.2f)",
                best_K, stage1_df$ICL[1], stage1_df$BIC[1]))
  } else {
    stage1_df <- stage1_df[order(stage1_df$BIC), ]
    best_K <- stage1_df$K[1]
    message(sprintf("\nStage 1 Result: Optimal K = %d (BIC = %.2f)", best_K, stage1_df$BIC[1]))
  }

  # ===== Stage 2: Fix K at best_K, search for r =====
  # Stage 2 always uses BIC (r selection is less problematic)
  message(sprintf("\nStage 2: Searching for optimal r with K fixed at %d (using BIC)", best_K))

  stage2_results <- list()
  for(j in seq_along(rvec)) {
    r <- rvec[j]
    message(sprintf("  Fitting K=%d, r=%d...", best_K, r))

    fit_pca <- mixture_pca_em_fit(
      list_of_data = list_of_data,
      K = best_K,
      r = r,
      max_iter = max_iter,
      nIterPCA = nIterPCA,
      tol      = tol,
      method   = method,
      n_init   = n_init,
      use_kmeans_init = use_kmeans_init,
      subject_rdim_for_kmeans = subject_rdim_for_kmeans,
      mc_cores = mc_cores,
      n_threads = n_threads
    )

    loglik_val <- fit_pca$logLik
    bic_val <- compute_BIC_mpca(loglik_val, best_K, r, M, N_total_rows)

    stage2_results[[j]] <- list(
      K = best_K,
      r = r,
      logLik = loglik_val,
      BIC = bic_val,
      model = fit_pca
    )
  }

  # Find best r
  stage2_df <- do.call(rbind, lapply(stage2_results, function(res) {
    data.frame(K=res$K, r=res$r, logLik=res$logLik, BIC=res$BIC)
  }))
  stage2_df <- stage2_df[order(stage2_df$BIC), ]
  best_r <- stage2_df$r[1]

  message(sprintf("\nStage 2 Result: Optimal r = %d (BIC = %.2f)", best_r, stage2_df$BIC[1]))
  message(sprintf("\n=== Final Model: K=%d, r=%d ===", best_K, best_r))

  # Get the best model
  best_model_idx <- which(sapply(stage2_results, function(x) x$r) == best_r)
  best_model_info <- stage2_results[[best_model_idx[1]]]

  list(
    stage1_summary = stage1_df,
    stage2_summary = stage2_df,
    best_K = best_K,
    best_r = best_r,
    best_model_info = best_model_info,
    best_model = best_model_info$model,
    criterion = criterion
  )
}
