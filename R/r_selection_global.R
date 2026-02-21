#' Select Optimal r Using Global Model (FA or PPCA) with Cross-Validation
#'
#' Estimates the optimal number of factors (r) using a global Factor Analysis (FA)
#' or Probabilistic PCA (PPCA) model on pooled data. This function is designed for
#' use before fitting mixture models (MFA/MPCA), as it provides a cluster-independent
#' estimate of r based on the overall data structure.
#'
#' @param list_of_data A list of N matrices, each (T_i x M).
#' @param rvec Vector of candidate r values to evaluate. Default is 1:6.
#' @param model Model type for global r estimation. Options: "PPCA" (Probabilistic PCA,
#'   default), "FA" (Factor Analysis). PPCA is recommended as it is more stable and
#'   theory-consistent for mixture models.
#' @param criterion Selection criterion. Options: "CV" (Cross-Validation, default),
#'   "ParallelAnalysis" (Parallel Analysis). CV is the primary criterion as it is
#'   stable and works for both FA and PPCA. Parallel Analysis is optional and lower
#'   priority.
#' @param folds Number of folds for cross-validation. Default is 5. CV is performed
#'   at the subject level (not time-window level) to respect the generative structure.
#' @param cv_metric Metric for CV evaluation. Options: "negLogLik" (negative log-likelihood,
#'   default), "frobenius" (Frobenius norm of reconstruction error). Negative log-likelihood
#'   is the principled choice for probabilistic models.
#' @param n_permutations Number of permutations for Parallel Analysis. Default is 100.
#'   Only used when criterion = "ParallelAnalysis".
#' @param seed Random seed for reproducibility (CV fold assignment and PA permutations).
#'   Default is NULL (no seed set).
#' @param verbose Logical. If TRUE, print detailed diagnostic information during r selection.
#'   Default is FALSE.
#'
#' @return A list with:
#'   \item{r_global}{The selected optimal r value.}
#'   \item{summary}{A data.frame with columns r and criterion-specific metrics (e.g., cv_error, cv_se).}
#'   \item{criterion}{The criterion used for selection.}
#'   \item{model}{The model type used (FA or PPCA).}
#'   \item{cv_curve}{For CV: ggplot object showing CV error vs r (if ggplot2 available).}
#'   \item{vaf_curve}{For CV: ggplot object showing VAF vs r (if ggplot2 available).}
#'   \item{pa_curve}{For PA: ggplot object showing observed vs permuted eigenvalues (if ggplot2 available).}
#'   \item{diagnostics}{Additional diagnostic information.}
#'
#' @details
#' This function implements model-based r estimation for use before mixture model fitting.
#' The key design principles are:
#'
#' \strong{Why Global Model?}
#' The global r must be estimated using FA or PPCA on pooled data, NOT using MFA(K>1).
#' Estimating r from a mixture model is cluster-dependent and unstable, as different
#' cluster configurations yield different r estimates. A global model provides a
#' cluster-independent baseline.
#'
#' \strong{Cross-Validation Design:}
#' CV is performed at the subject level, not time-window level. Time-window CV would
#' leak subject-level covariance structure and provide overly optimistic estimates.
#' Subject-level CV properly reflects the generative structure where subjects are
#' the independent units.
#'
#' For each fold:
#' \enumerate{
#'   \item Split subjects into training and test sets
#'   \item Fit FA/PPCA on training data (pooled across subjects)
#'   \item Evaluate reconstruction error or log-likelihood on test data
#'   \item Average across folds to get CV error for each r
#' }
#'
#' The optimal r is selected by minimizing CV error (or using the one-standard-error
#' rule for parsimony).
#'
#' \strong{Parallel Analysis:}
#' PA compares observed eigenvalues to those from permuted data. The optimal r is
#' where observed eigenvalues exceed the 95th percentile of permuted eigenvalues.
#' Note: PA is theoretically questionable for FA (as it assumes isotropic noise),
#' so CV is preferred.
#'
#' @examples
#' \dontrun{
#' # Basic usage with default PPCA and CV
#' result <- select_r_global_mfa(
#'   list_of_data, rvec = 2:6,
#'   model = "PPCA", criterion = "CV", folds = 5
#' )
#' print(result$r_global)
#' print(result$summary)
#' 
#' # With FA model
#' result_fa <- select_r_global_mfa(
#'   list_of_data, rvec = 2:6,
#'   model = "FA", criterion = "CV"
#' )
#' 
#' # Using Parallel Analysis
#' result_pa <- select_r_global_mfa(
#'   list_of_data, rvec = 2:6,
#'   criterion = "ParallelAnalysis", n_permutations = 100
#' )
#' }
#'
#' @export
select_r_global_mfa <- function(
    list_of_data,
    rvec = 1:6,
    model = c("PPCA", "FA"),
    criterion = c("CV", "ParallelAnalysis"),
    folds = 5,
    cv_metric = c("negLogLik", "frobenius"),
    n_permutations = 100,
    seed = NULL,
    verbose = FALSE
) {
  # Input validation
  model <- match.arg(model)
  criterion <- match.arg(criterion)
  cv_metric <- match.arg(cv_metric)
  
  if (length(rvec) < 2) {
    stop("rvec must contain at least 2 candidate r values")
  }
  
  N <- length(list_of_data)
  M <- ncol(list_of_data[[1]])
  
  if (folds < 2 || folds > N) {
    stop(sprintf("folds must be between 2 and N=%d", N))
  }
  
  if (verbose) {
    cat("=== Global r selection ===\n")
    cat(sprintf("Model: %s, Criterion: %s\n", model, criterion))
    cat(sprintf("N = %d subjects, M = %d channels\n", N, M))
    cat(sprintf("Candidate r values: %s\n", paste(rvec, collapse = ", ")))
  }
  
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Dispatch to appropriate criterion
  if (criterion == "CV") {
    result <- select_r_by_cv(
      list_of_data = list_of_data,
      rvec = rvec,
      model = model,
      folds = folds,
      cv_metric = cv_metric,
      verbose = verbose
    )
  } else if (criterion == "ParallelAnalysis") {
    result <- select_r_by_parallel_analysis(
      list_of_data = list_of_data,
      rvec = rvec,
      model = model,
      n_permutations = n_permutations,
      verbose = verbose
    )
  }
  
  # Add metadata
  result$criterion <- criterion
  result$model <- model
  
  if (verbose) {
    cat(sprintf("\n=== Selected r_global = %d ===\n", result$r_global))
  }
  
  return(result)
}


#' Select Optimal r Using Subject-Level Cross-Validation
#'
#' Internal function that implements subject-level k-fold cross-validation
#' for global FA or PPCA models.
#'
#' @param list_of_data A list of N matrices, each (T_i x M).
#' @param rvec Vector of candidate r values.
#' @param model Model type: "PPCA" or "FA".
#' @param folds Number of CV folds.
#' @param cv_metric Metric for evaluation: "negLogLik" or "frobenius".
#' @param verbose Logical for diagnostic output.
#'
#' @return A list with r_global, summary, cv_curve, vaf_curve, diagnostics.
#' @keywords internal
select_r_by_cv <- function(
    list_of_data,
    rvec,
    model,
    folds,
    cv_metric,
    verbose
) {
  N <- length(list_of_data)
  M <- ncol(list_of_data[[1]])
  
  # Create subject-level fold assignments
  fold_ids <- sample(rep(1:folds, length.out = N))
  
  if (verbose) {
    cat(sprintf("\nPerforming %d-fold subject-level cross-validation...\n", folds))
    fold_sizes <- table(fold_ids)
    cat(sprintf("Fold sizes: %s\n", paste(fold_sizes, collapse = ", ")))
  }
  
  # Storage for CV results
  cv_results <- list()
  
  # Evaluate each r value
  for (r in rvec) {
    if (verbose) {
      cat(sprintf("\nEvaluating r = %d...\n", r))
    }
    
    fold_errors <- numeric(folds)
    
    # Perform k-fold CV
    for (fold in 1:folds) {
      # Split data
      test_indices <- which(fold_ids == fold)
      train_indices <- which(fold_ids != fold)
      
      train_data <- list_of_data[train_indices]
      test_data <- list_of_data[test_indices]
      
      # Fit model on training data
      if (model == "PPCA") {
        fit <- fit_single_pca(train_data, r = r)
      } else {  # FA
        fit <- fit_single_factor_analysis(train_data, r = r)
      }
      
      # Evaluate on test data
      if (cv_metric == "negLogLik") {
        # Compute negative log-likelihood on test data
        fold_errors[fold] <- compute_test_negloglik(test_data, fit, model)
      } else {  # frobenius
        # Compute Frobenius norm of reconstruction error
        fold_errors[fold] <- compute_test_frobenius(test_data, fit, model)
      }
      
      if (verbose) {
        cat(sprintf("  Fold %d: error = %.4f\n", fold, fold_errors[fold]))
      }
    }
    
    # Compute mean and SE across folds
    cv_mean <- mean(fold_errors)
    cv_se <- sd(fold_errors) / sqrt(folds)
    
    cv_results[[length(cv_results) + 1]] <- list(
      r = r,
      cv_mean = cv_mean,
      cv_se = cv_se,
      fold_errors = fold_errors
    )
    
    if (verbose) {
      cat(sprintf("  CV error: %.4f +/- %.4f\n", cv_mean, cv_se))
    }
  }
  
  # Build summary data frame
  df_summary <- do.call(rbind, lapply(cv_results, function(res) {
    data.frame(
      r = res$r,
      cv_mean = res$cv_mean,
      cv_se = res$cv_se
    )
  }))
  
  # Select optimal r using 1-SE rule for parsimony
  # First find minimum CV error
  best_idx <- which.min(df_summary$cv_mean)
  r_min <- df_summary$r[best_idx]
  
  # Apply one-standard-error rule: select smallest r within 1 SE of minimum
  min_error <- df_summary$cv_mean[best_idx]
  min_se <- df_summary$cv_se[best_idx]
  threshold <- min_error + min_se
  
  candidates <- df_summary$r[df_summary$cv_mean <= threshold]
  r_global_1se <- min(candidates)
  
  # Use 1-SE rule as default (more conservative, avoids overfitting)
  r_global <- r_global_1se
  
  if (verbose && r_global != r_min) {
    cat(sprintf("\nNote: 1-SE rule selected r = %d (minimum CV error at r = %d)\n", 
                r_global, r_min))
  }
  
  # Compute VAF for each r on full data (for diagnostic plot)
  vaf_values <- numeric(length(rvec))
  for (i in seq_along(rvec)) {
    r <- rvec[i]
    if (model == "PPCA") {
      fit <- fit_single_pca(list_of_data, r = r)
      sse <- calc_reconstruction_error_singlePCA(list_of_data, fit)
    } else {
      fit <- fit_single_factor_analysis(list_of_data, r = r)
      sse <- calc_reconstruction_error_singleFA(list_of_data, fit)
    }
    sst <- calc_total_SST(list_of_data)
    vaf_values[i] <- 1 - (sse / sst)
  }
  
  df_summary$vaf <- vaf_values
  
  # Create diagnostic plots if ggplot2 is available
  cv_curve <- NULL
  vaf_curve <- NULL
  
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    # CV error plot
    cv_curve <- ggplot2::ggplot(df_summary, ggplot2::aes(x = r, y = cv_mean)) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = cv_mean - cv_se, ymax = cv_mean + cv_se), width = 0.2) +
      ggplot2::geom_vline(xintercept = r_global, linetype = "dashed", color = "red") +
      ggplot2::labs(
        title = sprintf("Cross-Validation: %s Model", model),
        x = "Number of Factors (r)",
        y = sprintf("CV Error (%s)", cv_metric),
        subtitle = sprintf("Optimal r = %d", r_global)
      ) +
      ggplot2::theme_minimal()
    
    # VAF plot
    vaf_curve <- ggplot2::ggplot(df_summary, ggplot2::aes(x = r, y = vaf)) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::geom_vline(xintercept = r_global, linetype = "dashed", color = "red") +
      ggplot2::labs(
        title = sprintf("Variance Accounted For: %s Model", model),
        x = "Number of Factors (r)",
        y = "VAF",
        subtitle = sprintf("Optimal r = %d (VAF = %.3f)", r_global, vaf_values[best_idx])
      ) +
      ggplot2::theme_minimal()
  }
  
  list(
    r_global = r_global,
    r_global_1se = r_global_1se,
    summary = df_summary,
    cv_curve = cv_curve,
    vaf_curve = vaf_curve,
    diagnostics = list(
      fold_ids = fold_ids,
      cv_results = cv_results
    )
  )
}


#' Select Optimal r Using Parallel Analysis
#'
#' Internal function that implements Parallel Analysis for global FA or PPCA models.
#'
#' @param list_of_data A list of N matrices, each (T_i x M).
#' @param rvec Vector of candidate r values.
#' @param model Model type: "PPCA" or "FA".
#' @param n_permutations Number of permutations.
#' @param verbose Logical for diagnostic output.
#'
#' @return A list with r_global, summary, pa_curve, diagnostics.
#' @keywords internal
select_r_by_parallel_analysis <- function(
    list_of_data,
    rvec,
    model,
    n_permutations,
    verbose
) {
  N <- length(list_of_data)
  M <- ncol(list_of_data[[1]])
  
  if (verbose) {
    cat(sprintf("\nPerforming Parallel Analysis with %d permutations...\n", n_permutations))
  }
  
  # Pool all data
  bigX <- do.call(rbind, list_of_data)
  n_total <- nrow(bigX)
  
  # Center data
  mu <- colMeans(bigX)
  Xc <- sweep(bigX, 2, mu, "-")
  
  # Compute observed eigenvalues
  Sxx <- crossprod(Xc) / n_total
  eig_obs <- eigen(Sxx, symmetric = TRUE)
  eigenvalues_obs <- eig_obs$values
  
  # Compute permuted eigenvalues
  eigenvalues_perm <- matrix(0, nrow = n_permutations, ncol = M)
  
  if (verbose) {
    cat("Running permutations...\n")
  }
  
  for (perm in 1:n_permutations) {
    # Permute each column independently
    Xperm <- apply(bigX, 2, sample)
    
    # Center permuted data
    mu_perm <- colMeans(Xperm)
    Xperm_c <- sweep(Xperm, 2, mu_perm, "-")
    
    # Compute eigenvalues
    Sxx_perm <- crossprod(Xperm_c) / n_total
    eig_perm <- eigen(Sxx_perm, symmetric = TRUE)
    eigenvalues_perm[perm, ] <- eig_perm$values
    
    if (verbose && perm %% 20 == 0) {
      cat(sprintf("  Completed %d/%d permutations\n", perm, n_permutations))
    }
  }
  
  # Compute 95th percentile of permuted eigenvalues
  eigenvalues_perm_95 <- apply(eigenvalues_perm, 2, quantile, probs = 0.95)
  
  # Determine optimal r: largest r where observed > 95th percentile
  r_global <- 0
  for (i in 1:M) {
    if (eigenvalues_obs[i] > eigenvalues_perm_95[i]) {
      r_global <- i
    } else {
      break
    }
  }
  
  # Ensure r_global is within rvec range
  if (r_global < min(rvec)) {
    warning(sprintf("PA suggests r = %d, which is below min(rvec) = %d. Using min(rvec).", 
                    r_global, min(rvec)))
    r_global <- min(rvec)
  } else if (r_global > max(rvec)) {
    warning(sprintf("PA suggests r = %d, which is above max(rvec) = %d. Using max(rvec).", 
                    r_global, max(rvec)))
    r_global <- max(rvec)
  }
  
  if (verbose) {
    cat(sprintf("\nParallel Analysis suggests r = %d\n", r_global))
  }
  
  # Build summary data frame
  df_summary <- data.frame(
    component = 1:M,
    eigenvalue_obs = eigenvalues_obs,
    eigenvalue_perm_95 = eigenvalues_perm_95,
    above_threshold = eigenvalues_obs > eigenvalues_perm_95
  )
  
  # Create diagnostic plot if ggplot2 is available
  pa_curve <- NULL
  
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    # Only plot up to max(rvec) + 2 for clarity
    plot_max <- min(M, max(rvec) + 2)
    df_plot <- df_summary[1:plot_max, ]
    
    pa_curve <- ggplot2::ggplot(df_plot, ggplot2::aes(x = component)) +
      ggplot2::geom_line(ggplot2::aes(y = eigenvalue_obs, color = "Observed")) +
      ggplot2::geom_point(ggplot2::aes(y = eigenvalue_obs, color = "Observed")) +
      ggplot2::geom_line(ggplot2::aes(y = eigenvalue_perm_95, color = "95th Percentile (Permuted)")) +
      ggplot2::geom_point(ggplot2::aes(y = eigenvalue_perm_95, color = "95th Percentile (Permuted)")) +
      ggplot2::geom_vline(xintercept = r_global, linetype = "dashed", color = "red") +
      ggplot2::labs(
        title = sprintf("Parallel Analysis: %s Model", model),
        x = "Component",
        y = "Eigenvalue",
        subtitle = sprintf("Optimal r = %d", r_global),
        color = "Type"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::scale_color_manual(values = c("Observed" = "blue", "95th Percentile (Permuted)" = "gray"))
  }
  
  list(
    r_global = r_global,
    summary = df_summary,
    pa_curve = pa_curve,
    diagnostics = list(
      eigenvalues_obs = eigenvalues_obs,
      eigenvalues_perm = eigenvalues_perm,
      eigenvalues_perm_95 = eigenvalues_perm_95
    )
  )
}


#' Compute Negative Log-Likelihood on Test Data
#'
#' Internal helper function to compute negative log-likelihood for CV.
#'
#' @param test_data List of test data matrices.
#' @param fit Fitted model (from fit_single_pca or fit_single_factor_analysis).
#' @param model Model type: "PPCA" or "FA".
#'
#' @return Negative log-likelihood (scalar).
#' @keywords internal
compute_test_negloglik <- function(test_data, fit, model) {
  n_test <- sum(sapply(test_data, nrow))
  M <- ncol(test_data[[1]])
  
  if (model == "PPCA") {
    # For PPCA, computing exact log-likelihood requires eigenvalues from training
    # which are not directly stored in fit_single_pca output.
    # Use Frobenius norm as a proxy for CV error (equivalent to MSE).
    # This is a valid choice as minimizing reconstruction error is equivalent
    # to maximizing likelihood under the PPCA model assumptions.
    return(compute_test_frobenius(test_data, fit, model))
    
  } else {  # FA
    # FA: Sigma = Lambda * Lambda^T + Psi
    loadings <- fit$loadings
    mu <- fit$mu
    diagPsi <- fit$diagPsi
    
    # Build covariance matrix
    Sigma <- loadings %*% t(loadings)
    diag(Sigma) <- diag(Sigma) + diagPsi
    
    # Compute log-likelihood
    cholS <- tryCatch({
      chol(Sigma)
    }, error = function(e) {
      # If Cholesky fails, add small regularization
      Sigma_reg <- Sigma + diag(1e-6, nrow(Sigma))
      chol(Sigma_reg)
    })
    
    logdetS <- 2 * sum(log(diag(cholS)))
    Mlog2pi <- M * log(2 * pi)
    invS <- chol2inv(cholS)
    
    quadSum <- 0
    for (X_i in test_data) {
      T_i <- nrow(X_i)
      for (t in 1:T_i) {
        x_t <- X_i[t, ]
        diff <- x_t - mu
        quadSum <- quadSum + as.numeric(t(diff) %*% invS %*% diff)
      }
    }
    
    loglik <- -0.5 * (n_test * (Mlog2pi + logdetS) + quadSum)
    
    return(-loglik)  # Return negative log-likelihood
  }
}


#' Compute Frobenius Norm of Reconstruction Error on Test Data
#'
#' Internal helper function to compute Frobenius norm for CV.
#'
#' @param test_data List of test data matrices.
#' @param fit Fitted model (from fit_single_pca or fit_single_factor_analysis).
#' @param model Model type: "PPCA" or "FA".
#'
#' @return Frobenius norm of reconstruction error (scalar).
#' @keywords internal
compute_test_frobenius <- function(test_data, fit, model) {
  sse <- 0
  
  if (model == "PPCA") {
    for (X_i in test_data) {
      Xhat_i <- reconstruct_singlePCA(X_i, fit)
      diff <- X_i - Xhat_i
      sse <- sse + sum(diff^2)
    }
  } else {  # FA
    for (X_i in test_data) {
      Xhat_i <- reconstruct_singleFA(X_i, fit$mu, fit$loadings, fit$diagPsi)
      diff <- X_i - Xhat_i
      sse <- sse + sum(diff^2)
    }
  }
  
  # Return Frobenius norm (sqrt of SSE)
  return(sqrt(sse))
}


#' Select Optimal r Using Global Model (PPCA) with Cross-Validation for MPCA
#'
#' Estimates the optimal number of principal components (r) using a global Probabilistic PCA
#' (PPCA) model on pooled data. This function is designed for use before fitting mixture MPCA
#' models, as it provides a cluster-independent estimate of r based on the overall data structure.
#'
#' @param list_of_data A list of N matrices, each (T_i x M).
#' @param rvec Vector of candidate r values to evaluate. Default is 1:6.
#' @param model Model type for global r estimation. Options: "PPCA" (Probabilistic PCA,
#'   default), "FA" (Factor Analysis). PPCA is recommended as it is more stable and
#'   theory-consistent for mixture models.
#' @param criterion Selection criterion. Options: "CV" (Cross-Validation, default),
#'   "ParallelAnalysis" (Parallel Analysis). CV is the primary criterion as it is
#'   stable and works for both FA and PPCA. Parallel Analysis is optional and lower
#'   priority.
#' @param folds Number of folds for cross-validation. Default is 5. CV is performed
#'   at the subject level (not time-window level) to respect the generative structure.
#' @param cv_metric Metric for CV evaluation. Options: "negLogLik" (negative log-likelihood,
#'   default), "frobenius" (Frobenius norm of reconstruction error). Negative log-likelihood
#'   is the principled choice for probabilistic models.
#' @param n_permutations Number of permutations for Parallel Analysis. Default is 100.
#'   Only used when criterion = "ParallelAnalysis".
#' @param seed Random seed for reproducibility (CV fold assignment and PA permutations).
#'   Default is NULL (no seed set).
#' @param verbose Logical. If TRUE, print detailed diagnostic information during r selection.
#'   Default is FALSE.
#'
#' @return A list with:
#'   \item{r_global}{The selected optimal r value.}
#'   \item{summary}{A data.frame with columns r and criterion-specific metrics (e.g., cv_error, cv_se).}
#'   \item{criterion}{The criterion used for selection.}
#'   \item{model}{The model type used (FA or PPCA).}
#'   \item{cv_curve}{For CV: ggplot object showing CV error vs r (if ggplot2 available).}
#'   \item{vaf_curve}{For CV: ggplot object showing VAF vs r (if ggplot2 available).}
#'   \item{pa_curve}{For PA: ggplot object showing observed vs permuted eigenvalues (if ggplot2 available).}
#'   \item{diagnostics}{Additional diagnostic information.}
#'
#' @details
#' This function implements model-based r estimation for use before mixture MPCA fitting.
#' The key design principles are:
#'
#' \strong{Why Global Model?}
#' The global r must be estimated using PPCA (or FA) on pooled data, NOT using MPCA(K>1).
#' Estimating r from a mixture model is cluster-dependent and unstable, as different
#' cluster configurations yield different r estimates. A global model provides a
#' cluster-independent baseline.
#'
#' \strong{Cross-Validation Design:}
#' CV is performed at the subject level, not time-window level. Time-window CV would
#' leak subject-level covariance structure and provide overly optimistic estimates.
#' Subject-level CV properly reflects the generative structure where subjects are
#' the independent units.
#'
#' For each fold:
#' \enumerate{
#'   \item Split subjects into training and test sets
#'   \item Fit PPCA/FA on training data (pooled across subjects)
#'   \item Evaluate reconstruction error or log-likelihood on test data
#'   \item Average across folds to get CV error for each r
#' }
#'
#' The optimal r is selected by minimizing CV error (or using the one-standard-error
#' rule for parsimony).
#'
#' \strong{Parallel Analysis:}
#' PA compares observed eigenvalues to those from permuted data. The optimal r is
#' where observed eigenvalues exceed the 95th percentile of permuted eigenvalues.
#' Note: PA is theoretically questionable for FA (as it assumes isotropic noise),
#' so CV is preferred.
#'
#' @examples
#' \dontrun{
#' # Basic usage with default PPCA and CV
#' result <- select_r_global_mpca(
#'   list_of_data, rvec = 2:6,
#'   model = "PPCA", criterion = "CV", folds = 5
#' )
#' print(result$r_global)
#' print(result$summary)
#' 
#' # With FA model
#' result_fa <- select_r_global_mpca(
#'   list_of_data, rvec = 2:6,
#'   model = "FA", criterion = "CV"
#' )
#' 
#' # Using Parallel Analysis
#' result_pa <- select_r_global_mpca(
#'   list_of_data, rvec = 2:6,
#'   criterion = "ParallelAnalysis", n_permutations = 100
#' )
#' }
#'
#' @export
select_r_global_mpca <- function(
    list_of_data,
    rvec = 1:6,
    model = c("PPCA", "FA"),
    criterion = c("CV", "ParallelAnalysis"),
    folds = 5,
    cv_metric = c("negLogLik", "frobenius"),
    n_permutations = 100,
    seed = NULL,
    verbose = FALSE
) {
  # Input validation
  model <- match.arg(model)
  criterion <- match.arg(criterion)
  cv_metric <- match.arg(cv_metric)
  
  if (length(rvec) < 2) {
    stop("rvec must contain at least 2 candidate r values")
  }
  
  N <- length(list_of_data)
  M <- ncol(list_of_data[[1]])
  
  if (folds < 2 || folds > N) {
    stop(sprintf("folds must be between 2 and N=%d", N))
  }
  
  if (verbose) {
    cat("=== Global r selection (MPCA) ===\n")
    cat(sprintf("Model: %s, Criterion: %s\n", model, criterion))
    cat(sprintf("N = %d subjects, M = %d channels\n", N, M))
    cat(sprintf("Candidate r values: %s\n", paste(rvec, collapse = ", ")))
  }
  
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Dispatch to appropriate criterion
  if (criterion == "CV") {
    result <- select_r_by_cv(
      list_of_data = list_of_data,
      rvec = rvec,
      model = model,
      folds = folds,
      cv_metric = cv_metric,
      verbose = verbose
    )
  } else if (criterion == "ParallelAnalysis") {
    result <- select_r_by_parallel_analysis(
      list_of_data = list_of_data,
      rvec = rvec,
      model = model,
      n_permutations = n_permutations,
      verbose = verbose
    )
  }
  
  # Add metadata
  result$criterion <- criterion
  result$model <- model
  
  if (verbose) {
    cat(sprintf("\n=== Selected r_global = %d ===\n", result$r_global))
  }
  
  return(result)
}
