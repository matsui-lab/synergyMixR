# Bootstrap Confidence Intervals for Mixture Model Loadings
#
# This file provides functions for computing bootstrap confidence intervals
# for synergy loadings from mixture factor analysis and mixture PCA models.

#' Bootstrap Confidence Intervals for Synergy Loadings
#'
#' Computes bootstrap confidence intervals for the loading matrices in mixture
#' PCA or mixture FA models. Uses Procrustes alignment to ensure bootstrap
#' samples are comparable across resamples.
#'
#' @param list_of_data A list of N matrices, each of dimension (T_i x M),
#'   containing the observed EMG data.
#' @param K Integer. Number of clusters.
#' @param r Integer. Number of synergies (factors/principal components).
#' @param B Integer. Number of bootstrap iterations. Default: 200.
#' @param alpha Numeric. Significance level for confidence intervals.
#'   Default: 0.05 (95\% CI).
#' @param mc_cores Integer. Number of cores for parallel computation.
#'   Default: 1 (sequential).
#' @param model_type Character. Either "MPCA" or "MFA". Default: "MPCA".
#' @param max_iter Integer. Maximum EM iterations per fit. Default: 50.
#' @param n_init Integer. Number of random initializations. Default: 1.
#' @param use_kmeans_init Logical. Whether to use k-means initialization.
#'   Default: TRUE.
#' @param seed Integer. Random seed for reproducibility. Default: NULL.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{loadings_median}{List of K matrices (M x r), median loadings for
#'     each cluster}
#'   \item{loadings_lower}{List of K matrices (M x r), lower CI bound for
#'     each cluster}
#'   \item{loadings_upper}{List of K matrices (M x r), upper CI bound for
#'     each cluster}
#'   \item{loadings_se}{List of K matrices (M x r), standard errors for
#'     each cluster}
#'   \item{bootstrap_samples}{Array of dimension (B x K x M x r) containing
#'     all aligned bootstrap samples (optional, if return_samples=TRUE)}
#'   \item{reference_fit}{The reference model fit used for alignment}
#'   \item{n_successful}{Number of successful bootstrap iterations}
#'   \item{alpha}{The significance level used}
#' }
#'
#' @details
#' The bootstrap procedure:
#' 1. Fit the model to the full dataset to obtain reference loadings
#' 2. For each bootstrap iteration b = 1, ..., B:
#'    a. Resample subjects with replacement
#'    b. Fit the model to the bootstrap sample
#'    c. Align bootstrap loadings to reference using Procrustes rotation
#' 3. Compute percentile confidence intervals from aligned samples
#'
#' Procrustes alignment is essential because factor loadings are only
#' identifiable up to rotation and sign. Without alignment, bootstrap samples
#' would not be comparable, leading to inflated variance estimates.
#'
#' @section Procrustes Alignment:
#' For each cluster k, the bootstrap loading matrix W_k(b) is aligned to
#' the reference W_k using orthogonal Procrustes:
#' W_aligned(b) = W_k(b) R*
#' where R* minimizes ||W_k - W_k(b) R||_F subject to R'R = I.
#'
#' The optimal R is found via SVD: if W_k' W_k(b) = UDV', then R* = VU'.
#'
#' @examples
#' \dontrun{
#' # Simulate some data
#' sim <- simulate_dynamic_synergy_data(N = 50, K = 2, r = 3, M = 8, T_each = 100)
#'
#' # Compute bootstrap CIs
#' boot_ci <- bootstrap_loading_ci(
#'   list_of_data = sim$X_list,
#'   K = 2,
#'   r = 3,
#'   B = 100,
#'   alpha = 0.05,
#'   mc_cores = 4
#' )
#'
#' # View median loadings for cluster 1
#' print(boot_ci$loadings_median[[1]])
#'
#' # View 95% CI width for first synergy, first muscle, cluster 1
#' width <- boot_ci$loadings_upper[[1]][1,1] - boot_ci$loadings_lower[[1]][1,1]
#' print(width)
#' }
#'
#' @seealso \code{\link{mixture_pca_em_fit}}, \code{\link{mfa_em_fit}}
#'
#' @export
bootstrap_loading_ci <- function(list_of_data,
                                  K,
                                  r,
                                  B = 200,
                                  alpha = 0.05,
                                  mc_cores = 1,
                                  model_type = c("MPCA", "MFA"),
                                  max_iter = 50,
                                  n_init = 1,
                                  use_kmeans_init = TRUE,
                                  seed = NULL,
                                  verbose = TRUE) {

  model_type <- match.arg(model_type)

  # Input validation
  validate_list_of_data(list_of_data)
  validate_model_params(K = K, r = r, max_iter = max_iter, n_init = n_init)

  N <- length(list_of_data)
  M <- ncol(list_of_data[[1]])

  validate_K_vs_N(K, N)
  validate_r_vs_M(r, M)

  if (!is.numeric(B) || B < 1) {
    stop("B must be a positive integer")
  }
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be between 0 and 1")
  }

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (verbose) {
    cat("=== Bootstrap Loading CIs ===\n")
    cat(sprintf("Model: %s, N=%d, K=%d, r=%d, M=%d\n", model_type, N, K, r, M))
    cat(sprintf("Bootstrap iterations: %d, alpha: %.2f\n", B, alpha))
    cat(sprintf("Parallel cores: %d\n\n", mc_cores))
  }

  # =========================================
  # Step 1: Fit reference model
  # =========================================
  if (verbose) cat("Fitting reference model...\n")

  fit_model <- function(data, verbose_fit = FALSE) {
    if (model_type == "MPCA") {
      mixture_pca_em_fit(
        list_of_data = data,
        K = K,
        r = r,
        max_iter = max_iter,
        n_init = n_init,
        use_kmeans_init = use_kmeans_init,
        verbose = verbose_fit
      )
    } else {
      mfa_em_fit(
        list_of_data = data,
        K = K,
        r = r,
        max_iter = max_iter,
        n_init = n_init,
        use_kmeans_init = use_kmeans_init,
        verbose = verbose_fit
      )
    }
  }

  reference_fit <- fit_model(list_of_data, verbose_fit = verbose)

  # Extract reference loadings
  if (model_type == "MPCA") {
    reference_W <- reference_fit$W
    if (is.null(reference_W)) {
      # Reconstruct W from P and D if needed
      reference_W <- lapply(seq_len(K), function(k) {
        reference_fit$P[[k]] %*% reference_fit$D[[k]]
      })
    }
  } else {
    reference_W <- reference_fit$Lambda
  }

  if (verbose) cat("Reference model fitted.\n\n")

  # =========================================
  # Step 2: Bootstrap iterations
  # =========================================
  if (verbose) cat(sprintf("Running %d bootstrap iterations...\n", B))

  # Pre-generate bootstrap indices for reproducibility
  boot_indices <- lapply(1:B, function(b) sample(N, N, replace = TRUE))

  # Function to run single bootstrap iteration
  run_one_bootstrap <- function(b) {
    # Resample subjects
    boot_idx <- boot_indices[[b]]
    boot_data <- list_of_data[boot_idx]

    # Fit model
    boot_fit <- tryCatch({
      fit_model(boot_data, verbose_fit = FALSE)
    }, error = function(e) NULL)

    if (is.null(boot_fit)) {
      return(NULL)
    }

    # Extract loadings
    if (model_type == "MPCA") {
      boot_W <- boot_fit$W
      if (is.null(boot_W)) {
        boot_W <- lapply(seq_len(K), function(k) {
          boot_fit$P[[k]] %*% boot_fit$D[[k]]
        })
      }
    } else {
      boot_W <- boot_fit$Lambda
    }

    # Handle cluster label switching
    # Match bootstrap clusters to reference clusters using Hungarian algorithm
    boot_W_aligned <- align_bootstrap_loadings(boot_W, reference_W, boot_fit$z)

    return(boot_W_aligned)
  }

  # Run bootstrap in parallel or sequentially
  if (mc_cores > 1 && .Platform$OS.type != "windows") {
    boot_results <- parallel::mclapply(
      1:B,
      run_one_bootstrap,
      mc.cores = mc_cores
    )
  } else {
    boot_results <- lapply(1:B, function(b) {
      if (verbose && b %% 20 == 0) {
        cat(sprintf("  Bootstrap %d/%d\n", b, B))
      }
      run_one_bootstrap(b)
    })
  }

  # Filter out failed iterations
  boot_results <- boot_results[!sapply(boot_results, is.null)]
  n_successful <- length(boot_results)

  if (verbose) {
    cat(sprintf("Completed %d/%d bootstrap iterations successfully.\n\n",
                n_successful, B))
  }

  if (n_successful < 10) {
    warning("Very few successful bootstrap iterations. Results may be unreliable.")
  }

  # =========================================
  # Step 3: Compute confidence intervals
  # =========================================
  if (verbose) cat("Computing confidence intervals...\n")

  # Convert boot_results to array for easier computation
  # boot_results is a list of length n_successful
  # Each element is a list of K matrices (M x r)

  loadings_median <- vector("list", K)
  loadings_lower <- vector("list", K)
  loadings_upper <- vector("list", K)
  loadings_se <- vector("list", K)

  for (k in 1:K) {
    # Extract all bootstrap samples for cluster k
    # boot_samples_k is an array of dimension (n_successful x M x r)
    boot_samples_k <- array(NA, dim = c(n_successful, M, r))

    for (b in 1:n_successful) {
      if (!is.null(boot_results[[b]][[k]])) {
        boot_samples_k[b, , ] <- boot_results[[b]][[k]]
      }
    }

    # Compute quantiles for each element
    loadings_median[[k]] <- matrix(NA, M, r)
    loadings_lower[[k]] <- matrix(NA, M, r)
    loadings_upper[[k]] <- matrix(NA, M, r)
    loadings_se[[k]] <- matrix(NA, M, r)

    for (i in 1:M) {
      for (j in 1:r) {
        vals <- boot_samples_k[, i, j]
        vals <- vals[!is.na(vals)]

        if (length(vals) > 0) {
          loadings_median[[k]][i, j] <- median(vals)
          loadings_lower[[k]][i, j] <- quantile(vals, alpha / 2)
          loadings_upper[[k]][i, j] <- quantile(vals, 1 - alpha / 2)
          loadings_se[[k]][i, j] <- sd(vals)
        }
      }
    }

    # Copy row/column names from reference
    if (!is.null(rownames(reference_W[[k]]))) {
      rownames(loadings_median[[k]]) <- rownames(reference_W[[k]])
      rownames(loadings_lower[[k]]) <- rownames(reference_W[[k]])
      rownames(loadings_upper[[k]]) <- rownames(reference_W[[k]])
      rownames(loadings_se[[k]]) <- rownames(reference_W[[k]])
    }
  }

  if (verbose) cat("Done.\n")

  # Return results
  list(
    loadings_median = loadings_median,
    loadings_lower = loadings_lower,
    loadings_upper = loadings_upper,
    loadings_se = loadings_se,
    reference_fit = reference_fit,
    reference_loadings = reference_W,
    n_successful = n_successful,
    B = B,
    alpha = alpha,
    K = K,
    r = r,
    model_type = model_type
  )
}


#' Align Bootstrap Loadings to Reference Using Procrustes
#'
#' Internal function that aligns bootstrap loading matrices to reference
#' loadings using orthogonal Procrustes rotation and handles cluster
#' label switching.
#'
#' @param boot_W List of K loading matrices from bootstrap sample
#' @param reference_W List of K loading matrices from reference fit
#' @param boot_z Bootstrap cluster assignments (for label switching detection)
#'
#' @return List of K aligned loading matrices
#'
#' @keywords internal
align_bootstrap_loadings <- function(boot_W, reference_W, boot_z) {
  K <- length(reference_W)
  r <- ncol(reference_W[[1]])

  # Step 1: Match clusters between bootstrap and reference
  # Use loading similarity to find best cluster matching
  similarity_matrix <- matrix(0, K, K)

  for (k_boot in 1:K) {
    for (k_ref in 1:K) {
      # Compute similarity as trace(reference' * bootstrap)
      # after normalizing columns
      W_ref <- reference_W[[k_ref]]
      W_boot <- boot_W[[k_boot]]

      # Normalize columns
      W_ref_norm <- apply(W_ref, 2, function(x) x / sqrt(sum(x^2) + 1e-10))
      W_boot_norm <- apply(W_boot, 2, function(x) x / sqrt(sum(x^2) + 1e-10))

      # Compute Frobenius norm of difference
      similarity_matrix[k_boot, k_ref] <- sum(abs(W_ref_norm * W_boot_norm))
    }
  }

  # Use Hungarian algorithm for optimal matching
  if (requireNamespace("clue", quietly = TRUE)) {
    assignment <- clue::solve_LSAP(similarity_matrix, maximum = TRUE)
    cluster_mapping <- as.integer(assignment)
  } else {
    # Greedy matching fallback
    cluster_mapping <- integer(K)
    used <- logical(K)
    for (k in 1:K) {
      available <- which(!used)
      best <- available[which.max(similarity_matrix[k, available])]
      cluster_mapping[k] <- best
      used[best] <- TRUE
    }
  }

  # Step 2: Align each cluster's loadings using Procrustes
  aligned_W <- vector("list", K)

  for (k_ref in 1:K) {
    k_boot <- which(cluster_mapping == k_ref)
    if (length(k_boot) == 0) {
      k_boot <- k_ref  # Fallback
    }

    W_ref <- reference_W[[k_ref]]
    W_boot <- boot_W[[k_boot]]

    # Orthogonal Procrustes: find R that minimizes ||W_ref - W_boot R||_F
    # Solution: R = V U' where W_ref' W_boot = U D V'
    svd_result <- svd(t(W_ref) %*% W_boot)
    R_opt <- svd_result$v %*% t(svd_result$u)

    # Apply rotation
    aligned_W[[k_ref]] <- W_boot %*% R_opt

    # Sign correction: ensure positive correlation with reference for each column
    for (j in 1:r) {
      if (cor(W_ref[, j], aligned_W[[k_ref]][, j]) < 0) {
        aligned_W[[k_ref]][, j] <- -aligned_W[[k_ref]][, j]
      }
    }
  }

  aligned_W
}


#' Print Summary of Bootstrap CI Results
#'
#' @param x Result from bootstrap_loading_ci
#' @param cluster Which cluster to summarize (default: 1)
#' @param ... Additional arguments (ignored)
#'
#' @export
print.bootstrap_ci <- function(x, cluster = 1, ...) {
  cat("Bootstrap Confidence Intervals for Synergy Loadings\n")
  cat("====================================================\n")
  cat(sprintf("Model: %s\n", x$model_type))
  cat(sprintf("Clusters (K): %d\n", x$K))
  cat(sprintf("Synergies (r): %d\n", x$r))
  cat(sprintf("Bootstrap iterations: %d/%d successful\n", x$n_successful, x$B))
  cat(sprintf("Confidence level: %d%%\n", round((1 - x$alpha) * 100)))
  cat("\n")
  cat(sprintf("Showing cluster %d:\n", cluster))
  cat("\nMedian Loadings:\n")
  print(round(x$loadings_median[[cluster]], 3))
  cat("\nStandard Errors:\n")
  print(round(x$loadings_se[[cluster]], 3))
  invisible(x)
}
