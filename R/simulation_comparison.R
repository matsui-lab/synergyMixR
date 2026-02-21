#' Simulation Comparison Functions
#'
#' Functions for comparing different methods (SingleFA, SinglePCA, MixtureFA, MixturePCA)
#' on simulated data.
#'
#' @name simulation_comparison
NULL

# ============================================================
# Single Factor Analysis (FA)
# ============================================================

#' Fit a Single Factor Analysis (Ignoring Clusters)
#'
#' This function concatenates all subjects' data into one large matrix,
#' then performs a naive factor analysis. It returns a simple BIC measure
#' and the fitted loadings, mean, and diagonal Psi.
#'
#' @param list_of_data A list of (T_i x M) matrices.
#' @param r The factor dimension.
#'
#' @return A list with elements \code{(loadings, mu, diagPsi, logLik, BIC)}.
#' @export
#' @examples
#' \dontrun{
#' sim <- generate_synergy_data(N = 50, K = 3, r = 2, M = 6, T_each = 100, seed = 123)
#' fa_res <- fit_single_factor_analysis(sim$list_of_data, r = 2)
#' print(fa_res$BIC)
#' }
fit_single_factor_analysis <- function(list_of_data, r = 2) {
  bigX <- do.call(rbind, list_of_data)
  n_total <- nrow(bigX)
  M <- ncol(bigX)

  mu <- colMeans(bigX)
  Xc <- sweep(bigX, 2, mu, "-")
  Sxx <- crossprod(Xc) / n_total

  eig_out <- eigen(Sxx, symmetric = TRUE)
  d_r  <- eig_out$values[1:r]
  V_r  <- eig_out$vectors[, 1:r, drop = FALSE]
  loadings <- V_r %*% diag(sqrt(d_r), r, r)

  recon <- loadings %*% t(loadings)
  diagPsi <- diag(Sxx) - diag(recon)
  diagPsi[diagPsi < 1e-12] <- 1e-12

  Sigma <- recon
  diag(Sigma) <- diag(Sigma) + diagPsi
  cholS <- chol(Sigma)
  logdetS <- 2 * sum(log(diag(cholS)))
  Mlog2pi <- M * log(2 * pi)
  invS <- solve(Sigma)
  quadSum <- 0
  for (i in seq_len(n_total)) {
    xi <- bigX[i, ]
    diff <- xi - mu
    quadSum <- quadSum + (t(diff) %*% invS %*% diff)
  }
  logLik_val <- -0.5 * (n_total * (Mlog2pi + logdetS) + quadSum)

  # Naive param count
  param_count <- M * r + M
  BIC_val <- -2 * logLik_val + param_count * log(n_total)

  list(
    loadings = loadings,
    mu = mu,
    diagPsi = diagPsi,
    logLik = logLik_val,
    BIC    = BIC_val
  )
}

#' Compute Total Reconstruction SSE for Single FA
#'
#' @param list_of_data The original data list.
#' @param fa_res A result from \code{fit_single_factor_analysis}.
#'
#' @return A numeric value of the sum of squared errors (SSE).
#' @export
calc_reconstruction_error_singleFA <- function(list_of_data, fa_res) {
  loadings <- fa_res$loadings
  mu <- fa_res$mu
  diagPsi <- fa_res$diagPsi
  SSE <- 0
  for (X_i in list_of_data) {
    Xhat_i <- reconstruct_singleFA(X_i, mu, loadings, diagPsi)
    diff <- X_i - Xhat_i
    SSE <- SSE + sum(diff^2)
  }
  SSE
}

# ============================================================
# Single PCA (Ignoring Clusters)
# ============================================================

#' Fit a Single PCA (Ignoring Clusters)
#'
#' Concatenates all data into one matrix, performs SVD, keeps \code{r}
#' principal components, then computes a naive PPCA logLik and BIC.
#'
#' @param list_of_data List of (T_i x M) matrices.
#' @param r Number of principal components.
#'
#' @return A list with \code{(P, mu, sigma2, logLik, BIC)}.
#' @export
#' @examples
#' \dontrun{
#' sim <- generate_synergy_data(N = 50, K = 3, r = 2, M = 6, T_each = 100, seed = 123)
#' pca_res <- fit_single_pca(sim$list_of_data, r = 2)
#' print(pca_res$BIC)
#' }
fit_single_pca <- function(list_of_data, r = 2) {
  bigX <- do.call(rbind, list_of_data)
  n_total <- nrow(bigX)
  M <- ncol(bigX)

  mu <- colMeans(bigX)
  Xc <- sweep(bigX, 2, mu, "-")

  svd_out <- svd(Xc)
  D_r <- svd_out$d[1:r]
  V_r <- svd_out$v[, 1:r, drop = FALSE]
  P <- V_r  # principal directions (M x r)

  lam_j <- (D_r^2) / n_total
  if (r < M) {
    leftover <- 0
    for (j in seq.int(r + 1, M)) {
      leftover <- leftover + (svd_out$d[j]^2) / n_total
    }
    sigma2 <- leftover / (M - r)
  } else {
    sigma2 <- 1e-12
  }

  # rank-r => Sigma = W W^T + sigma^2 I
  W <- matrix(0, nrow = M, ncol = r)
  for (j in seq_len(r)) {
    val_j <- lam_j[j] - sigma2
    if (val_j < 0) val_j <- 0
    W[, j] <- P[, j] * sqrt(val_j)
  }
  Sigma <- W %*% t(W)
  diag(Sigma) <- diag(Sigma) + sigma2

  cholS <- chol(Sigma)
  logdetS <- 2 * sum(log(diag(cholS)))
  Mlog2pi <- M * log(2 * pi)
  invS <- solve(Sigma)
  quadSum <- 0
  for (i in seq_len(n_total)) {
    xi <- bigX[i, ]
    diff <- xi - mu
    quadSum <- quadSum + (t(diff) %*% invS %*% diff)
  }
  logLik_val <- -0.5 * (n_total * (Mlog2pi + logdetS) + quadSum)

  # param count ~ (M*r) + M + 1
  param_count <- (M * r) + M + 1
  BIC_val <- -2 * logLik_val + param_count * log(n_total)

  list(
    P = P,
    mu = mu,
    sigma2 = sigma2,
    logLik = logLik_val,
    BIC    = BIC_val
  )
}

#' Compute Total Reconstruction SSE for Single PCA
#'
#' @param list_of_data The data list.
#' @param pca_res The result from \code{fit_single_pca}.
#'
#' @return SSE (sum of squared errors).
#' @export
calc_reconstruction_error_singlePCA <- function(list_of_data, pca_res) {
  SSE <- 0
  for (X_i in list_of_data) {
    Xhat_i <- reconstruct_singlePCA(X_i, pca_res)
    diff <- X_i - Xhat_i
    SSE <- SSE + sum(diff^2)
  }
  SSE
}

# ============================================================
# Mixture FA / Mixture PCA Reconstruction
# ============================================================

#' Compute SSE for a MixtureFA model
#'
#' Uses the final hard assignments \code{mfa_fit$z}, reconstruct each subject's data.
#'
#' @param list_of_data The data list.
#' @param mfa_fit The mixture FA result, which contains \code{$z}, \code{$Lambda[[k]]}, \code{$mu[[k]]}, \code{$Psi[[k]]}.
#'
#' @return SSE (sum of squared errors).
#' @export
calc_reconstruction_error_mixtureFA <- function(list_of_data, mfa_fit) {
  z_vec <- mfa_fit$z
  SSE <- 0
  for (i in seq_along(list_of_data)) {
    k_i <- z_vec[i]
    X_i <- list_of_data[[i]]
    Xhat_i <- reconstruct_mixtureFA_oneCluster(X_i, k_i, mfa_fit)
    SSE <- SSE + sum((X_i - Xhat_i)^2)
  }
  SSE
}

#' Compute SSE for a MixturePCA model
#'
#' @param list_of_data The data list.
#' @param mpca_fit The mixture PCA result, containing \code{$z} and each cluster's \code{P, mu}.
#'
#' @return SSE (sum of squared errors).
#' @export
calc_reconstruction_error_mixturePCA <- function(list_of_data, mpca_fit) {
  z_vec <- mpca_fit$z
  SSE <- 0
  for (i in seq_along(list_of_data)) {
    k_i <- z_vec[i]
    X_i <- list_of_data[[i]]
    Xhat_i <- reconstruct_mixturePCA_oneCluster(X_i, k_i, mpca_fit)
    SSE <- SSE + sum((X_i - Xhat_i)^2)
  }
  SSE
}

#' Compute the Total Sum of Squares (SST) for a List of Data Matrices
#'
#' This function concatenates all matrices in \code{list_of_data} into one big matrix,
#' calculates the global mean, and then sums the squared differences from that mean.
#' The result can be used to compute VAF (Variance Accounted For).
#'
#' @param list_of_data A list of \code{(T_i x M)} matrices.
#' @return A numeric scalar, the total sum of squares (SST).
#' @export
calc_total_SST <- function(list_of_data) {
  # 1) Concatenate everything
  bigX <- do.call(rbind, list_of_data)  # shape: (sum T_i) x M

  # 2) Global mean across all rows
  mu_global <- colMeans(bigX)

  # 3) Sum of squared diffs
  Xc <- sweep(bigX, 2, mu_global, "-")  # center
  sst <- sum(Xc^2)  # sum of squares
  sst
}

# ============================================================
# Two-step Baseline Algorithm
# ============================================================

#' Two-step Baseline Algorithm for PCA/FA
#'
#' This function implements a two-step baseline approach for comparison with
#' mixture-based methods:
#' \enumerate{
#'   \item \strong{Step 1: Individual Decomposition} - Apply PCA or FA to each
#'         subject's data separately to extract \code{r} basis vectors per subject.
#'   \item \strong{Step 1.5: Basis Alignment} - Align basis vectors across subjects
#'         using the Hungarian algorithm to ensure corresponding components are in
#'         the same order, improving clustering accuracy.
#'   \item \strong{Step 2: Clustering} - Normalize the basis vectors, concatenate
#'         them across all subjects, and apply k-means clustering. Assign each
#'         subject to a cluster based on majority vote of their basis vectors.
#' }
#'
#' @param data_list A list of \code{(T_i x M)} matrices, one per subject.
#' @param r Number of components (factors or principal components) to extract per subject.
#' @param K Number of clusters for k-means.
#' @param method Either \code{"PCA"} or \code{"FA"} for the decomposition method.
#' @param seed Random seed for reproducibility.
#' @param z_true Optional true cluster labels for computing ARI. If \code{NULL}, ARI is set to \code{NA}.
#' @param refine_assignments Logical. If \code{TRUE}, performs a single reassignment
#'   pass after fitting cluster-specific models in Step 4. Default is \code{FALSE}.
#' @param align_basis Logical. If \code{TRUE} (default), applies Hungarian algorithm
#'   alignment in Step 1.5 to ensure corresponding components across subjects are in
#'   the same order before clustering. Set to \code{FALSE} to disable alignment.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{z_est}}{Integer vector of estimated cluster assignments (1..K).}
#'   \item{\code{W_list}}{List of normalized basis matrices, one per subject.}
#'   \item{\code{W_centroids}}{Matrix of cluster centroids in basis space.}
#'   \item{\code{mu_list}}{List of centering vectors used for standardization, one per subject.}
#'   \item{\code{sd_list}}{List of scaling vectors used for standardization, one per subject.}
#'   \item{\code{cluster_models}}{List of cluster-specific models (each with \code{W} and \code{n_members}).}
#'   \item{\code{method}}{Character string indicating the method used (e.g., "TwoStep_PCA").}
#'   \item{\code{BIC}}{Set to \code{NA} (k-means does not have a canonical BIC).}
#'   \item{\code{ARI}}{Adjusted Rand Index if \code{z_true} is provided, otherwise \code{NA}.}
#' }
#'
#' @details
#' The algorithm addresses sign indeterminacy by normalizing each basis vector
#' so that its maximum absolute value element is positive. Tie-breaking in
#' majority vote is handled by computing the mean distance from a subject's
#' basis vectors to each tied cluster centroid and selecting the closer one.
#'
#' Step 1.5 (basis alignment) uses the Hungarian algorithm to align basis vectors
#' across subjects. This ensures that corresponding components (e.g., the first
#' component of each subject) represent similar patterns, which is crucial for
#' accurate clustering. Without alignment, the arbitrary ordering of components
#' from individual PCA/FA can lead to poor clustering performance.
#'
#' Step 4 (re-optimization) fits cluster-specific models in standardized space
#' to ensure consistency with Step 1. This improves reconstruction metrics (SSE/VAF)
#' but does not change ARI unless \code{refine_assignments = TRUE}.
#'
#' If \code{refine_assignments = TRUE}, a single reassignment pass is performed
#' after fitting cluster-specific models. Each subject is reassigned to the cluster
#' that minimizes reconstruction error. This can improve ARI while maintaining the
#' two-step nature of the baseline (unlike iterative EM which would drift toward
#' the mixture model).
#'
#' For FA decomposition, the function uses \code{psych::fa()} if available.
#' If the \code{psych} package is not installed or if FA fails for a subject,
#' the function falls back to PCA with a warning.
#'
#' @examples
#' \dontrun{
#' sim <- generate_synergy_data(N = 50, K = 3, r = 2, M = 6, T_each = 100, seed = 123)
#' result <- fit_two_step_baseline(
#'   sim$list_of_data, r = 2, K = 3, method = "PCA",
#'   seed = 123, z_true = sim$true_cluster
#' )
#' print(result$ARI)
#' }
#'
#' @export
fit_two_step_baseline <- function(data_list, r, K, method = c("PCA", "FA"),
                                   seed = 123, z_true = NULL, refine_assignments = FALSE,
                                   align_basis = TRUE) {
  set.seed(seed)
  method <- match.arg(method)
  N <- length(data_list)
  M <- ncol(data_list[[1]])
  W_list <- vector("list", N)
  mu_list <- vector("list", N)
  sd_list <- vector("list", N)
  
  # Helper function to normalize basis vectors with sign alignment
  normalize_basis <- function(W) {
    W_norm <- apply(W, 2, function(col) {
      # L2 normalize
      col_norm <- col / sqrt(sum(col^2))
      # Align sign: make max absolute value element positive
      max_idx <- which.max(abs(col_norm))
      if (col_norm[max_idx] < 0) {
        col_norm <- -col_norm
      }
      col_norm
    })
    W_norm
  }
  
  # Step 1: Individual decomposition
  for (i in seq_len(N)) {
    Xi <- data_list[[i]]
    # Center and scale
    Xi_scaled <- scale(Xi, center = TRUE, scale = TRUE)
    
    # Store centering and scaling parameters for reconstruction
    mu_i <- attr(Xi_scaled, "scaled:center")
    sd_i <- attr(Xi_scaled, "scaled:scale")
    
    # Guard against zero-variance columns
    if (any(sd_i <= 0 | is.na(sd_i))) {
      zero_var_idx <- which(sd_i <= 0 | is.na(sd_i))
      warning(sprintf("Subject %d has zero-variance columns: %s. Setting sd=1 for these columns.",
                      i, paste(zero_var_idx, collapse=", ")))
      sd_i[sd_i <= 0 | is.na(sd_i)] <- 1
    }
    
    mu_list[[i]] <- mu_i
    sd_list[[i]] <- sd_i
    
    if (method == "PCA") {
      # Use prcomp for numerical stability
      pca_res <- tryCatch({
        prcomp(Xi_scaled, center = FALSE, scale. = FALSE)
      }, error = function(e) {
        warning(sprintf("PCA failed for subject %d: %s. Using fallback.", i, e$message))
        NULL
      })
      
      if (!is.null(pca_res)) {
        Wi <- pca_res$rotation[, 1:r, drop = FALSE]
      } else {
        # Fallback: use identity-based basis
        Wi <- diag(M)[, 1:min(r, M), drop = FALSE]
        if (r > M) {
          Wi <- cbind(Wi, matrix(0, nrow = M, ncol = r - M))
        }
      }
    } else {
      # method == "FA"
      if (!requireNamespace("psych", quietly = TRUE)) {
        warning("Package 'psych' not available. Falling back to PCA for FA method.")
        pca_res <- tryCatch({
          prcomp(Xi_scaled, center = FALSE, scale. = FALSE)
        }, error = function(e) {
          warning(sprintf("PCA fallback failed for subject %d: %s", i, e$message))
          NULL
        })
        
        if (!is.null(pca_res)) {
          Wi <- pca_res$rotation[, 1:r, drop = FALSE]
        } else {
          Wi <- diag(M)[, 1:min(r, M), drop = FALSE]
          if (r > M) {
            Wi <- cbind(Wi, matrix(0, nrow = M, ncol = r - M))
          }
        }
      } else {
        fa_res <- tryCatch({
          psych::fa(Xi_scaled, nfactors = r, rotate = "none", fm = "ml")
        }, error = function(e) {
          warning(sprintf("FA failed for subject %d: %s. Falling back to PCA.", i, e$message))
          NULL
        })
        
        if (!is.null(fa_res)) {
          Wi <- as.matrix(fa_res$loadings[, 1:r, drop = FALSE])
        } else {
          # Fallback to PCA
          pca_res <- tryCatch({
            prcomp(Xi_scaled, center = FALSE, scale. = FALSE)
          }, error = function(e) {
            warning(sprintf("PCA fallback failed for subject %d: %s", i, e$message))
            NULL
          })
          
          if (!is.null(pca_res)) {
            Wi <- pca_res$rotation[, 1:r, drop = FALSE]
          } else {
            Wi <- diag(M)[, 1:min(r, M), drop = FALSE]
            if (r > M) {
              Wi <- cbind(Wi, matrix(0, nrow = M, ncol = r - M))
            }
          }
        }
      }
    }
    
    # Normalize with sign alignment
    Wi <- normalize_basis(Wi)
    W_list[[i]] <- Wi
  }
  
  # Step 1.5: Basis Alignment (optional, controlled by align_basis parameter)
  # Align basis vectors across subjects using Hungarian algorithm
  # This ensures that corresponding components are in the same order across subjects,
  # which improves clustering accuracy
  if (align_basis && N > 1) {
    W_list <- align_basis_across_subjects(W_list, reference = 1, method = "correlation", verbose = FALSE)
    message("Performed basis alignment using Hungarian matching")
  } else if (!align_basis && N > 1) {
    message("Skipping basis alignment (align_basis = FALSE)")
  }
  
  # Step 2: Clustering
  # Concatenate all basis vectors: each subject contributes r vectors
  W_all <- do.call(rbind, lapply(W_list, t))  # (N*r) x M
  
  # Run k-means
  km <- tryCatch({
    stats::kmeans(W_all, centers = K, nstart = 20)
  }, error = function(e) {
    stop(sprintf("k-means clustering failed: %s", e$message))
  })
  
  labels <- km$cluster
  centroids <- km$centers  # K x M
  
  # Assign subject clusters by majority vote with tie-breaking
  z_est <- integer(N)
  for (i in seq_len(N)) {
    seg <- ((i - 1) * r + 1):(i * r)
    subject_labels <- labels[seg]
    
    # Count votes
    vote_counts <- table(subject_labels)
    max_count <- max(vote_counts)
    tied_clusters <- as.integer(names(vote_counts[vote_counts == max_count]))
    
    if (length(tied_clusters) == 1) {
      z_est[i] <- tied_clusters[1]
    } else {
      # Tie-breaking: compute mean distance to each tied centroid
      Wi <- W_list[[i]]  # M x r
      mean_dists <- sapply(tied_clusters, function(k) {
        centroid_k <- centroids[k, ]
        dists <- apply(Wi, 2, function(col) {
          sqrt(sum((col - centroid_k)^2))
        })
        mean(dists)
      })
      # Choose cluster with smallest mean distance
      min_idx <- which.min(mean_dists)
      if (length(min_idx) > 1) {
        # Still tied, choose smallest cluster label
        z_est[i] <- min(tied_clusters[min_idx])
      } else {
        z_est[i] <- tied_clusters[min_idx]
      }
    }
  }
  
  # Step 4: Re-optimize FA/PCA with fixed cluster assignments
  # Fit cluster-specific models in standardized space for consistency with Step 1
  cluster_models <- vector("list", K)
  for (k in 1:K) {
    cluster_idx <- which(z_est == k)
    if (length(cluster_idx) > 0) {
      # Stack standardized data for this cluster
      cluster_data_standardized <- list()
      for (idx in cluster_idx) {
        Xi <- data_list[[idx]]
        mu_i <- mu_list[[idx]]
        sd_i <- sd_list[[idx]]
        Xi_centered <- sweep(Xi, 2, mu_i, "-")
        Xi_standardized <- sweep(Xi_centered, 2, sd_i, "/")
        cluster_data_standardized[[length(cluster_data_standardized) + 1]] <- Xi_standardized
      }
      
      # Concatenate cluster data
      cluster_bigX <- do.call(rbind, cluster_data_standardized)
      
      # Fit model on standardized data
      if (method == "PCA") {
        # Use prcomp for numerical stability
        pca_res <- tryCatch({
          prcomp(cluster_bigX, center = FALSE, scale. = FALSE)
        }, error = function(e) {
          warning(sprintf("PCA failed for cluster %d: %s. Using fallback.", k, e$message))
          NULL
        })
        
        if (!is.null(pca_res)) {
          W_k <- pca_res$rotation[, 1:r, drop = FALSE]
        } else {
          # Fallback: use identity-based basis
          W_k <- diag(M)[, 1:min(r, M), drop = FALSE]
          if (r > M) {
            W_k <- cbind(W_k, matrix(0, nrow = M, ncol = r - M))
          }
        }
      } else {
        # method == "FA"
        if (!requireNamespace("psych", quietly = TRUE)) {
          warning("Package 'psych' not available. Falling back to PCA for FA method.")
          pca_res <- tryCatch({
            prcomp(cluster_bigX, center = FALSE, scale. = FALSE)
          }, error = function(e) {
            warning(sprintf("PCA fallback failed for cluster %d: %s", k, e$message))
            NULL
          })
          
          if (!is.null(pca_res)) {
            W_k <- pca_res$rotation[, 1:r, drop = FALSE]
          } else {
            W_k <- diag(M)[, 1:min(r, M), drop = FALSE]
            if (r > M) {
              W_k <- cbind(W_k, matrix(0, nrow = M, ncol = r - M))
            }
          }
        } else {
          fa_res <- tryCatch({
            psych::fa(cluster_bigX, nfactors = r, rotate = "none", fm = "ml")
          }, error = function(e) {
            warning(sprintf("FA failed for cluster %d: %s. Falling back to PCA.", k, e$message))
            NULL
          })
          
          if (!is.null(fa_res)) {
            W_k <- as.matrix(fa_res$loadings[, 1:r, drop = FALSE])
          } else {
            # Fallback to PCA
            pca_res <- tryCatch({
              prcomp(cluster_bigX, center = FALSE, scale. = FALSE)
            }, error = function(e) {
              warning(sprintf("PCA fallback failed for cluster %d: %s", k, e$message))
              NULL
            })
            
            if (!is.null(pca_res)) {
              W_k <- pca_res$rotation[, 1:r, drop = FALSE]
            } else {
              W_k <- diag(M)[, 1:min(r, M), drop = FALSE]
              if (r > M) {
                W_k <- cbind(W_k, matrix(0, nrow = M, ncol = r - M))
              }
            }
          }
        }
      }
      
      # Normalize with sign alignment (same as Step 1)
      W_k <- normalize_basis(W_k)
      cluster_models[[k]] <- list(W = W_k, n_members = length(cluster_idx))
    } else {
      cluster_models[[k]] <- NULL  # Empty cluster
    }
  }
  
  # Step 5 (optional): Reassign subjects based on reconstruction error
  if (refine_assignments) {
    z_refined <- integer(N)
    for (i in seq_len(N)) {
      Xi <- data_list[[i]]
      mu_i <- mu_list[[i]]
      sd_i <- sd_list[[i]]
      best_k <- z_est[i]  # Default to current assignment
      # Guard against invalid default
      if (is.na(best_k) || best_k < 1 || best_k > K) {
        best_k <- 1L
      }
      best_sse <- Inf
      
      for (k in 1:K) {
        if (!is.null(cluster_models[[k]])) {
          # Compute reconstruction SSE for this subject with cluster k's model
          W_k <- cluster_models[[k]]$W
          sse_k <- compute_subject_sse_with_basis(Xi, W_k, mu_i, sd_i)
          if (sse_k < best_sse) {
            best_sse <- sse_k
            best_k <- k
          }
        }
      }
      z_refined[i] <- best_k
    }
    z_est <- z_refined
  }
  
  # Ensure z_est is always valid (integer, in range 1:K)
  storage.mode(z_est) <- "integer"
  z_est[is.na(z_est) | z_est < 1 | z_est > K] <- 1L
  
  # Compute ARI if z_true is provided
  ARI_val <- NA
  if (!is.null(z_true)) {
    if (!requireNamespace("mclust", quietly = TRUE)) {
      warning("Package 'mclust' not available. Cannot compute ARI.")
    } else {
      ARI_val <- mclust::adjustedRandIndex(z_est, z_true)
    }
  }
  
  list(
    z_est = z_est,
    W_list = W_list,
    W_centroids = t(centroids),  # M x K
    mu_list = mu_list,
    sd_list = sd_list,
    cluster_models = cluster_models,
    method = paste0("TwoStep_", method),
    BIC = NA,
    ARI = ARI_val
  )
}

#' Helper function to compute SSE for a single subject with a given basis
#'
#' @param Xi Data matrix for subject i (T x M)
#' @param Wi Basis matrix (M x r)
#' @param mu_i Mean vector for subject i (length M)
#' @param sd_i Standard deviation vector for subject i (length M)
#' @return SSE value for this subject
#' @keywords internal
compute_subject_sse_with_basis <- function(Xi, Wi, mu_i, sd_i) {
  # Standardize using the same parameters from fitting
  Xi_centered <- sweep(Xi, 2, mu_i, "-")
  Xi_standardized <- sweep(Xi_centered, 2, sd_i, "/")
  
  # Project in standardized space
  H_i <- Xi_standardized %*% Wi  # T x r
  
  # Reconstruct in standardized space
  Xi_standardized_hat <- H_i %*% t(Wi)  # T x M
  
  # Unstandardize back to original units
  Xi_hat <- sweep(Xi_standardized_hat, 2, sd_i, "*")
  Xi_hat <- sweep(Xi_hat, 2, mu_i, "+")
  
  # Compute SSE in original units
  sum((Xi - Xi_hat)^2)
}

#' Helper function to compute SSE for Two-step baseline
#'
#' @param data_list List of data matrices
#' @param two_step_res Result from fit_two_step_baseline
#' @return SSE value
#' @keywords internal
calc_reconstruction_error_twostep <- function(data_list, two_step_res) {
  z_est <- two_step_res$z_est
  W_list <- two_step_res$W_list
  cluster_models <- two_step_res$cluster_models
  mu_list <- two_step_res$mu_list
  sd_list <- two_step_res$sd_list
  SSE <- 0
  
  for (i in seq_along(data_list)) {
    Xi <- data_list[[i]]
    k <- z_est[i]
    
    # Use cluster-specific model if available, otherwise fall back to individual W_i
    # Guard against invalid k values (NA, 0, or out of bounds)
    K <- length(cluster_models)
    use_cluster_model <- !is.null(cluster_models) && 
                         !is.na(k) && 
                         k >= 1 && 
                         k <= K && 
                         !is.null(cluster_models[[k]])
    
    if (use_cluster_model) {
      Wi <- cluster_models[[k]]$W  # Use cluster basis
    } else {
      Wi <- W_list[[i]]  # Fall back to individual basis
    }
    
    mu_i <- mu_list[[i]]
    sd_i <- sd_list[[i]]
    
    # Use the helper function for consistent SSE computation
    SSE <- SSE + compute_subject_sse_with_basis(Xi, Wi, mu_i, sd_i)
  }
  SSE
}
