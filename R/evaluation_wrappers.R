#' Extract Cluster Labels from a Fitted Model
#'
#' This function extracts cluster assignment labels from a fitted MFA or MPCA
#' model object. It searches through common field names used by different
#' model implementations.
#'
#' @param fit A fitted model object (MFA or MPCA).
#'
#' @return An integer vector of cluster assignments (1 to K), or NA if
#'   cluster labels cannot be extracted.
#'
#' @details
#' The function searches for cluster labels in the following fields (in order):
#' \itemize{
#'   \item \code{z}: Standard field for MFA/MPCA models
#'   \item \code{cluster}: Alternative naming
#'   \item \code{clusters}: Plural form
#'   \item \code{class}: Classification result
#'   \item \code{z_hat}: Estimated assignments
#'   \item \code{cl}: Short form
#'   \item \code{hard_cluster}: Hard assignment from soft clustering
#'   \item \code{labels}: Generic label field
#' }
#'
#' If the model has a \code{resp} (responsibility) matrix but no hard assignments,
#' the function computes hard assignments using MAP (maximum a posteriori).
#'
#' @examples
#' \dontrun{
#' fit <- mfa_em_fit(list_of_data, K = 3, r = 2)
#' pred_cluster <- extract_cluster_labels(fit)
#' }
#'
#' @export
extract_cluster_labels <- function(fit) {
  if (is.null(fit)) {
    return(NA)
  }
  
  # List of candidate field names for cluster assignments
  candidate_fields <- c("z", "cluster", "clusters", "class", "z_hat", 
                        "cl", "hard_cluster", "labels")
  
  # Search for cluster labels
  for (field in candidate_fields) {
    if (!is.null(fit[[field]])) {
      labels <- fit[[field]]
      if (is.numeric(labels) || is.integer(labels)) {
        return(as.integer(labels))
      }
    }
  }
  
  # If no hard assignments found, try to compute from responsibilities
  if (!is.null(fit$resp)) {
    resp <- fit$resp
    if (is.matrix(resp) && nrow(resp) > 0 && ncol(resp) > 0) {
      # MAP assignment: argmax over clusters for each subject
      labels <- apply(resp, 1, which.max)
      return(as.integer(labels))
    }
  }
  
  # Could not extract cluster labels
  return(NA)
}


#' Compute Adjusted Rand Index (ARI)
#'
#' Computes the Adjusted Rand Index between true and predicted cluster
#' assignments. This is a NA-safe wrapper around \code{mclust::adjustedRandIndex}.
#'
#' @param true_cluster Integer vector of true cluster assignments.
#' @param pred_cluster Integer vector of predicted cluster assignments.
#'
#' @return Numeric scalar: ARI value in [-1, 1], where 1 indicates perfect
#'   agreement. Returns NA if inputs are invalid or mclust is not available.
#'
#' @details
#' The Adjusted Rand Index (ARI) measures the similarity between two clusterings,
#' adjusted for chance. It is commonly used to evaluate clustering accuracy
#' when ground truth is available.
#'
#' @examples
#' \dontrun{
#' true_cluster <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
#' pred_cluster <- c(1, 1, 2, 2, 2, 2, 3, 3, 3)
#' ari <- compute_ari(true_cluster, pred_cluster)
#' }
#'
#' @export
compute_ari <- function(true_cluster, pred_cluster) {
  # Handle NA inputs
  if (any(is.na(true_cluster)) || any(is.na(pred_cluster))) {
    return(NA_real_)
  }
  
  # Check lengths match
  if (length(true_cluster) != length(pred_cluster)) {
    warning("Length mismatch between true_cluster and pred_cluster")
    return(NA_real_)
  }
  
  # Check for empty inputs
  if (length(true_cluster) == 0) {
    return(NA_real_)
  }
  
  # Use mclust::adjustedRandIndex if available
  if (requireNamespace("mclust", quietly = TRUE)) {
    tryCatch({
      ari <- mclust::adjustedRandIndex(true_cluster, pred_cluster)
      return(ari)
    }, error = function(e) {
      warning("Error computing ARI: ", e$message)
      return(NA_real_)
    })
  } else {
    warning("Package 'mclust' is required for ARI computation")
    return(NA_real_)
  }
}


#' Compute VAF (Variance Accounted For) for a Fitted Model
#'
#' Computes the global Variance Accounted For (VAF) for a fitted MFA or MPCA
#' model. This is a unified wrapper that dispatches to the appropriate
#' method-specific function.
#'
#' @param method Model type: "MFA" or "MPCA".
#' @param list_of_data A list of data matrices, each (T_i x M).
#' @param fit A fitted model object.
#'
#' @return Numeric scalar: VAF value in [0, 1], where 1 indicates perfect
#'   reconstruction. Returns NA if computation fails.
#'
#' @details
#' VAF is computed as \eqn{1 - SSE/SST}, where SSE is the sum of squared
#' reconstruction errors and SST is the total sum of squares.
#'
#' @examples
#' \dontrun{
#' fit <- mfa_em_fit(list_of_data, K = 3, r = 2)
#' vaf <- compute_vaf_wrapper("MFA", list_of_data, fit)
#' }
#'
#' @export
compute_vaf_wrapper <- function(method, list_of_data, fit) {
  if (is.null(fit)) {
    return(NA_real_)
  }
  
  tryCatch({
    if (method == "MFA") {
      vaf <- compute_global_vaf_mfa(list_of_data, fit)
    } else if (method == "MPCA") {
      vaf <- compute_global_vaf_mpca(list_of_data, fit)
    } else {
      warning("Unknown method: ", method)
      return(NA_real_)
    }
    return(vaf)
  }, error = function(e) {
    warning("Error computing VAF: ", e$message)
    return(NA_real_)
  })
}


#' Compute BIC for a Fitted Model
#'
#' Computes the Bayesian Information Criterion (BIC) for a fitted MFA or MPCA
#' model. This is a unified wrapper that dispatches to the appropriate
#' method-specific function.
#'
#' @param method Model type: "MFA" or "MPCA".
#' @param fit A fitted model object.
#' @param K Number of clusters.
#' @param r Factor dimension.
#' @param M Number of observed dimensions (muscles).
#' @param N_obs Number of observations (typically sum of T_i).
#'
#' @return Numeric scalar: BIC value. Lower is better. Returns NA if
#'   computation fails.
#'
#' @details
#' BIC is computed as \eqn{-2 \log L + p \log(N)}, where L is the likelihood,
#' p is the number of parameters, and N is the sample size.
#'
#' @examples
#' \dontrun{
#' fit <- mfa_em_fit(list_of_data, K = 3, r = 2)
#' N_obs <- sum(sapply(list_of_data, nrow))
#' M <- ncol(list_of_data[[1]])
#' bic <- compute_bic_wrapper("MFA", fit, K = 3, r = 2, M = M, N_obs = N_obs)
#' }
#'
#' @export
compute_bic_wrapper <- function(method, fit, K, r, M, N_obs) {
  if (is.null(fit)) {
    return(NA_real_)
  }
  
  tryCatch({
    loglik <- fit$logLik
    if (is.null(loglik) || !is.finite(loglik)) {
      return(NA_real_)
    }
    
    if (method == "MFA") {
      bic <- compute_BIC_mfa(loglik, K, r, M, N_obs)
    } else if (method == "MPCA") {
      bic <- compute_BIC_mpca(loglik, K, r, M, N_obs)
    } else {
      warning("Unknown method: ", method)
      return(NA_real_)
    }
    return(bic)
  }, error = function(e) {
    warning("Error computing BIC: ", e$message)
    return(NA_real_)
  })
}


#' Evaluate Model Selection Results
#'
#' Comprehensive evaluation of model selection results, computing all relevant
#' metrics including ARI, VAF, BIC, logLik, and selection accuracy.
#'
#' @param list_of_data A list of data matrices, each (T_i x M).
#' @param true_cluster Integer vector of true cluster assignments.
#' @param K_true True number of clusters.
#' @param r_true True factor dimension.
#' @param select_result Result from \code{\link{select_K_r}}.
#' @param runtime_sec Runtime in seconds (optional).
#'
#' @return A list with evaluation metrics:
#' \describe{
#'   \item{\code{K_hat}}{Estimated number of clusters.}
#'   \item{\code{r_hat}}{Estimated factor dimension.}
#'   \item{\code{K_correct}}{Logical: K_hat == K_true.}
#'   \item{\code{r_correct}}{Logical: r_hat == r_true.}
#'   \item{\code{Kr_correct}}{Logical: both K and r correct.}
#'   \item{\code{ARI}}{Adjusted Rand Index.}
#'   \item{\code{VAF}}{Variance Accounted For.}
#'   \item{\code{logLik}}{Log-likelihood.}
#'   \item{\code{BIC}}{Bayesian Information Criterion.}
#'   \item{\code{runtime_sec}}{Runtime in seconds.}
#' }
#'
#' @examples
#' \dontrun{
#' # Generate data
#' data <- generate_synergy_data(N = 30, K = 3, r = 3, seed = 123)
#'
#' # Run model selection
#' start_time <- Sys.time()
#' select_result <- select_K_r(data$list_of_data, method = "MFA", strategy = "grid")
#' runtime <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
#'
#' # Evaluate
#' eval_result <- evaluate_model_selection(
#'   data$list_of_data, data$true_cluster,
#'   K_true = 3, r_true = 3,
#'   select_result, runtime_sec = runtime
#' )
#' }
#'
#' @export
evaluate_model_selection <- function(
    list_of_data,
    true_cluster,
    K_true,
    r_true,
    select_result,
    runtime_sec = NA_real_
) {
  # Extract selection results
  K_hat <- select_result$best_K
  r_hat <- select_result$best_r
  best_model <- select_result$best_model
  method <- select_result$meta$method
  
  # Selection accuracy
  K_correct <- !is.na(K_hat) && K_hat == K_true
  r_correct <- !is.na(r_hat) && r_hat == r_true
  Kr_correct <- K_correct && r_correct
  
  # Extract predicted cluster labels
  pred_cluster <- extract_cluster_labels(best_model)
  
  # Compute ARI
  ari <- compute_ari(true_cluster, pred_cluster)
  
  # Compute VAF
  vaf <- compute_vaf_wrapper(method, list_of_data, best_model)
  
  # Extract logLik
  loglik <- if (!is.null(best_model$logLik)) best_model$logLik else NA_real_
  
  # Compute BIC
  M <- ncol(list_of_data[[1]])
  N_obs <- sum(sapply(list_of_data, nrow))
  bic <- compute_bic_wrapper(method, best_model, K_hat, r_hat, M, N_obs)
  
  list(
    K_hat = K_hat,
    r_hat = r_hat,
    K_correct = K_correct,
    r_correct = r_correct,
    Kr_correct = Kr_correct,
    ARI = ari,
    VAF = vaf,
    logLik = loglik,
    BIC = bic,
    runtime_sec = runtime_sec
  )
}


# =============================================================================
# Subspace Recovery Metrics (Step 4)
# =============================================================================

#' Compute Subspace Similarity Between Two Basis Matrices
#'
#' Computes the similarity between two subspaces spanned by basis matrices
#' using QR decomposition and SVD. This metric measures how well the estimated
#' subspace recovers the true subspace, independent of rotation and scaling.
#'
#' @param L_true A (M x r_true) matrix representing the true basis.
#' @param L_hat A (M x r_hat) matrix representing the estimated basis.
#'
#' @return Numeric scalar in [0, 1], where 1 indicates identical subspaces.
#'   Returns NA if inputs are invalid (NULL, empty, or contain NA/NaN/Inf).
#'
#' @details
#' The subspace similarity is computed as follows:
#' \enumerate{
#'   \item Orthonormalize both bases using QR decomposition
#'   \item Compute the SVD of Q_true' * Q_hat
#'   \item Return the mean of squared singular values (up to min_dim)
#' }
#'
#' This metric is invariant to rotation and scaling of the basis vectors,
#' making it suitable for comparing factor analysis solutions where
#' rotation indeterminacy exists.
#'
#' When r_true != r_hat, the metric uses min(r_true, r_hat) dimensions,
#' which allows comparison even when the estimated rank differs from truth.
#'
#' @examples
#' \dontrun{
#' # Same subspace (rotated)
#' L_true <- matrix(c(1, 0, 0, 1), 2, 2)
#' R <- matrix(c(cos(pi/4), -sin(pi/4), sin(pi/4), cos(pi/4)), 2, 2)
#' L_hat <- L_true %*% R
#' compute_subspace_similarity(L_true, L_hat)  # Should be ~1
#'
#' # Orthogonal subspaces
#' L_true <- matrix(c(1, 0, 0, 0), 2, 2)
#' L_hat <- matrix(c(0, 0, 0, 1), 2, 2)
#' compute_subspace_similarity(L_true, L_hat)  # Should be ~0
#' }
#'
#' @export
compute_subspace_similarity <- function(L_true, L_hat) {
  # NA-safe input validation
  if (is.null(L_true) || is.null(L_hat)) {
    return(NA_real_)
  }
  
  if (!is.matrix(L_true) || !is.matrix(L_hat)) {
    return(NA_real_)
  }
  
  if (nrow(L_true) == 0 || ncol(L_true) == 0 ||
      nrow(L_hat) == 0 || ncol(L_hat) == 0) {
    return(NA_real_)
  }
  
  if (anyNA(L_true) || anyNA(L_hat) ||
      any(is.infinite(L_true)) || any(is.infinite(L_hat)) ||
      any(is.nan(L_true)) || any(is.nan(L_hat))) {
    return(NA_real_)
  }
  
  # Check dimension compatibility (must have same number of rows = M)
  if (nrow(L_true) != nrow(L_hat)) {
    warning("Dimension mismatch: L_true has ", nrow(L_true), 
            " rows, L_hat has ", nrow(L_hat), " rows")
    return(NA_real_)
  }
  
  tryCatch({
    # Orthonormalize both bases using QR decomposition
    qr_true <- qr(L_true)
    qr_hat <- qr(L_hat)
    
    r_true <- ncol(L_true)
    r_hat <- ncol(L_hat)
    
    # Extract Q matrices (orthonormal bases)
    Q_true <- qr.Q(qr_true)
    Q_hat <- qr.Q(qr_hat)
    
    # Compute SVD of Q_true' * Q_hat
    # The singular values measure the principal angles between subspaces
    cross_product <- t(Q_true) %*% Q_hat
    sv <- svd(cross_product)$d
    
    # Use min_dim singular values
    min_dim <- min(r_true, r_hat)
    
    # Subspace similarity: mean of squared singular values
    # sv^2 are the squared cosines of principal angles
    similarity <- mean(sv[seq_len(min_dim)]^2)
    
    return(similarity)
  }, error = function(e) {
    warning("Error computing subspace similarity: ", e$message)
    return(NA_real_)
  })
}


#' Align Cluster Labels to Ground Truth Using Hungarian Algorithm
#'
#' Finds the optimal mapping from estimated cluster labels to true cluster
#' labels using the Hungarian algorithm. This resolves the label switching
#' problem inherent in clustering.
#'
#' @param z_true Integer vector of true cluster assignments (1 to K_true).
#' @param z_hat Integer vector of estimated cluster assignments (1 to K_hat).
#'
#' @return A list with:
#' \describe{
#'   \item{\code{mapping}}{Named integer vector where names are estimated labels
#'     and values are corresponding true labels.}
#'   \item{\code{z_aligned}}{Integer vector of aligned estimated labels
#'     (same length as z_hat, with labels remapped to match z_true).}
#'   \item{\code{accuracy}}{Numeric scalar: proportion of correctly matched
#'     assignments after alignment.}
#' }
#' Returns a list with all NA values if inputs are invalid.
#'
#' @details
#' The Hungarian algorithm (via \code{clue::solve_LSAP}) finds the assignment
#' that maximizes the overlap between estimated and true clusters. This is
#' essential for fair evaluation of clustering methods where label identities
#' are arbitrary.
#'
#' When K_hat != K_true, the algorithm handles the mismatch by:
#' \itemize{
#'   \item If K_hat < K_true: Some true clusters may not have a match
#'   \item If K_hat > K_true: Some estimated clusters map to the same true cluster
#' }
#'
#' @examples
#' \dontrun{
#' z_true <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
#' z_hat <- c(3, 3, 3, 1, 1, 1, 2, 2, 2)  # Labels are permuted
#' result <- align_clusters_to_truth(z_true, z_hat)
#' result$accuracy  # Should be 1.0 (perfect after alignment)
#' }
#'
#' @export
align_clusters_to_truth <- function(z_true, z_hat) {
  # Return structure for NA/error cases
  na_result <- list(
    mapping = NA,
    z_aligned = NA,
    accuracy = NA_real_
  )
  
  # NA-safe input validation
  if (is.null(z_true) || is.null(z_hat)) {
    return(na_result)
  }
  
  if (length(z_true) == 0 || length(z_hat) == 0) {
    return(na_result)
  }
  
  if (any(is.na(z_true)) || any(is.na(z_hat))) {
    return(na_result)
  }
  
  if (length(z_true) != length(z_hat)) {
    warning("Length mismatch: z_true has ", length(z_true),
            " elements, z_hat has ", length(z_hat), " elements")
    return(na_result)
  }
  
  # Check for clue package
  if (!requireNamespace("clue", quietly = TRUE)) {
    warning("Package 'clue' is required for cluster alignment. ",
            "Please install it with: install.packages('clue')")
    return(na_result)
  }
  
  tryCatch({
    # Get unique cluster labels
    unique_true <- sort(unique(z_true))
    unique_hat <- sort(unique(z_hat))
    
    K_true <- length(unique_true)
    K_hat <- length(unique_hat)
    
    # Build contingency matrix (K_hat x K_true)
    # Entry (i, j) = number of samples where z_hat == unique_hat[i] and z_true == unique_true[j]
    contingency <- matrix(0, nrow = K_hat, ncol = K_true)
    for (i in seq_len(K_hat)) {
      for (j in seq_len(K_true)) {
        contingency[i, j] <- sum(z_hat == unique_hat[i] & z_true == unique_true[j])
      }
    }
    
    # Handle case where K_hat != K_true by padding the matrix
    max_K <- max(K_hat, K_true)
    if (K_hat < max_K || K_true < max_K) {
      padded <- matrix(0, nrow = max_K, ncol = max_K)
      padded[seq_len(K_hat), seq_len(K_true)] <- contingency
      contingency <- padded
    }
    
    # Use Hungarian algorithm to find optimal assignment (maximize overlap)
    assignment <- clue::solve_LSAP(contingency, maximum = TRUE)
    
    # Build mapping from estimated to true labels
    mapping <- integer(K_hat)
    names(mapping) <- as.character(unique_hat)
    
    for (i in seq_len(K_hat)) {
      assigned_idx <- assignment[i]
      if (assigned_idx <= K_true) {
        mapping[i] <- unique_true[assigned_idx]
      } else {
        # Estimated cluster maps to a "phantom" true cluster (K_hat > K_true case)
        # Assign to the most overlapping true cluster
        mapping[i] <- unique_true[which.max(contingency[i, seq_len(K_true)])]
      }
    }
    
    # Apply mapping to get aligned labels
    z_aligned <- z_hat
    for (i in seq_len(K_hat)) {
      z_aligned[z_hat == unique_hat[i]] <- mapping[i]
    }
    
    # Compute accuracy after alignment
    accuracy <- mean(z_aligned == z_true)
    
    list(
      mapping = mapping,
      z_aligned = z_aligned,
      accuracy = accuracy
    )
  }, error = function(e) {
    warning("Error in cluster alignment: ", e$message)
    return(na_result)
  })
}


#' Compute Cluster-wise Subspace Recovery Score
#'
#' Computes subspace similarity scores for each cluster, using the aligned
#' mapping between estimated and true clusters. This provides a detailed
#' assessment of how well each cluster's subspace is recovered.
#'
#' @param true_Lambda_list List of K_true matrices, each (M x r) representing
#'   true cluster-specific basis matrices.
#' @param est_Lambda_list List of K_hat matrices, each (M x r_hat) representing
#'   estimated cluster-specific basis matrices.
#' @param mapping Named integer vector from \code{\link{align_clusters_to_truth}},
#'   mapping estimated cluster labels to true cluster labels.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{scores}}{Named numeric vector of subspace similarity scores
#'     for each true cluster (names are cluster indices as characters).}
#'   \item{\code{mean_score}}{Numeric scalar: mean of non-NA scores.}
#'   \item{\code{n_matched}}{Integer: number of true clusters with a valid
#'     estimated match.}
#' }
#' Returns a list with NA values if inputs are invalid.
#'
#' @details
#' For each true cluster k:
#' \enumerate{
#'   \item Find the estimated cluster that maps to k (using the mapping)
#'   \item Compute subspace similarity between true_Lambda_list[[k]] and
#'     the corresponding estimated basis
#'   \item If no estimated cluster maps to k, the score is NA
#' }
#'
#' This function is NA-safe: if estimation failed for some clusters
#' (NULL or invalid matrices), those clusters receive NA scores and
#' are excluded from the mean calculation.
#'
#' @examples
#' \dontrun{
#' # Generate data and fit model
#' data <- generate_synergy_data(N = 30, K = 3, r = 2, seed = 123)
#' fit <- mfa_em_fit(data$list_of_data, K = 3, r = 2)
#'
#' # Align clusters
#' alignment <- align_clusters_to_truth(data$true_cluster, fit$z)
#'
#' # Compute cluster-wise scores
#' scores <- compute_clusterwise_subspace_score(
#'   data$true_Lambda_list,
#'   fit$Lambda,
#'   alignment$mapping
#' )
#' }
#'
#' @export
compute_clusterwise_subspace_score <- function(true_Lambda_list, est_Lambda_list, mapping) {
  # Return structure for NA/error cases
  na_result <- list(
    scores = NA,
    mean_score = NA_real_,
    n_matched = NA_integer_
  )
  
  # NA-safe input validation
  if (is.null(true_Lambda_list) || is.null(est_Lambda_list)) {
    return(na_result)
  }
  
  if (!is.list(true_Lambda_list) || !is.list(est_Lambda_list)) {
    return(na_result)
  }
  
  if (length(true_Lambda_list) == 0 || length(est_Lambda_list) == 0) {
    return(na_result)
  }
  
  # Handle NA mapping (from failed alignment)
  if (length(mapping) == 1 && is.na(mapping)) {
    return(na_result)
  }
  
  if (is.null(mapping) || length(mapping) == 0) {
    return(na_result)
  }
  
  tryCatch({
    K_true <- length(true_Lambda_list)
    K_hat <- length(est_Lambda_list)
    
    # Get unique true cluster labels from mapping values
    unique_true <- sort(unique(as.integer(mapping)))
    
    # Initialize scores vector
    scores <- rep(NA_real_, K_true)
    names(scores) <- as.character(seq_len(K_true))
    
    # For each true cluster, find the estimated cluster that maps to it
    # and compute subspace similarity
    for (k in seq_len(K_true)) {
      # Find which estimated cluster(s) map to true cluster k
      est_indices <- which(as.integer(mapping) == k)
      
      if (length(est_indices) == 0) {
        # No estimated cluster maps to this true cluster
        scores[k] <- NA_real_
        next
      }
      
      # Use the first matching estimated cluster
      # (in case multiple estimated clusters map to the same true cluster)
      est_idx <- est_indices[1]
      
      # Get the estimated cluster label from mapping names
      est_label <- as.integer(names(mapping)[est_idx])
      
      # Validate indices
      if (est_label < 1 || est_label > K_hat) {
        scores[k] <- NA_real_
        next
      }
      
      # Get basis matrices
      L_true <- true_Lambda_list[[k]]
      L_hat <- est_Lambda_list[[est_label]]
      
      # Compute subspace similarity
      scores[k] <- compute_subspace_similarity(L_true, L_hat)
    }
    
    # Compute mean of non-NA scores
    valid_scores <- scores[!is.na(scores)]
    mean_score <- if (length(valid_scores) > 0) mean(valid_scores) else NA_real_
    n_matched <- sum(!is.na(scores))
    
    list(
      scores = scores,
      mean_score = mean_score,
      n_matched = n_matched
    )
  }, error = function(e) {
    warning("Error computing cluster-wise subspace score: ", e$message)
    return(na_result)
  })
}
