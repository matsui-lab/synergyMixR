#' Select Optimal r Using Cluster-wise Evaluation for MFA
#'
#' Selects the optimal number of factors (r) by evaluating each cluster separately,
#' rather than averaging metrics across all clusters. This approach is useful when
#' different clusters may have different optimal r values.
#'
#' @param list_of_data A list of N matrices, each (T_i x M).
#' @param K Number of clusters (fixed).
#' @param rvec Vector of candidate r values to evaluate. Default is 1:6.
#' @param criteria Character vector specifying which criteria to use for selection.
#'   Options: "VAF", "VSS", "smoothness". Default is c("VAF", "VSS", "smoothness").
#' @param agg_method Method for aggregating cluster-specific r values. Options:
#'   "max" (conservative, takes maximum r across clusters),
#'   "mean" (average r, rounded),
#'   "min" (minimum r across clusters),
#'   "majority" (most common r value).
#'   Default is "max".
#' @param vaf_threshold Minimum VAF threshold for each cluster. Set to NULL to disable
#'   hard filtering. Default is NULL (disabled).
#' @param vss_threshold Threshold for VSS criterion. Default is 0.5.
#' @param vss_method Method for computing VSS. Options: "threshold" (count variables with 
#'   exactly one strong loading above threshold), "ratio" (compute concentration ratio 
#'   max(lambda^2)/sum(lambda^2) for each variable). Default is "threshold". The "ratio" 
#'   method is more robust when loadings are small in absolute magnitude.
#' @param vss_rotate Rotation method to apply before computing VSS. Options: "none" (no rotation),
#'   "varimax" (orthogonal rotation). Default is "none". VSS on unrotated FA loadings is
#'   rotation-dependent and often uninformative. Rotation is recommended for meaningful VSS.
#' @param lambda_smooth Smoothness penalty weight (deprecated, use weights instead). Default is 1e-3.
#' @param weights Named vector of weights for combining criteria. 
#'   Default is c(VAF=0.4, VSS=0.4, smoothness=0.2).
#' @param vss_normalize Method for normalizing VSS to account for r-dependence.
#'   Options: "logr" (divide by log(r)), "max" (divide by max VSS), "r" (divide by r),
#'   "r2" (divide by r^2), "none". Default is "logr".
#' @param w_elbow Weight for DeltaVAF (elbow detection). Higher values emphasize diminishing
#'   returns in VAF. Default is 0 (disabled).
#' @param w_r Weight for direct r penalty (AIC/BIC style). Higher values penalize larger r.
#'   Default is 0 (disabled).
#' @param vaf_degeneracy_threshold Threshold for detecting degenerate clusters. If max VAF
#'   across all r values is below this threshold, the cluster is considered degenerate and
#'   a warning is issued. Default is 0.6.
#' @param max_iter Maximum EM iterations for initial clustering.
#' @param nIterFA Sub-iterations for FA update.
#' @param tol Convergence tolerance.
#' @param n_init Number of random initializations for initial clustering.
#' @param use_kmeans_init Whether to use k-means initialization.
#' @param subject_rdim_for_kmeans Dimension for k-means features.
#' @param mc_cores Number of cores for parallelization.
#' @param n_threads Number of OpenMP threads.
#' @param verbose Logical. If TRUE, print detailed diagnostic information during r selection,
#'   including per-r metrics for each cluster and warnings about NA values or degenerate clusters.
#'   Default is FALSE.
#' @param mode Selection mode. Options: "clusterwise" (default, evaluate r separately per cluster),
#'   "global-first" (first estimate global r using FA/PPCA on pooled data, then use for MFA).
#'   The "global-first" mode is recommended as it provides a cluster-independent baseline.
#' @param global_model Model for global r estimation (only used when mode = "global-first").
#'   Options: "PPCA" (default), "FA". PPCA is more stable for mixture models.
#' @param global_criterion Criterion for global r selection (only used when mode = "global-first").
#'   Options: "CV" (default), "ParallelAnalysis". CV is the primary criterion.
#' @param global_folds Number of CV folds for global r selection (only used when mode = "global-first").
#'   Default is 5.
#' @param global_seed Random seed for global r selection (only used when mode = "global-first").
#'   Default is NULL.
#' @param r_init Initial r value for the first-stage clustering (only used when mode = "clusterwise").
#'   If NULL (default), uses median(rvec). When called from hybrid mode in select_r_mixture(),
#'   this is set to r_global to use the global estimate as initialization.
#'
#' @return A list with:
#'   \item{cluster_r}{Named vector of optimal r for each cluster.}
#'   \item{agg_r}{Aggregated r value across clusters.}
#'   \item{cluster_metrics}{List of data frames with metrics for each cluster.}
#'   \item{initial_fit}{The initial MFA fit used to determine cluster assignments.}
#'   \item{agg_method}{The aggregation method used.}
#'   \item{mode}{The selection mode used.}
#'   \item{global_result}{(Only for mode = "global-first") The result from select_r_global_mfa().}
#'
#' @details
#' This function implements a two-stage approach:
#' \enumerate{
#'   \item Fit an initial MFA model with K clusters and a reasonable r value
#'         (using the median of rvec) to obtain cluster assignments.
#'   \item For each cluster k, subset the data to subjects assigned to that cluster,
#'         and evaluate different r values using cluster-specific metrics.
#'   \item Select the optimal r_k for each cluster based on normalized scoring that
#'         combines VAF (higher is better), VSS (higher is better, normalized by log(r)),
#'         and smoothness (lower is better) with specified weights.
#'   \item Aggregate the cluster-specific r values using the specified method.
#' }
#'
#' The scoring approach addresses the issue where VAF monotonically increases with r
#' by normalizing all metrics per cluster and using weighted combination. VSS is
#' normalized by log(r) to account for the fact that more factors naturally lead to
#' more complex structures. Smoothness acts as a penalty (subtracted from score).
#'
#' The cluster-wise approach is particularly useful when:
#' \itemize{
#'   \item Different motor control strategies require different numbers of synergies
#'   \item Averaging metrics across clusters masks cluster-specific structure
#'   \item You want to understand the heterogeneity in synergy complexity
#' }
#'
#' \strong{Important Note:} The cluster-wise MFA(K=1) approach is experimental and can be
#' unstable. Each cluster is evaluated independently using a single-cluster MFA model, which
#' may not capture the full complexity of the data. Use with caution and carefully inspect
#' the verbose output and degeneracy warnings.
#'
#' @examples
#' \dontrun{
#' # Basic usage with default weights (clusterwise mode)
#' result <- select_r_clusterwise_mfa(
#'   list_of_data, K = 3, rvec = 2:6,
#'   criteria = c("VAF", "VSS", "smoothness")
#' )
#' print(result$cluster_r)  # r for each cluster
#' print(result$agg_r)      # aggregated r
#' 
#' # Global-first mode (recommended)
#' result <- select_r_clusterwise_mfa(
#'   list_of_data, K = 3, rvec = 2:6,
#'   mode = "global-first",
#'   global_model = "PPCA",
#'   global_criterion = "CV"
#' )
#' print(result$agg_r)  # global r estimate
#' 
#' # Custom weights emphasizing VAF (clusterwise mode)
#' result <- select_r_clusterwise_mfa(
#'   list_of_data, K = 3, rvec = 2:6,
#'   mode = "clusterwise",
#'   criteria = c("VAF", "VSS", "smoothness"),
#'   weights = c(VAF = 0.6, VSS = 0.3, smoothness = 0.1)
#' )
#' }
#'
#' @export
select_r_clusterwise_mfa <- function(
    list_of_data,
    K,
    rvec = 1:6,
    criteria = c("VAF", "VSS", "smoothness"),
    agg_method = c("max", "mean", "min", "majority"),
    vaf_threshold = NULL,
    vss_threshold = 0.5,
    vss_method = c("threshold", "ratio"),
    vss_rotate = c("none", "varimax"),
    lambda_smooth = 1e-3,
    weights = c(VAF = 0.4, VSS = 0.4, smoothness = 0.2),
    vss_normalize = c("logr", "max", "r", "r2", "none"),
    w_elbow = 0,
    w_r = 0,
    vaf_degeneracy_threshold = 0.6,
    max_iter = 50,
    nIterFA = 5,
    tol = 1e-3,
    n_init = 1,
    use_kmeans_init = TRUE,
    subject_rdim_for_kmeans = 2,
    mc_cores = 1,
    n_threads = 1,
    verbose = FALSE,
    mode = c("clusterwise", "global-first"),
    global_model = c("PPCA", "FA"),
    global_criterion = c("CV", "ParallelAnalysis"),
    global_folds = 5,
    global_seed = NULL,
    r_init = NULL
) {
  # Input validation
  mode <- match.arg(mode)
  agg_method <- match.arg(agg_method)
  vss_method <- match.arg(vss_method)
  vss_rotate <- match.arg(vss_rotate)
  criteria <- match.arg(criteria, several.ok = TRUE)
  global_model <- match.arg(global_model)
  global_criterion <- match.arg(global_criterion)
  
  if (length(rvec) < 2) {
    stop("rvec must contain at least 2 candidate r values")
  }
  
  N <- length(list_of_data)
  M <- ncol(list_of_data[[1]])
  
  cat("=== Cluster-wise r selection for MFA ===\n")
  cat(sprintf("Mode: %s\n", mode))
  cat(sprintf("K = %d, rvec = %s\n", K, paste(rvec, collapse = ", ")))
  
  # Handle global-first mode
  if (mode == "global-first") {
    cat(sprintf("\nUsing global-first mode: estimating r using %s with %s\n", 
                global_model, global_criterion))
    
    global_result <- select_r_global_mfa(
      list_of_data = list_of_data,
      rvec = rvec,
      model = global_model,
      criterion = global_criterion,
      folds = global_folds,
      seed = global_seed,
      verbose = verbose
    )
    
    agg_r <- global_result$r_global
    
    cat(sprintf("\nGlobal r estimate: %d\n", agg_r))
    cat("Returning global r as the final estimate.\n")
    
    # Return early with global result
    return(list(
      cluster_r = NULL,
      agg_r = agg_r,
      cluster_metrics = NULL,
      initial_fit = NULL,
      agg_method = "global",
      mode = mode,
      global_result = global_result
    ))
  }
  
  # Continue with clusterwise mode
  cat(sprintf("Criteria: %s\n", paste(criteria, collapse = ", ")))
  cat(sprintf("Aggregation method: %s\n", agg_method))
  
  # Step 1: Fit initial MFA model to get cluster assignments
  # Use provided r_init or default to median r value for initial clustering
  if (is.null(r_init)) {
    r_init <- round(median(rvec))
  }
  cat(sprintf("\nStep 1: Fitting initial MFA with r = %d to determine cluster assignments...\n", r_init))
  
  initial_fit <- mfa_em_fit(
    list_of_data = list_of_data,
    K = K,
    r = r_init,
    max_iter = max_iter,
    nIterFA = nIterFA,
    tol = tol,
    n_init = n_init,
    use_kmeans_init = use_kmeans_init,
    subject_rdim_for_kmeans = subject_rdim_for_kmeans,
    mc_cores = mc_cores,
    n_threads = n_threads
  )
  
  z_vec <- initial_fit$z
  cluster_sizes <- table(z_vec)
  cat("Cluster sizes: [", paste(cluster_sizes, collapse = ", "), "]\n", sep = "")
  
  # Step 2: For each cluster, evaluate different r values
  cat("\nStep 2: Evaluating r values for each cluster...\n")
  
  cluster_r <- numeric(K)
  names(cluster_r) <- paste0("Cluster_", 1:K)
  cluster_metrics_list <- vector("list", K)
  
  for (k in 1:K) {
    cat(sprintf("\n--- Cluster %d (n = %d subjects) ---\n", k, sum(z_vec == k)))
    
    # Subset data for this cluster
    cluster_indices <- which(z_vec == k)
    
    if (length(cluster_indices) < 2) {
      warning(sprintf("Cluster %d has fewer than 2 subjects. Using minimum r.", k))
      cluster_r[k] <- min(rvec)
      cluster_metrics_list[[k]] <- data.frame(
        r = min(rvec),
        VAF = NA,
        VSS = NA,
        smoothness = NA,
        note = "Too few subjects"
      )
      next
    }
    
    cluster_data <- list_of_data[cluster_indices]
    
    # Evaluate each r value for this cluster
    cluster_results <- list()
    
    for (r in rvec) {
      # Fit single-cluster MFA (K=1) with this r
      fit_single <- mfa_em_fit(
        list_of_data = cluster_data,
        K = 1,
        r = r,
        max_iter = max_iter,
        nIterFA = nIterFA,
        tol = tol,
        n_init = 1,
        use_kmeans_init = FALSE,
        mc_cores = 1,
        n_threads = n_threads
      )
      
      # Compute cluster-specific metrics
      vaf_val <- NA
      vss_val <- NA
      smoothness_val <- NA
      
      if ("VAF" %in% criteria) {
        vaf_val <- compute_global_vaf_mfa(cluster_data, fit_single)
      }
      
      if ("VSS" %in% criteria) {
        vss_val <- compute_vss_criterion(fit_single$Lambda, threshold = vss_threshold, method = vss_method, rotate = vss_rotate)
      }
      
      if ("smoothness" %in% criteria) {
        tryCatch({
          fit_single <- compute_factor_scores_mfa(fit_single, cluster_data)
          if (!is.null(fit_single$factor_scores)) {
            smoothness_val <- compute_smoothness_penalty(fit_single$factor_scores)
            # Check for NA smoothness and provide diagnostic
            if (is.na(smoothness_val) && verbose) {
              n_na <- sum(sapply(fit_single$factor_scores, function(x) any(is.na(x))))
              cat(sprintf("    WARNING: Smoothness is NA for r=%d (factor scores contain %d NA subjects)\n", r, n_na))
            }
          } else {
            smoothness_val <- NA
            if (verbose) {
              cat(sprintf("    WARNING: Factor scores are NULL for r=%d\n", r))
            }
          }
        }, error = function(e) {
          smoothness_val <<- NA
          if (verbose) {
            cat(sprintf("    WARNING: Smoothness computation failed for r=%d: %s\n", r, e$message))
          }
        })
      }
      
      cluster_results[[length(cluster_results) + 1]] <- list(
        r = r,
        VAF = vaf_val,
        VSS = vss_val,
        smoothness = smoothness_val
      )
    }
    
    # Build summary for this cluster
    df_cluster <- do.call(rbind, lapply(cluster_results, function(res) {
      data.frame(
        r = res$r,
        VAF = res$VAF,
        VSS = res$VSS,
        smoothness = res$smoothness
      )
    }))
    
    # Check for cluster degeneracy: if max VAF across r is too low, the cluster is degenerate
    if ("VAF" %in% criteria && !all(is.na(df_cluster$VAF))) {
      max_vaf <- max(df_cluster$VAF, na.rm = TRUE)
      if (max_vaf < vaf_degeneracy_threshold) {
        warning(sprintf("Cluster %d is degenerate: max VAF = %.3f < threshold = %.3f. R-selection may be unreliable.", 
                        k, max_vaf, vaf_degeneracy_threshold))
        if (verbose) {
          cat(sprintf("  WARNING: Cluster %d has low VAF (max = %.3f), r-selection unreliable\n", k, max_vaf))
        }
      }
    }
    
    # Verbose logging: print per-r summaries for this cluster
    if (verbose) {
      cat(sprintf("\n  Cluster %d metrics:\n", k))
      for (i in seq_len(nrow(df_cluster))) {
        cat(sprintf("    r=%d: VAF=%.3f, VSS=%.3f, smoothness=%.3f\n",
                    df_cluster$r[i],
                    ifelse(is.na(df_cluster$VAF[i]), NA, df_cluster$VAF[i]),
                    ifelse(is.na(df_cluster$VSS[i]), NA, df_cluster$VSS[i]),
                    ifelse(is.na(df_cluster$smoothness[i]), NA, df_cluster$smoothness[i])))
      }
    }
    
    # Select best r for this cluster based on criteria
    # Use normalized scoring with proper weighting
    best_r_k <- select_best_r_from_metrics(
      df_cluster, 
      criteria, 
      vaf_threshold, 
      lambda_smooth,
      weights = weights,
      vss_normalize = vss_normalize[1],
      w_elbow = w_elbow,
      w_r = w_r,
      verbose = verbose
    )
    
    cluster_r[k] <- best_r_k
    cluster_metrics_list[[k]] <- df_cluster
    
    cat(sprintf("Selected r = %d for Cluster %d\n", best_r_k, k))
  }
  
  # Step 3: Aggregate cluster-specific r values
  cat("\n--- Aggregating cluster-specific r values ---\n")
  cat("Cluster-specific r values:", paste(cluster_r, collapse = ", "), "\n")
  
  agg_r <- aggregate_r_values(cluster_r, agg_method)
  
  cat(sprintf("Aggregated r (method = %s): %d\n", agg_method, agg_r))
  
  list(
    cluster_r = cluster_r,
    agg_r = agg_r,
    cluster_metrics = cluster_metrics_list,
    initial_fit = initial_fit,
    agg_method = agg_method,
    mode = mode,
    global_result = NULL
  )
}


#' Select Optimal r Using Cluster-wise Evaluation for MPCA
#'
#' Selects the optimal number of principal components (r) by evaluating each cluster
#' separately, rather than averaging metrics across all clusters.
#'
#' @param list_of_data A list of N matrices, each (T_i x M).
#' @param K Number of clusters (fixed).
#' @param rvec Vector of candidate r values to evaluate. Default is 1:6.
#' @param criteria Character vector specifying which criteria to use for selection.
#'   Options: "VAF", "VSS", "smoothness". Default is c("VAF", "VSS", "smoothness").
#' @param agg_method Method for aggregating cluster-specific r values. Options:
#'   "max", "mean", "min", "majority". Default is "max".
#' @param vaf_threshold Minimum VAF threshold for each cluster. Set to NULL to disable
#'   hard filtering. Default is NULL (disabled).
#' @param vss_threshold Threshold for VSS criterion. Default is 0.5.
#' @param vss_method Method for computing VSS. Options: "threshold" (count variables with 
#'   exactly one strong loading above threshold), "ratio" (compute concentration ratio 
#'   max(lambda^2)/sum(lambda^2) for each variable). Default is "threshold". The "ratio" 
#'   method is more robust when loadings are small in absolute magnitude.
#' @param vss_rotate Rotation method to apply before computing VSS. Options: "none" (no rotation),
#'   "varimax" (orthogonal rotation). Default is "none". VSS on unrotated FA loadings is
#'   rotation-dependent and often uninformative. Rotation is recommended for meaningful VSS.
#' @param lambda_smooth Smoothness penalty weight (deprecated, use weights instead). Default is 1e-3.
#' @param weights Named vector of weights for combining criteria. 
#'   Default is c(VAF=0.4, VSS=0.4, smoothness=0.2).
#' @param vss_normalize Method for normalizing VSS to account for r-dependence.
#'   Options: "logr" (divide by log(r)), "max" (divide by max VSS), "r" (divide by r),
#'   "r2" (divide by r^2), "none". Default is "logr".
#' @param w_elbow Weight for DeltaVAF (elbow detection). Higher values emphasize diminishing
#'   returns in VAF. Default is 0 (disabled).
#' @param w_r Weight for direct r penalty (AIC/BIC style). Higher values penalize larger r.
#'   Default is 0 (disabled).
#' @param vaf_degeneracy_threshold Threshold for detecting degenerate clusters. If max VAF
#'   across all r values is below this threshold, the cluster is considered degenerate and
#'   a warning is issued. Default is 0.6.
#' @param max_iter Maximum EM iterations for initial clustering.
#' @param nIterPCA Sub-iterations for PCA update.
#' @param tol Convergence tolerance.
#' @param method "EM" or "closed_form".
#' @param n_init Number of random initializations for initial clustering.
#' @param use_kmeans_init Whether to use k-means initialization.
#' @param subject_rdim_for_kmeans Dimension for k-means features.
#' @param mc_cores Number of cores for parallelization.
#' @param n_threads Number of OpenMP threads.
#' @param verbose Logical. If TRUE, print detailed diagnostic information during r selection,
#'   including per-r metrics for each cluster and warnings about NA values or degenerate clusters.
#'   Default is FALSE.
#' @param mode Selection mode. Options: "clusterwise" (default, evaluate r separately per cluster),
#'   "global-first" (first estimate global r using FA/PPCA on pooled data, then use for MPCA).
#'   The "global-first" mode is recommended as it provides a cluster-independent baseline.
#' @param global_model Model for global r estimation (only used when mode = "global-first").
#'   Options: "PPCA" (default), "FA". PPCA is more stable for mixture models.
#' @param global_criterion Criterion for global r selection (only used when mode = "global-first").
#'   Options: "CV" (default), "ParallelAnalysis". CV is the primary criterion.
#' @param global_folds Number of CV folds for global r selection (only used when mode = "global-first").
#'   Default is 5.
#' @param global_seed Random seed for global r selection (only used when mode = "global-first").
#'   Default is NULL.
#' @param r_init Initial r value for the first-stage clustering (only used when mode = "clusterwise").
#'   If NULL (default), uses median(rvec). When called from hybrid mode in select_r_mixture(),
#'   this is set to r_global to use the global estimate as initialization.
#'
#' @return A list with:
#'   \item{cluster_r}{Named vector of optimal r for each cluster.}
#'   \item{agg_r}{Aggregated r value across clusters.}
#'   \item{cluster_metrics}{List of data frames with metrics for each cluster.}
#'   \item{initial_fit}{The initial MPCA fit used to determine cluster assignments.}
#'   \item{agg_method}{The aggregation method used.}
#'   \item{mode}{The selection mode used.}
#'   \item{global_result}{(Only for mode = "global-first") The result from select_r_global_mfa().}
#'
#' @details
#' This function implements the same two-stage approach as
#' \code{\link{select_r_clusterwise_mfa}}, but for Mixture PCA models.
#' See \code{\link{select_r_clusterwise_mfa}} for details on the scoring approach.
#'
#' \strong{Important Note:} The cluster-wise MPCA(K=1) approach is experimental and can be
#' unstable. Each cluster is evaluated independently using a single-cluster MPCA model, which
#' may not capture the full complexity of the data. Use with caution and carefully inspect
#' the verbose output and degeneracy warnings.
#'
#' @examples
#' \dontrun{
#' result <- select_r_clusterwise_mpca(
#'   list_of_data, K = 3, rvec = 2:6,
#'   criteria = c("VAF", "VSS", "smoothness")
#' )
#' print(result$cluster_r)
#' print(result$agg_r)
#' }
#'
#' @export
select_r_clusterwise_mpca <- function(
    list_of_data,
    K,
    rvec = 1:6,
    criteria = c("VAF", "VSS", "smoothness"),
    agg_method = c("max", "mean", "min", "majority"),
    vaf_threshold = NULL,
    vss_threshold = 0.5,
    vss_method = c("threshold", "ratio"),
    vss_rotate = c("none", "varimax"),
    lambda_smooth = 1e-3,
    weights = c(VAF = 0.4, VSS = 0.4, smoothness = 0.2),
    vss_normalize = c("logr", "max", "r", "r2", "none"),
    w_elbow = 0,
    w_r = 0,
    vaf_degeneracy_threshold = 0.6,
    max_iter = 50,
    nIterPCA = 5,
    tol = 1e-3,
    method = "EM",
    n_init = 1,
    use_kmeans_init = TRUE,
    subject_rdim_for_kmeans = 2,
    mc_cores = 1,
    n_threads = 1,
    verbose = FALSE,
    mode = c("clusterwise", "global-first"),
    global_model = c("PPCA", "FA"),
    global_criterion = c("CV", "ParallelAnalysis"),
    global_folds = 5,
    global_seed = NULL,
    r_init = NULL
) {
  # Input validation
  mode <- match.arg(mode)
  agg_method <- match.arg(agg_method)
  vss_method <- match.arg(vss_method)
  vss_rotate <- match.arg(vss_rotate)
  criteria <- match.arg(criteria, several.ok = TRUE)
  global_model <- match.arg(global_model)
  global_criterion <- match.arg(global_criterion)
  
  if (length(rvec) < 2) {
    stop("rvec must contain at least 2 candidate r values")
  }
  
  N <- length(list_of_data)
  M <- ncol(list_of_data[[1]])
  
  cat("=== Cluster-wise r selection for MPCA ===\n")
  cat(sprintf("Mode: %s\n", mode))
  cat(sprintf("K = %d, rvec = %s\n", K, paste(rvec, collapse = ", ")))
  
  # Handle global-first mode
  if (mode == "global-first") {
    cat(sprintf("\nUsing global-first mode: estimating r using %s with %s\n", 
                global_model, global_criterion))
    
    global_result <- select_r_global_mfa(
      list_of_data = list_of_data,
      rvec = rvec,
      model = global_model,
      criterion = global_criterion,
      folds = global_folds,
      seed = global_seed,
      verbose = verbose
    )
    
    agg_r <- global_result$r_global
    
    cat(sprintf("\nGlobal r estimate: %d\n", agg_r))
    cat("Returning global r as the final estimate.\n")
    
    # Return early with global result
    return(list(
      cluster_r = NULL,
      agg_r = agg_r,
      cluster_metrics = NULL,
      initial_fit = NULL,
      agg_method = "global",
      mode = mode,
      global_result = global_result
    ))
  }
  
  # Continue with clusterwise mode
  cat(sprintf("Criteria: %s\n", paste(criteria, collapse = ", ")))
  cat(sprintf("Aggregation method: %s\n", agg_method))
  
  # Step 1: Fit initial MPCA model to get cluster assignments
  # Use provided r_init or default to median r value for initial clustering
  if (is.null(r_init)) {
    r_init <- round(median(rvec))
  }
  cat(sprintf("\nStep 1: Fitting initial MPCA with r = %d to determine cluster assignments...\n", r_init))
  
  initial_fit <- mixture_pca_em_fit(
    list_of_data = list_of_data,
    K = K,
    r = r_init,
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
  
  # Ensure W field exists
  initial_fit <- ensure_mpca_has_W(initial_fit)
  
  z_vec <- initial_fit$z
  cluster_sizes <- table(z_vec)
  cat("Cluster sizes: [", paste(cluster_sizes, collapse = ", "), "]\n", sep = "")
  
  # Step 2: For each cluster, evaluate different r values
  cat("\nStep 2: Evaluating r values for each cluster...\n")
  
  cluster_r <- numeric(K)
  names(cluster_r) <- paste0("Cluster_", 1:K)
  cluster_metrics_list <- vector("list", K)
  
  for (k in 1:K) {
    cat(sprintf("\n--- Cluster %d (n = %d subjects) ---\n", k, sum(z_vec == k)))
    
    # Subset data for this cluster
    cluster_indices <- which(z_vec == k)
    
    if (length(cluster_indices) < 2) {
      warning(sprintf("Cluster %d has fewer than 2 subjects. Using minimum r.", k))
      cluster_r[k] <- min(rvec)
      cluster_metrics_list[[k]] <- data.frame(
        r = min(rvec),
        VAF = NA,
        VSS = NA,
        smoothness = NA,
        note = "Too few subjects"
      )
      next
    }
    
    cluster_data <- list_of_data[cluster_indices]
    
    # Evaluate each r value for this cluster
    cluster_results <- list()
    
    for (r in rvec) {
      # Fit single-cluster MPCA (K=1) with this r
      fit_single <- mixture_pca_em_fit(
        list_of_data = cluster_data,
        K = 1,
        r = r,
        max_iter = max_iter,
        nIterPCA = nIterPCA,
        tol = tol,
        method = method,
        n_init = 1,
        use_kmeans_init = FALSE,
        mc_cores = 1,
        n_threads = n_threads
      )
      
      # Ensure W field exists
      fit_single <- ensure_mpca_has_W(fit_single)
      
      # Compute cluster-specific metrics
      vaf_val <- NA
      vss_val <- NA
      smoothness_val <- NA
      
      if ("VAF" %in% criteria) {
        vaf_val <- compute_global_vaf_mpca(cluster_data, fit_single)
      }
      
      if ("VSS" %in% criteria) {
        vss_val <- compute_vss_criterion(fit_single$W, threshold = vss_threshold, method = vss_method, rotate = vss_rotate)
      }
      
      if ("smoothness" %in% criteria) {
        tryCatch({
          fit_single <- compute_factor_scores_mpca(fit_single, cluster_data)
          if (!is.null(fit_single$factor_scores)) {
            smoothness_val <- compute_smoothness_penalty(fit_single$factor_scores)
            # Check for NA smoothness and provide diagnostic
            if (is.na(smoothness_val) && verbose) {
              n_na <- sum(sapply(fit_single$factor_scores, function(x) any(is.na(x))))
              cat(sprintf("    WARNING: Smoothness is NA for r=%d (factor scores contain %d NA subjects)\n", r, n_na))
            }
          } else {
            smoothness_val <- NA
            if (verbose) {
              cat(sprintf("    WARNING: Factor scores are NULL for r=%d\n", r))
            }
          }
        }, error = function(e) {
          smoothness_val <<- NA
          if (verbose) {
            cat(sprintf("    WARNING: Smoothness computation failed for r=%d: %s\n", r, e$message))
          }
        })
      }
      
      cluster_results[[length(cluster_results) + 1]] <- list(
        r = r,
        VAF = vaf_val,
        VSS = vss_val,
        smoothness = smoothness_val
      )
    }
    
    # Build summary for this cluster
    df_cluster <- do.call(rbind, lapply(cluster_results, function(res) {
      data.frame(
        r = res$r,
        VAF = res$VAF,
        VSS = res$VSS,
        smoothness = res$smoothness
      )
    }))
    
    # Check for cluster degeneracy: if max VAF across r is too low, the cluster is degenerate
    if ("VAF" %in% criteria && !all(is.na(df_cluster$VAF))) {
      max_vaf <- max(df_cluster$VAF, na.rm = TRUE)
      if (max_vaf < vaf_degeneracy_threshold) {
        warning(sprintf("Cluster %d is degenerate: max VAF = %.3f < threshold = %.3f. R-selection may be unreliable.", 
                        k, max_vaf, vaf_degeneracy_threshold))
        if (verbose) {
          cat(sprintf("  WARNING: Cluster %d has low VAF (max = %.3f), r-selection unreliable\n", k, max_vaf))
        }
      }
    }
    
    # Verbose logging: print per-r summaries for this cluster
    if (verbose) {
      cat(sprintf("\n  Cluster %d metrics:\n", k))
      for (i in seq_len(nrow(df_cluster))) {
        cat(sprintf("    r=%d: VAF=%.3f, VSS=%.3f, smoothness=%.3f\n",
                    df_cluster$r[i],
                    ifelse(is.na(df_cluster$VAF[i]), NA, df_cluster$VAF[i]),
                    ifelse(is.na(df_cluster$VSS[i]), NA, df_cluster$VSS[i]),
                    ifelse(is.na(df_cluster$smoothness[i]), NA, df_cluster$smoothness[i])))
      }
    }
    
    # Select best r for this cluster
    best_r_k <- select_best_r_from_metrics(
      df_cluster, 
      criteria, 
      vaf_threshold, 
      lambda_smooth,
      weights = weights,
      vss_normalize = vss_normalize[1],
      w_elbow = w_elbow,
      w_r = w_r,
      verbose = verbose
    )
    
    cluster_r[k] <- best_r_k
    cluster_metrics_list[[k]] <- df_cluster
    
    cat(sprintf("Selected r = %d for Cluster %d\n", best_r_k, k))
  }
  
  # Step 3: Aggregate cluster-specific r values
  cat("\n--- Aggregating cluster-specific r values ---\n")
  cat("Cluster-specific r values:", paste(cluster_r, collapse = ", "), "\n")
  
  agg_r <- aggregate_r_values(cluster_r, agg_method)
  
  cat(sprintf("Aggregated r (method = %s): %d\n", agg_method, agg_r))
  
  list(
    cluster_r = cluster_r,
    agg_r = agg_r,
    cluster_metrics = cluster_metrics_list,
    initial_fit = initial_fit,
    agg_method = agg_method,
    mode = mode,
    global_result = NULL
  )
}


#' Select Best r from Metrics Data Frame
#'
#' Helper function to select the best r value from a data frame of metrics
#' using normalized scoring with proper weighting.
#'
#' @param df_metrics Data frame with columns r, VAF, VSS, smoothness.
#' @param criteria Character vector of criteria to use.
#' @param vaf_threshold Minimum VAF threshold (optional, NULL to disable).
#' @param lambda_smooth Smoothness penalty weight (deprecated, use weights instead).
#' @param weights Named vector of weights for each criterion. 
#'   Default: c(VAF=0.4, VSS=0.4, smoothness=0.2).
#' @param vss_normalize Method for normalizing VSS: "logr" (divide by log(r)), 
#'   "max" (divide by max), "r" (divide by r), "r2" (divide by r^2), or "none".
#'   Default: "logr".
#' @param w_elbow Weight for DeltaVAF (elbow detection). Default: 0 (disabled).
#' @param w_r Weight for direct r penalty. Default: 0 (disabled).
#' @param verbose Logical. If TRUE, print diagnostic information. Default: FALSE.
#'
#' @return The selected r value.
#' @keywords internal
select_best_r_from_metrics <- function(df_metrics, criteria, vaf_threshold = NULL, 
                                       lambda_smooth = 1e-3,
                                       weights = c(VAF = 0.4, VSS = 0.4, smoothness = 0.2),
                                       vss_normalize = "logr",
                                       w_elbow = 0,
                                       w_r = 0,
                                       verbose = FALSE) {
  
  # Min-max normalization with safe fallback for constant values
  normalize_minmax <- function(x) {
    x_clean <- x[!is.na(x)]
    if (length(x_clean) == 0) return(rep(NA, length(x)))
    if (length(unique(x_clean)) == 1) return(rep(0.5, length(x)))
    x_min <- min(x_clean, na.rm = TRUE)
    x_max <- max(x_clean, na.rm = TRUE)
    (x - x_min) / (x_max - x_min)
  }
  
  # Apply optional VAF threshold as soft filter (warn but don't hard-filter)
  if (!is.null(vaf_threshold) && "VAF" %in% criteria && !all(is.na(df_metrics$VAF))) {
    n_above_threshold <- sum(df_metrics$VAF >= vaf_threshold, na.rm = TRUE)
    if (n_above_threshold == 0) {
      warning(sprintf("No r values satisfy VAF threshold (%.3f). Proceeding with all r values.", 
                      vaf_threshold))
    }
  }
  
  # If only one criterion, use it directly (simple case)
  if (length(criteria) == 1) {
    if (criteria == "VAF") {
      best_idx <- which.max(df_metrics$VAF)
      return(df_metrics$r[best_idx])
    } else if (criteria == "VSS") {
      best_idx <- which.max(df_metrics$VSS)
      return(df_metrics$r[best_idx])
    } else if (criteria == "smoothness") {
      best_idx <- which.min(df_metrics$smoothness)
      return(df_metrics$r[best_idx])
    }
  }
  
  # Multiple criteria: normalize and combine with weighted scoring
  # Initialize score vector
  score <- rep(0, nrow(df_metrics))
  n_criteria <- 0
  
  # Degeneracy detector: check if non-VAF metrics are flat/NA
  # If so, automatically enable fallbacks even when w_elbow=0 and w_r=0
  vss_degenerate <- FALSE
  smooth_degenerate <- FALSE
  
  if ("VSS" %in% criteria && !all(is.na(df_metrics$VSS))) {
    vss_clean <- df_metrics$VSS[!is.na(df_metrics$VSS)]
    if (length(vss_clean) < 2 || (max(vss_clean) - min(vss_clean)) < 1e-10) {
      vss_degenerate <- TRUE
    }
  } else if ("VSS" %in% criteria) {
    vss_degenerate <- TRUE  # All NA
  }
  
  if ("smoothness" %in% criteria && !all(is.na(df_metrics$smoothness))) {
    smooth_clean <- df_metrics$smoothness[!is.na(df_metrics$smoothness)]
    if (length(smooth_clean) < 2 || (max(smooth_clean) - min(smooth_clean)) < 1e-10) {
      smooth_degenerate <- TRUE
    }
  } else if ("smoothness" %in% criteria) {
    smooth_degenerate <- TRUE  # All NA
  }
  
  # Check if all non-VAF metrics are degenerate
  non_vaf_criteria <- setdiff(criteria, "VAF")
  all_degenerate <- FALSE
  if (length(non_vaf_criteria) > 0) {
    degenerate_count <- sum(c(
      if ("VSS" %in% non_vaf_criteria) vss_degenerate else FALSE,
      if ("smoothness" %in% non_vaf_criteria) smooth_degenerate else FALSE
    ))
    all_degenerate <- (degenerate_count == length(non_vaf_criteria))
  }
  
  # If all non-VAF metrics are degenerate, enable fallbacks automatically
  w_elbow_effective <- w_elbow
  w_r_effective <- w_r
  if (all_degenerate && "VAF" %in% criteria) {
    if (w_elbow == 0) {
      w_elbow_effective <- 0.5  # Auto-enable elbow detection
      message("Degeneracy detected: all non-VAF metrics are flat/NA. Auto-enabling DeltaVAF fallback (w_elbow=0.5).")
    }
    if (w_r == 0) {
      w_r_effective <- 0.2  # Auto-enable r penalty
      message("Degeneracy detected: all non-VAF metrics are flat/NA. Auto-enabling r penalty fallback (w_r=0.2).")
    }
  }
  
  # Process VAF (higher is better)
  if ("VAF" %in% criteria && !all(is.na(df_metrics$VAF))) {
    vaf_norm <- normalize_minmax(df_metrics$VAF)
    w <- if ("VAF" %in% names(weights)) weights["VAF"] else 0.4
    score <- score + w * vaf_norm
    n_criteria <- n_criteria + 1
  }
  
  # Process VSS (higher is better, but normalize by r first)
  if ("VSS" %in% criteria && !all(is.na(df_metrics$VSS))) {
    # Adjust VSS to account for r-dependence
    if (vss_normalize == "logr") {
      # Normalize by log(r) to account for increasing complexity
      vss_adj <- df_metrics$VSS / log(pmax(df_metrics$r, 2))
    } else if (vss_normalize == "max") {
      # Normalize by maximum VSS value
      max_vss <- max(df_metrics$VSS, na.rm = TRUE)
      vss_adj <- if (max_vss > 0) df_metrics$VSS / max_vss else df_metrics$VSS
    } else if (vss_normalize == "r") {
      # Normalize by r (stronger penalty for larger r)
      vss_adj <- df_metrics$VSS / pmax(df_metrics$r, 1)
    } else if (vss_normalize == "r2") {
      # Normalize by r^2 (even stronger penalty for larger r)
      vss_adj <- df_metrics$VSS / pmax(df_metrics$r^2, 1)
    } else if (vss_normalize == "none") {
      # No normalization
      vss_adj <- df_metrics$VSS
    } else {
      vss_adj <- df_metrics$VSS
    }
    
    vss_norm <- normalize_minmax(vss_adj)
    w <- if ("VSS" %in% names(weights)) weights["VSS"] else 0.4
    score <- score + w * vss_norm
    n_criteria <- n_criteria + 1
  }
  
  # Process smoothness (lower is better, so subtract)
  if ("smoothness" %in% criteria && !all(is.na(df_metrics$smoothness))) {
    smooth_norm <- normalize_minmax(df_metrics$smoothness)
    w <- if ("smoothness" %in% names(weights)) weights["smoothness"] else 0.2
    score <- score - w * smooth_norm  # Subtract because lower smoothness is better
    n_criteria <- n_criteria + 1
  }
  
  # Add DeltaVAF (elbow detection) if w_elbow_effective > 0
  if (w_elbow_effective > 0 && "VAF" %in% names(df_metrics) && !all(is.na(df_metrics$VAF))) {
    # Sort by r to ensure correct ordering
    df_sorted <- df_metrics[order(df_metrics$r), ]
    
    # Compute DeltaVAF = VAF(r) - VAF(r-1)
    delta_vaf <- c(0, diff(df_sorted$VAF))  # First r gets 0
    
    # Normalize DeltaVAF (higher is better - indicates elbow)
    delta_vaf_norm <- normalize_minmax(delta_vaf)
    
    # Reorder to match original df_metrics order
    delta_vaf_norm_reordered <- delta_vaf_norm[match(df_metrics$r, df_sorted$r)]
    
    # Add to score
    score <- score + w_elbow_effective * delta_vaf_norm_reordered
  }
  
  # Add direct r penalty if w_r_effective > 0
  if (w_r_effective > 0) {
    # Normalize r values (0 to 1)
    r_norm <- normalize_minmax(df_metrics$r)
    
    # Subtract penalty (higher r gets penalized more)
    score <- score - w_r_effective * r_norm
  }
  
  # Handle case where all criteria are NA
  if (n_criteria == 0 && w_elbow_effective == 0 && w_r_effective == 0) {
    warning("All criteria resulted in NA values. Selecting minimum r.")
    return(min(df_metrics$r))
  }
  
  # Select r with maximum score (higher score is better)
  # In case of ties, choose the smallest r (parsimony principle)
  max_score <- max(score, na.rm = TRUE)
  tied_indices <- which(abs(score - max_score) < 1e-10)
  
  if (length(tied_indices) > 1) {
    # Multiple r values have the same score, choose smallest
    best_idx <- tied_indices[which.min(df_metrics$r[tied_indices])]
  } else {
    best_idx <- tied_indices[1]
  }
  
  return(df_metrics$r[best_idx])
}


#' Aggregate Cluster-specific r Values
#'
#' Helper function to aggregate r values across clusters using different methods.
#'
#' @param cluster_r Named numeric vector of r values for each cluster.
#' @param method Aggregation method: "max", "mean", "min", or "majority".
#'
#' @return The aggregated r value (integer).
#' @keywords internal
aggregate_r_values <- function(cluster_r, method) {
  if (method == "max") {
    return(as.integer(max(cluster_r)))
  } else if (method == "mean") {
    return(as.integer(round(mean(cluster_r))))
  } else if (method == "min") {
    return(as.integer(min(cluster_r)))
  } else if (method == "majority") {
    # Most common r value (mode)
    freq_table <- table(cluster_r)
    mode_r <- as.integer(names(freq_table)[which.max(freq_table)])
    return(mode_r)
  } else {
    stop("Unknown aggregation method: ", method)
  }
}
