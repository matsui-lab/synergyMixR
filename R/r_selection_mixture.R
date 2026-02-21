# ============================================================
# Unified r Selection Interface (Phase 2 API Unification)
# ============================================================

#' Select Optimal Number of Synergies (r)
#'
#' Main entry point for selecting the optimal number of synergies (latent factors)
#' in mixture models. This unified interface dispatches to specialized selection
#' methods based on the specified parameters.
#'
#' @param list_of_data A list of N matrices, each of dimension (T_i x M)
#' @param model_type Character, either "MFA" or "MPCA"
#' @param method Character, selection method. One of:
#'   \itemize{
#'     \item "vaf" - Variance Accounted For threshold
#'     \item "vss" - Very Simple Structure criterion
#'     \item "smoothness" - Smoothness-based selection
#'     \item "cv" - Cross-validation
#'     \item "multicriteria" - Combined multiple criteria
#'   }
#' @param scope Character, scope of selection. One of:
#'   \itemize{
#'     \item "global" - Estimate r from pooled data
#'     \item "clusterwise" - Estimate r separately for each cluster
#'     \item "hybrid" - Combine global and clusterwise estimates
#'   }
#' @param K Integer, number of clusters (required for clusterwise/hybrid scope)
#' @param rvec Integer vector, candidate r values to evaluate. Default: 1:6
#' @param ... Additional arguments passed to underlying selection functions
#'
#' @return A list containing:
#'   \item{best_r}{Selected optimal r value}
#'   \item{summary}{Data frame with evaluation results for all candidate r values}
#'   \item{best_model}{Fitted model with optimal r (if refit_final = TRUE)}
#'   \item{details}{Additional details from the selection method}
#'
#' @details
#' This function provides a simplified interface to the various r selection methods
#' available in the package. For more fine-grained control, use the specialized
#' functions directly (e.g., \code{select_r_by_vaf_mfa}, \code{select_r_mixture}).
#'
#' @examples
#' \dontrun{
#' # Simple VAF-based selection for MFA
#' result <- select_r(list_of_data, model_type = "MFA", method = "vaf",
#'                    scope = "global", K = 2, rvec = 1:5)
#'
#' # Multicriteria selection with clusterwise scope
#' result <- select_r(list_of_data, model_type = "MFA", method = "multicriteria",
#'                    scope = "clusterwise", K = 3, rvec = 2:6)
#'
#' # Cross-validation for MPCA
#' result <- select_r(list_of_data, model_type = "MPCA", method = "cv",
#'                    scope = "global", rvec = 1:4)
#' }
#'
#' @seealso \code{\link{select_r_mixture}}, \code{\link{select_r_by_vaf_mfa}},
#'   \code{\link{select_r_by_vss_mfa}}, \code{\link{select_r_by_smoothness_mfa}}
#'
#' @export
select_r <- function(
    list_of_data,
    model_type = c("MFA", "MPCA"),
    method = c("vaf", "vss", "smoothness", "cv", "multicriteria"),
    scope = c("global", "clusterwise", "hybrid"),
    K = NULL,
    rvec = 1:6,
    ...
) {
  # Input validation
  model_type <- match.arg(model_type)
  method <- match.arg(method)
  scope <- match.arg(scope)

  validate_list_of_data(list_of_data)

  # Check K requirement
  if (scope %in% c("clusterwise", "hybrid") && is.null(K)) {
    stop("K (number of clusters) is required for scope = '", scope, "'",
         call. = FALSE)
  }

  # Dispatch based on method and scope
  result <- switch(
    paste(method, scope, sep = "_"),

    # VAF-based methods
    "vaf_global" = {
      if (model_type == "MFA") {
        select_r_by_vaf_mfa(list_of_data, K = K, rvec = rvec, ...)
      } else {
        select_r_by_vaf_mpca(list_of_data, K = K, rvec = rvec, ...)
      }
    },
    "vaf_clusterwise" = ,
    "vaf_hybrid" = {
      select_r_mixture(
        list_of_data, K = K, rvec = rvec,
        mode = scope, model_type = model_type,
        criteria = "VAF", ...
      )
    },

    # VSS-based methods
    "vss_global" = {
      if (model_type == "MFA") {
        select_r_by_vss_mfa(list_of_data, K = K, rvec = rvec, ...)
      } else {
        select_r_by_vss_mpca(list_of_data, K = K, rvec = rvec, ...)
      }
    },
    "vss_clusterwise" = ,
    "vss_hybrid" = {
      select_r_mixture(
        list_of_data, K = K, rvec = rvec,
        mode = scope, model_type = model_type,
        criteria = "VSS", ...
      )
    },

    # Smoothness-based methods
    "smoothness_global" = {
      if (model_type == "MFA") {
        select_r_by_smoothness_mfa(list_of_data, K = K, rvec = rvec, ...)
      } else {
        select_r_by_smoothness_mpca(list_of_data, K = K, rvec = rvec, ...)
      }
    },
    "smoothness_clusterwise" = ,
    "smoothness_hybrid" = {
      select_r_mixture(
        list_of_data, K = K, rvec = rvec,
        mode = scope, model_type = model_type,
        criteria = "smoothness", ...
      )
    },

    # Cross-validation (global only)
    "cv_global" = {
      select_r_by_cv(list_of_data, rvec = rvec, model = model_type, ...)
    },
    "cv_clusterwise" = ,
    "cv_hybrid" = {
      stop("Cross-validation is only supported for scope = 'global'",
           call. = FALSE)
    },

    # Multicriteria methods
    "multicriteria_global" = {
      if (model_type == "MFA") {
        select_r_multicriteria_mfa(list_of_data, K = K, rvec = rvec, ...)
      } else {
        select_r_multicriteria_mpca(list_of_data, K = K, rvec = rvec, ...)
      }
    },
    "multicriteria_clusterwise" = ,
    "multicriteria_hybrid" = {
      select_r_mixture(
        list_of_data, K = K, rvec = rvec,
        mode = scope, model_type = model_type,
        criteria = c("VAF", "VSS", "smoothness"), ...
      )
    },

    # Default case
    stop("Unknown combination of method='", method, "' and scope='", scope, "'",
         call. = FALSE)
  )

  # Standardize output format
  if (!is.list(result)) {
    result <- list(best_r = result)
  }

  # Add metadata
  result$meta <- list(
    model_type = model_type,
    method = method,
    scope = scope,
    rvec = rvec
  )

  result
}


# ============================================================
# Legacy Interface with Deprecation Warnings
# ============================================================

#' @name select_r
#' @rdname select_r
#' @section Deprecated Functions:
#' The following functions are deprecated and will call \code{select_r()}
#' with appropriate parameters:
#' \itemize{
#'   \item \code{select_r_by_vaf_mfa} - Use \code{select_r(method="vaf", model_type="MFA")}
#'   \item \code{select_r_by_vaf_mpca} - Use \code{select_r(method="vaf", model_type="MPCA")}
#'   \item \code{select_r_by_vss_mfa} - Use \code{select_r(method="vss", model_type="MFA")}
#'   \item \code{select_r_by_vss_mpca} - Use \code{select_r(method="vss", model_type="MPCA")}
#' }
#' @keywords internal
NULL


#' Unified r Selection for Mixture Models (MFA/MPCA)
#'
#' This function provides a unified interface for selecting the optimal number of
#' factors (r) in mixture models. It supports three modes: "clusterwise" (evaluate r
#' separately per cluster), "global-first" (estimate r from pooled data first), and
#' "hybrid" (combine both approaches).
#'
#' @param list_of_data A list of N matrices, each (T_i x M).
#' @param K Number of clusters (fixed).
#' @param rvec Vector of candidate r values to evaluate. Default is 1:6.
#' @param mode Selection mode. Options:
#'   \itemize{
#'     \item "clusterwise" (default): Evaluate r separately for each cluster, then aggregate
#'     \item "global-first": Estimate r using global FA/PPCA on pooled data
#'     \item "hybrid": Use global r as initialization for clusterwise evaluation
#'   }
#' @param model_type Model type. Options: "MFA" (Mixture Factor Analysis, default),
#'   "MPCA" (Mixture Principal Component Analysis).
#' @param criteria Character vector specifying which criteria to use for clusterwise selection.
#'   Options: "VAF", "VSS", "smoothness". Default is c("VAF", "VSS", "smoothness").
#'   Only used when mode = "clusterwise" or "hybrid".
#' @param agg_method Method for aggregating cluster-specific r values. Options:
#'   "max" (conservative, takes maximum r across clusters),
#'   "mean" (average r, rounded),
#'   "min" (minimum r across clusters),
#'   "majority" (most common r value),
#'   "global" (use the global r estimate r_global as the final r; only valid in hybrid mode).
#'   Default is "max". In clusterwise mode, only "max", "mean", "min", "majority" are allowed.
#'   In hybrid mode, you may additionally use "global" to set r_final = r_global.
#'   Ignored when mode = "global-first".
#' @param global_model Model for global r estimation. Options: "PPCA" (default), "FA".
#'   Only used when mode = "global-first" or "hybrid".
#' @param global_criterion Criterion for global r selection. Options: "CV" (default),
#'   "ParallelAnalysis". Only used when mode = "global-first" or "hybrid".
#' @param global_folds Number of CV folds for global r selection. Default is 5.
#'   Only used when mode = "global-first" or "hybrid" with criterion = "CV".
#' @param global_seed Random seed for global r selection. Default is NULL.
#' @param weights Named vector of weights for combining criteria in clusterwise mode.
#'   Default is c(VAF=0.4, VSS=0.4, smoothness=0.2).
#' @param vss_threshold Threshold for VSS criterion. Default is 0.5.
#' @param vss_method Method for computing VSS. Options: "threshold", "ratio". Default is "threshold".
#' @param vss_rotate Rotation method for VSS. Options: "none", "varimax". Default is "none".
#' @param vss_normalize Method for normalizing VSS. Options: "logr", "max", "r", "r2", "none".
#'   Default is "logr".
#' @param max_iter Maximum EM iterations for model fitting.
#' @param nIterFA Sub-iterations for FA update (MFA only).
#' @param nIterPCA Sub-iterations for PCA update (MPCA only).
#' @param tol Convergence tolerance.
#' @param n_init Number of random initializations for mixture model fitting.
#' @param use_kmeans_init Whether to use k-means initialization.
#' @param subject_rdim_for_kmeans Dimension for k-means features.
#' @param mc_cores Number of cores for parallelization.
#' @param n_threads Number of OpenMP threads.
#' @param verbose Logical. If TRUE, print detailed diagnostic information.
#' @param refit_final Logical. If TRUE, refit the final model with selected r_final.
#'   Default is TRUE.
#'
#' @return A list with class "mixture_r_selection" containing:
#'   \item{r_final}{The final selected r value.}
#'   \item{r_global}{Global r estimate (if mode = "global-first" or "hybrid").}
#'   \item{r_cluster}{Named vector of cluster-specific r values (if mode = "clusterwise" or "hybrid").}
#'   \item{initial_fit}{Initial mixture model fit used for cluster assignments (if applicable).}
#'   \item{final_fit}{Final mixture model fit with r_final (if refit_final = TRUE).}
#'   \item{diagnostics}{List containing detailed diagnostic information.}
#'   \item{mode}{The selection mode used.}
#'   \item{model_type}{The model type used (MFA or MPCA).}
#'   \item{agg_method}{The aggregation method used (if applicable).}
#'
#' @details
#' This function provides a unified interface for r selection in mixture models,
#' integrating three complementary approaches:
#'
#' \strong{Mode: "global-first"}
#' \enumerate{
#'   \item Pool all data and fit a global FA/PPCA model
#'   \item Select r using cross-validation or parallel analysis
#'   \item Return r_global as the final estimate
#' }
#' This mode is fast and provides a cluster-independent baseline. It's recommended
#' when you want a simple, stable estimate that doesn't depend on cluster structure.
#'
#' \strong{Mode: "clusterwise"}
#' \enumerate{
#'   \item Fit initial mixture model with median(rvec) to get cluster assignments
#'   \item For each cluster, evaluate different r values using cluster-specific metrics
#'   \item Aggregate cluster-specific r values using specified method
#' }
#' This mode is useful when different clusters may require different numbers of factors.
#'
#' \strong{Mode: "hybrid"}
#' \enumerate{
#'   \item Estimate r_global using global FA/PPCA (as in "global-first")
#'   \item Use r_global as initialization for clusterwise evaluation
#'   \item Evaluate r values around r_global for each cluster
#'   \item Aggregate cluster-specific r values
#' }
#' This mode combines the stability of global estimation with the flexibility of
#' cluster-specific evaluation.
#'
#' @examples
#' \dontrun{
#' # Global-first mode (recommended for simplicity)
#' result <- select_r_mixture(
#'   list_of_data, K = 3, rvec = 2:6,
#'   mode = "global-first",
#'   model_type = "MFA"
#' )
#' print(result$r_final)
#'
#' # Clusterwise mode
#' result <- select_r_mixture(
#'   list_of_data, K = 3, rvec = 2:6,
#'   mode = "clusterwise",
#'   model_type = "MFA",
#'   agg_method = "max"
#' )
#'
#' # Hybrid mode
#' result <- select_r_mixture(
#'   list_of_data, K = 3, rvec = 2:6,
#'   mode = "hybrid",
#'   model_type = "MPCA"
#' )
#' }
#'
#' @export
select_r_mixture <- function(
    list_of_data,
    K,
    rvec = 1:6,
    mode = c("clusterwise", "global-first", "hybrid"),
    model_type = c("MFA", "MPCA"),
    criteria = c("VAF", "VSS", "smoothness"),
    agg_method = c("max", "mean", "min", "majority", "global"),
    global_model = c("PPCA", "FA"),
    global_criterion = c("CV", "ParallelAnalysis"),
    global_folds = 5,
    global_seed = NULL,
    weights = c(VAF = 0.4, VSS = 0.4, smoothness = 0.2),
    vss_threshold = 0.5,
    vss_method = c("threshold", "ratio"),
    vss_rotate = c("none", "varimax"),
    vss_normalize = c("logr", "max", "r", "r2", "none"),
    max_iter = 50,
    nIterFA = 5,
    nIterPCA = 5,
    tol = 1e-3,
    n_init = 1,
    use_kmeans_init = TRUE,
    subject_rdim_for_kmeans = 2,
    mc_cores = 1,
    n_threads = 1,
    verbose = FALSE,
    refit_final = TRUE
) {
  # Input validation
  mode <- match.arg(mode)
  model_type <- match.arg(model_type)
  agg_method <- match.arg(agg_method)
  global_model <- match.arg(global_model)
  global_criterion <- match.arg(global_criterion)
  vss_method <- match.arg(vss_method)
  vss_rotate <- match.arg(vss_rotate)
  vss_normalize <- match.arg(vss_normalize)
  criteria <- match.arg(criteria, several.ok = TRUE)
  
  # Validate agg_method compatibility with mode
  if (mode == "clusterwise" && agg_method == "global") {
    stop(
      'agg_method = "global" is not supported for mode = "clusterwise". ',
      'Use mode = "hybrid" with agg_method = "global", or choose one of ',
      '"max", "mean", "min", "majority".'
    )
  }
  
  if (length(rvec) < 2) {
    stop("rvec must contain at least 2 candidate r values")
  }
  
  N <- length(list_of_data)
  M <- ncol(list_of_data[[1]])
  
  if (verbose) {
    cat("=== Unified r Selection for Mixture Models ===\n")
    cat(sprintf("Model: %s, Mode: %s\n", model_type, mode))
    cat(sprintf("K = %d, N = %d subjects, M = %d channels\n", K, N, M))
    cat(sprintf("Candidate r values: %s\n", paste(rvec, collapse = ", ")))
  }
  
  # Initialize result structure
  result <- list(
    r_final = NULL,
    r_global = NULL,
    r_cluster = NULL,
    initial_fit = NULL,
    final_fit = NULL,
    diagnostics = list(),
    mode = mode,
    model_type = model_type,
    agg_method = agg_method
  )
  
  # Dispatch based on mode
  if (mode == "global-first") {
    # Mode 1: Global-first
    if (verbose) {
      cat(sprintf("\n=== Mode: global-first ===\n"))
      cat(sprintf("Estimating r using global %s with %s\n", global_model, global_criterion))
    }
    
    # Call appropriate global r selection function
    if (model_type == "MFA") {
      global_result <- select_r_global_mfa(
        list_of_data = list_of_data,
        rvec = rvec,
        model = global_model,
        criterion = global_criterion,
        folds = global_folds,
        seed = global_seed,
        verbose = verbose
      )
    } else {  # MPCA
      global_result <- select_r_global_mpca(
        list_of_data = list_of_data,
        rvec = rvec,
        model = global_model,
        criterion = global_criterion,
        folds = global_folds,
        seed = global_seed,
        verbose = verbose
      )
    }
    
    result$r_global <- global_result$r_global
    result$r_final <- global_result$r_global
    result$diagnostics$global_result <- global_result
    
    if (verbose) {
      cat(sprintf("\nGlobal r estimate: %d\n", result$r_final))
    }
    
  } else if (mode == "clusterwise") {
    # Mode 2: Clusterwise
    if (verbose) {
      cat(sprintf("\n=== Mode: clusterwise ===\n"))
      cat(sprintf("Evaluating r separately for each cluster\n"))
    }
    
    # Call appropriate clusterwise r selection function
    if (model_type == "MFA") {
      clusterwise_result <- select_r_clusterwise_mfa(
        list_of_data = list_of_data,
        K = K,
        rvec = rvec,
        criteria = criteria,
        agg_method = agg_method,
        weights = weights,
        vss_threshold = vss_threshold,
        vss_method = vss_method,
        vss_rotate = vss_rotate,
        vss_normalize = vss_normalize,
        max_iter = max_iter,
        nIterFA = nIterFA,
        tol = tol,
        n_init = n_init,
        use_kmeans_init = use_kmeans_init,
        subject_rdim_for_kmeans = subject_rdim_for_kmeans,
        mc_cores = mc_cores,
        n_threads = n_threads,
        verbose = verbose,
        mode = "clusterwise"
      )
    } else {  # MPCA
      clusterwise_result <- select_r_clusterwise_mpca(
        list_of_data = list_of_data,
        K = K,
        rvec = rvec,
        criteria = criteria,
        agg_method = agg_method,
        weights = weights,
        vss_threshold = vss_threshold,
        vss_method = vss_method,
        vss_rotate = vss_rotate,
        vss_normalize = vss_normalize,
        max_iter = max_iter,
        nIterPCA = nIterPCA,
        tol = tol,
        n_init = n_init,
        use_kmeans_init = use_kmeans_init,
        subject_rdim_for_kmeans = subject_rdim_for_kmeans,
        mc_cores = mc_cores,
        n_threads = n_threads,
        verbose = verbose,
        mode = "clusterwise"
      )
    }
    
    result$r_cluster <- clusterwise_result$cluster_r
    result$r_final <- clusterwise_result$agg_r
    result$initial_fit <- clusterwise_result$initial_fit
    result$diagnostics$cluster_metrics <- clusterwise_result$cluster_metrics
    
    if (verbose) {
      cat(sprintf("\nCluster-specific r values: %s\n", paste(result$r_cluster, collapse = ", ")))
      cat(sprintf("Aggregated r (method = %s): %d\n", agg_method, result$r_final))
    }
    
  } else if (mode == "hybrid") {
    # Mode 3: Hybrid
    if (verbose) {
      cat(sprintf("\n=== Mode: hybrid ===\n"))
      cat(sprintf("Step 1: Estimating global r using %s with %s\n", global_model, global_criterion))
    }
    
    # Step 1: Get global r estimate
    if (model_type == "MFA") {
      global_result <- select_r_global_mfa(
        list_of_data = list_of_data,
        rvec = rvec,
        model = global_model,
        criterion = global_criterion,
        folds = global_folds,
        seed = global_seed,
        verbose = verbose
      )
    } else {  # MPCA
      global_result <- select_r_global_mpca(
        list_of_data = list_of_data,
        rvec = rvec,
        model = global_model,
        criterion = global_criterion,
        folds = global_folds,
        seed = global_seed,
        verbose = verbose
      )
    }
    
    result$r_global <- global_result$r_global
    result$diagnostics$global_result <- global_result
    
    if (verbose) {
      cat(sprintf("Global r estimate: %d\n", result$r_global))
      cat(sprintf("\nStep 2: Refining with clusterwise evaluation\n"))
    }
    
    # Step 2: Use global r as initialization for clusterwise evaluation
    # For hybrid mode, we can either:
    # (a) Use r_global directly and skip clusterwise (if agg_method = "global")
    # (b) Use r_global as r_init for clusterwise evaluation
    
    if (agg_method == "global") {
      # Simply use global r
      result$r_final <- result$r_global
      
      if (verbose) {
        cat(sprintf("Using global r directly (agg_method = 'global'): %d\n", result$r_final))
      }
    } else {
      # Perform clusterwise evaluation with r_global as initialization
      # Pass r_global as r_init to use the global estimate for initial clustering
      cat(sprintf("\nUsing r_init = r_global (= %d) for clusterwise evaluation.\n", result$r_global))
      
      if (model_type == "MFA") {
        clusterwise_result <- select_r_clusterwise_mfa(
          list_of_data = list_of_data,
          K = K,
          rvec = rvec,
          criteria = criteria,
          agg_method = agg_method,
          weights = weights,
          vss_threshold = vss_threshold,
          vss_method = vss_method,
          vss_rotate = vss_rotate,
          vss_normalize = vss_normalize,
          max_iter = max_iter,
          nIterFA = nIterFA,
          tol = tol,
          n_init = n_init,
          use_kmeans_init = use_kmeans_init,
          subject_rdim_for_kmeans = subject_rdim_for_kmeans,
          mc_cores = mc_cores,
          n_threads = n_threads,
          verbose = verbose,
          mode = "clusterwise",
          r_init = result$r_global
        )
      } else {  # MPCA
        clusterwise_result <- select_r_clusterwise_mpca(
          list_of_data = list_of_data,
          K = K,
          rvec = rvec,
          criteria = criteria,
          agg_method = agg_method,
          weights = weights,
          vss_threshold = vss_threshold,
          vss_method = vss_method,
          vss_rotate = vss_rotate,
          vss_normalize = vss_normalize,
          max_iter = max_iter,
          nIterPCA = nIterPCA,
          tol = tol,
          n_init = n_init,
          use_kmeans_init = use_kmeans_init,
          subject_rdim_for_kmeans = subject_rdim_for_kmeans,
          mc_cores = mc_cores,
          n_threads = n_threads,
          verbose = verbose,
          mode = "clusterwise",
          r_init = result$r_global
        )
      }
      
      result$r_cluster <- clusterwise_result$cluster_r
      result$r_final <- clusterwise_result$agg_r
      result$initial_fit <- clusterwise_result$initial_fit
      result$diagnostics$cluster_metrics <- clusterwise_result$cluster_metrics
      
      if (verbose) {
        cat(sprintf("Cluster-specific r values: %s\n", paste(result$r_cluster, collapse = ", ")))
        cat(sprintf("Final r (hybrid, method = %s): %d\n", agg_method, result$r_final))
      }
    }
  }
  
  # Step 3: Optionally refit final model with selected r
  if (refit_final && !is.null(result$r_final)) {
    if (verbose) {
      cat(sprintf("\n=== Refitting final model with r = %d ===\n", result$r_final))
    }
    
    if (model_type == "MFA") {
      result$final_fit <- mfa_em_fit(
        list_of_data = list_of_data,
        K = K,
        r = result$r_final,
        max_iter = max_iter,
        nIterFA = nIterFA,
        tol = tol,
        n_init = n_init,
        use_kmeans_init = use_kmeans_init,
        subject_rdim_for_kmeans = subject_rdim_for_kmeans,
        mc_cores = mc_cores,
        n_threads = n_threads,
        verbose = verbose
      )
    } else {  # MPCA
      result$final_fit <- mixture_pca_em_fit(
        list_of_data = list_of_data,
        K = K,
        r = result$r_final,
        max_iter = max_iter,
        nIterPCA = nIterPCA,
        tol = tol,
        n_init = n_init,
        use_kmeans_init = use_kmeans_init,
        subject_rdim_for_kmeans = subject_rdim_for_kmeans,
        mc_cores = mc_cores,
        n_threads = n_threads,
        verbose = verbose
      )
    }
    
    if (verbose) {
      cat("Final model fit complete.\n")
    }
  }
  
  # Add class for S3 methods
  class(result) <- c("mixture_r_selection", "list")
  
  if (verbose) {
    cat(sprintf("\n=== Final selected r = %d ===\n", result$r_final))
  }
  
  return(result)
}


#' Print Method for mixture_r_selection Objects
#'
#' @param x A mixture_r_selection object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.mixture_r_selection <- function(x, ...) {
  cat("=== Mixture Model r Selection Results ===\n")
  cat(sprintf("Model Type: %s\n", x$model_type))
  cat(sprintf("Selection Mode: %s\n", x$mode))
  cat(sprintf("Final r: %d\n", x$r_final))
  
  if (!is.null(x$r_global)) {
    cat(sprintf("Global r: %d\n", x$r_global))
  }
  
  if (!is.null(x$r_cluster)) {
    cat(sprintf("Cluster-specific r: %s\n", paste(x$r_cluster, collapse = ", ")))
    cat(sprintf("Aggregation method: %s\n", x$agg_method))
  }
  
  if (!is.null(x$final_fit)) {
    cat("\nFinal model fitted: Yes\n")
    if (!is.null(x$final_fit$logLik)) {
      cat(sprintf("Log-likelihood: %.2f\n", x$final_fit$logLik))
    }
  } else {
    cat("\nFinal model fitted: No\n")
  }
  
  invisible(x)
}


#' Summary Method for mixture_r_selection Objects
#'
#' @param object A mixture_r_selection object.
#' @param ... Additional arguments (ignored).
#'
#' @export
summary.mixture_r_selection <- function(object, ...) {
  cat("=== Mixture Model r Selection Summary ===\n\n")
  
  cat("Model Configuration:\n")
  cat(sprintf("  Model Type: %s\n", object$model_type))
  cat(sprintf("  Selection Mode: %s\n", object$mode))
  cat(sprintf("  Final r: %d\n", object$r_final))
  
  if (!is.null(object$r_global)) {
    cat(sprintf("\nGlobal Estimation:\n"))
    cat(sprintf("  Global r: %d\n", object$r_global))
    if (!is.null(object$diagnostics$global_result)) {
      gr <- object$diagnostics$global_result
      cat(sprintf("  Criterion: %s\n", gr$criterion))
      cat(sprintf("  Model: %s\n", gr$model))
    }
  }
  
  if (!is.null(object$r_cluster)) {
    cat(sprintf("\nCluster-wise Evaluation:\n"))
    for (i in seq_along(object$r_cluster)) {
      cat(sprintf("  Cluster %d: r = %d\n", i, object$r_cluster[i]))
    }
    cat(sprintf("  Aggregation method: %s\n", object$agg_method))
  }
  
  if (!is.null(object$final_fit)) {
    cat("\nFinal Model:\n")
    cat("  Status: Fitted\n")
    if (!is.null(object$final_fit$logLik)) {
      cat(sprintf("  Log-likelihood: %.2f\n", object$final_fit$logLik))
    }
    if (!is.null(object$final_fit$z)) {
      cluster_sizes <- table(object$final_fit$z)
      cat(sprintf("  Cluster sizes: %s\n", paste(cluster_sizes, collapse = ", ")))
    }
  }
  
  invisible(object)
}


#' Plot Method for mixture_r_selection Objects
#'
#' Visualize r-selection results with multiple plot types.
#'
#' @param x A mixture_r_selection object.
#' @param type Type of plot to generate. Options:
#'   \itemize{
#'     \item "overview" (default): Multi-panel overview showing all available plots
#'     \item "global_cv": Global CV curve (only for hybrid/global-first modes with CV criterion)
#'     \item "cluster_metrics": Heatmap of cluster-specific metrics across r values
#'     \item "r_comparison": Bar plot comparing r_global and cluster-specific r values
#'     \item "cluster_scores": Line plot of composite scores for each cluster across r values
#'   }
#' @param ... Additional arguments (ignored).
#'
#' @details
#' The plot method provides several visualization types depending on the selection mode:
#'
#' \strong{For global-first mode:}
#' - "global_cv": Shows CV error curve with selected r_global
#'
#' \strong{For clusterwise mode:}
#' - "cluster_metrics": Heatmap of VAF/VSS/smoothness for each cluster and r value
#' - "cluster_scores": Composite scores showing how each cluster's metrics vary with r
#'
#' \strong{For hybrid mode:}
#' - All of the above plus "r_comparison" showing r_global vs cluster-specific r values
#'
#' The "overview" type automatically generates a multi-panel plot with all relevant
#' visualizations for the given mode.
#'
#' @examples
#' \dontrun{
#' result <- select_r_mixture(list_of_data, K = 3, mode = "hybrid")
#' plot(result)  # Overview
#' plot(result, type = "global_cv")
#' plot(result, type = "cluster_metrics")
#' plot(result, type = "r_comparison")
#' }
#'
#' @export
plot.mixture_r_selection <- function(x, type = c("overview", "global_cv", "cluster_metrics", 
                                                   "r_comparison", "cluster_scores"), ...) {
  type <- match.arg(type)
  
  # Determine which plots are available based on mode
  has_global <- !is.null(x$r_global) && !is.null(x$diagnostics$global_result)
  has_cluster <- !is.null(x$r_cluster) && !is.null(x$diagnostics$cluster_metrics)
  has_cv <- has_global && !is.null(x$diagnostics$global_result$cv_errors)
  
  if (type == "overview") {
    # Multi-panel overview
    n_plots <- sum(c(has_cv, has_cluster, has_cluster && has_global))
    
    if (n_plots == 0) {
      stop("No diagnostic data available for plotting. Run with verbose=TRUE or ensure diagnostics are stored.")
    }
    
    # Set up multi-panel layout
    if (n_plots == 1) {
      par(mfrow = c(1, 1))
    } else if (n_plots == 2) {
      par(mfrow = c(1, 2))
    } else {
      par(mfrow = c(2, 2))
    }
    
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    
    # Plot available panels
    if (has_cv) {
      .plot_global_cv(x)
    }
    if (has_cluster) {
      .plot_cluster_metrics(x)
    }
    if (has_cluster && has_global) {
      .plot_r_comparison(x)
    }
    if (has_cluster) {
      .plot_cluster_scores(x)
    }
    
  } else if (type == "global_cv") {
    if (!has_cv) {
      stop("Global CV data not available. This plot requires mode='global-first' or 'hybrid' with criterion='CV'.")
    }
    .plot_global_cv(x)
    
  } else if (type == "cluster_metrics") {
    if (!has_cluster) {
      stop("Cluster metrics not available. This plot requires mode='clusterwise' or 'hybrid'.")
    }
    .plot_cluster_metrics(x)
    
  } else if (type == "r_comparison") {
    if (!has_global || !has_cluster) {
      stop("r comparison requires both global and cluster data. Use mode='hybrid'.")
    }
    .plot_r_comparison(x)
    
  } else if (type == "cluster_scores") {
    if (!has_cluster) {
      stop("Cluster scores not available. This plot requires mode='clusterwise' or 'hybrid'.")
    }
    .plot_cluster_scores(x)
  }
  
  invisible(x)
}


# Internal plotting functions

.plot_global_cv <- function(x) {
  gr <- x$diagnostics$global_result
  
  if (is.null(gr$cv_errors)) {
    return(invisible(NULL))
  }
  
  rvec <- as.numeric(names(gr$cv_errors))
  errors <- gr$cv_errors
  
  plot(rvec, errors, type = "b", pch = 19, col = "steelblue", lwd = 2,
       xlab = "Number of factors (r)", ylab = "CV Error",
       main = sprintf("Global CV Curve (%s)\nSelected r = %d", gr$model, x$r_global),
       xaxt = "n")
  axis(1, at = rvec)
  
  # Mark selected r
  points(x$r_global, errors[as.character(x$r_global)], 
         pch = 19, col = "red", cex = 2)
  
  # Add grid
  grid(col = "gray80", lty = "dotted")
}


.plot_cluster_metrics <- function(x) {
  metrics_list <- x$diagnostics$cluster_metrics
  K <- length(metrics_list)
  
  if (K == 0 || is.null(metrics_list[[1]])) {
    return(invisible(NULL))
  }
  
  # Extract r values from first cluster
  rvec <- metrics_list[[1]]$r
  n_r <- length(rvec)
  
  # Create matrix for heatmap: rows = clusters, cols = r values
  # We'll use composite scores if available, otherwise VAF
  score_matrix <- matrix(NA, nrow = K, ncol = n_r)
  rownames(score_matrix) <- paste0("C", 1:K)
  colnames(score_matrix) <- paste0("r=", rvec)
  
  for (k in 1:K) {
    df <- metrics_list[[k]]
    if ("composite_score" %in% names(df)) {
      score_matrix[k, ] <- df$composite_score
    } else if ("VAF" %in% names(df)) {
      score_matrix[k, ] <- df$VAF
    }
  }
  
  # Create heatmap
  image(1:n_r, 1:K, t(score_matrix), 
        col = colorRampPalette(c("white", "yellow", "orange", "red"))(100),
        xlab = "Number of factors (r)", ylab = "Cluster",
        main = "Cluster-Specific Metrics Heatmap",
        xaxt = "n", yaxt = "n")
  axis(1, at = 1:n_r, labels = rvec)
  axis(2, at = 1:K, labels = 1:K, las = 1)
  
  # Add text values
  for (k in 1:K) {
    for (i in 1:n_r) {
      if (!is.na(score_matrix[k, i])) {
        text(i, k, sprintf("%.2f", score_matrix[k, i]), cex = 0.8)
      }
    }
  }
  
  # Mark selected r for each cluster
  if (!is.null(x$r_cluster)) {
    for (k in 1:K) {
      r_k <- x$r_cluster[k]
      r_idx <- which(rvec == r_k)
      if (length(r_idx) > 0) {
        points(r_idx, k, pch = 22, col = "blue", cex = 2, lwd = 2)
      }
    }
  }
}


.plot_r_comparison <- function(x) {
  K <- length(x$r_cluster)
  
  # Create bar plot comparing r values
  r_values <- c(x$r_global, x$r_cluster)
  names_vec <- c("Global", paste0("C", 1:K))
  colors_vec <- c("steelblue", rep("coral", K))
  
  barplot(r_values, names.arg = names_vec, col = colors_vec,
          ylab = "Selected r", main = "r Selection Comparison",
          ylim = c(0, max(r_values) * 1.2))
  
  # Add final r line
  abline(h = x$r_final, col = "darkgreen", lwd = 2, lty = 2)
  legend("topright", legend = c("Global", "Cluster", sprintf("Final (r=%d)", x$r_final)),
         fill = c("steelblue", "coral", NA), border = c("black", "black", NA),
         lty = c(NA, NA, 2), lwd = c(NA, NA, 2), col = c(NA, NA, "darkgreen"),
         bg = "white")
}


.plot_cluster_scores <- function(x) {
  metrics_list <- x$diagnostics$cluster_metrics
  K <- length(metrics_list)
  
  if (K == 0 || is.null(metrics_list[[1]])) {
    return(invisible(NULL))
  }
  
  # Extract r values
  rvec <- metrics_list[[1]]$r
  
  # Determine which metric to plot
  metric_name <- if ("composite_score" %in% names(metrics_list[[1]])) {
    "composite_score"
  } else if ("VAF" %in% names(metrics_list[[1]])) {
    "VAF"
  } else {
    return(invisible(NULL))
  }
  
  # Set up plot
  plot(NULL, xlim = range(rvec), ylim = c(0, 1),
       xlab = "Number of factors (r)", ylab = metric_name,
       main = sprintf("Cluster-Specific %s Curves", metric_name),
       xaxt = "n")
  axis(1, at = rvec)
  grid(col = "gray80", lty = "dotted")
  
  # Plot each cluster
  colors <- rainbow(K)
  for (k in 1:K) {
    df <- metrics_list[[k]]
    if (metric_name %in% names(df)) {
      lines(df$r, df[[metric_name]], col = colors[k], lwd = 2, type = "b", pch = 19)
    }
  }
  
  # Add legend
  legend("bottomright", legend = paste0("Cluster ", 1:K),
         col = colors, lwd = 2, pch = 19, bg = "white")
  
  # Mark selected r for each cluster
  if (!is.null(x$r_cluster)) {
    for (k in 1:K) {
      r_k <- x$r_cluster[k]
      df <- metrics_list[[k]]
      r_idx <- which(df$r == r_k)
      if (length(r_idx) > 0) {
        points(r_k, df[[metric_name]][r_idx], pch = 8, col = colors[k], cex = 2, lwd = 2)
      }
    }
  }
}
