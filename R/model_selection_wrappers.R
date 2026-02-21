#' Unified Model Selection for K and r
#'
#' This function provides a unified interface for selecting optimal K (number of
#' clusters) and r (factor dimension) for mixture models (MFA or MPCA). It supports
#' multiple selection strategies including grid search and two-stage selection.
#'
#' @param list_of_data A list of data matrices, each (T_i x M).
#' @param method Model type: "MFA" or "MPCA".
#' @param strategy Selection strategy: "grid" for full grid search, "twostage" for
#'   two-stage sequential selection (K first, then r).
#' @param Kvec Integer vector of candidate K values. Default is 1:5.
#' @param rvec Integer vector of candidate r values. Default is 1:5.
#' @param r_fixed Fixed r value for Stage 1 of twostage selection. Default is 2.
#' @param criterion Model selection criterion: "ICL" (Integrated Completed Likelihood)
#'   or "BIC" (Bayesian Information Criterion). Default is "ICL".
#' @param icl_max_iter Maximum iterations for ICL entropy calculation. Default is NULL
#'   (uses package default).
#' @param icl_temperature Temperature parameter for ICL entropy scaling. Higher values
#'   favor more clusters. Default is 2.0.
#' @param max_iter Maximum EM iterations. Default is 50.
#' @param nIterFA Sub-iterations for FA update (MFA only). Default is 5.
#' @param nIterPCA Sub-iterations for PCA update (MPCA only). Default is 5.
#' @param tol Convergence tolerance. Default is 1e-3.
#' @param mpca_method Method for MPCA: "EM" or "closed_form". Default is "EM".
#' @param n_init Number of random initializations. Default is 1.
#' @param use_kmeans_init Whether to use k-means initialization. Default is TRUE.
#' @param subject_rdim_for_kmeans Dimension for k-means features. Default is 2.
#' @param mc_cores_grid Number of cores for grid-level parallelization. Default is 1.
#' @param mc_cores Number of cores for initialization-level parallelization. Default is 1.
#' @param n_threads Number of OpenMP threads. Default is 1.
#' @param verbose Logical; if TRUE, print progress messages. Default is TRUE.
#'
#' @return A list with unified output structure:
#' \describe{
#'   \item{\code{best_K}}{Optimal number of clusters (scalar).}
#'   \item{\code{best_r}}{Optimal factor dimension (scalar).}
#'   \item{\code{best_model}}{Fitted model object with optimal (K, r).}
#'   \item{\code{summary_df}}{Data frame with selection summary. For grid strategy,
#'     contains all (K, r) combinations with BIC. For twostage, contains stage1 and
#'     stage2 summaries.}
#'   \item{\code{meta}}{List with metadata: method, strategy, Kvec, rvec, etc.}
#' }
#'
#' @details
#' This function wraps the existing model selection functions to provide a unified
#' interface for parameter sweep experiments. It handles:
#' \itemize{
#'   \item Parallel safety: Prevents nested parallelization by forcing n_threads=1
#'     when mc_cores_grid > 1
#'   \item Windows compatibility: Falls back to sequential execution on Windows
#'   \item Consistent output format: Returns the same structure regardless of
#'     method or strategy
#' }
#'
#' \strong{Strategy: "grid"}
#' Performs full grid search over all (K, r) combinations using BIC for selection.
#' This is thorough but computationally expensive.
#'
#' \strong{Strategy: "twostage"}
#' Two-stage sequential selection:
#' \enumerate{
#'   \item Stage 1: Fix r at r_fixed, search for optimal K using BIC
#'   \item Stage 2: Fix K at optimal value, search for optimal r using BIC
#' }
#' This is more stable and faster than full grid search.
#'
#' @examples
#' \dontrun{
#' # Grid search for MFA
#' result <- select_K_r(
#'   list_of_data, method = "MFA", strategy = "grid",
#'   Kvec = 2:4, rvec = 2:4
#' )
#' print(result$best_K)
#' print(result$best_r)
#'
#' # Two-stage selection for MPCA
#' result <- select_K_r(
#'   list_of_data, method = "MPCA", strategy = "twostage",
#'   Kvec = 2:4, rvec = 2:4, r_fixed = 2
#' )
#' }
#'
#' @seealso \code{\link{select_optimal_K_r_mfa}}, \code{\link{select_optimal_K_r_mpca}},
#'   \code{\link{select_optimal_K_r_mfa_twostage}}, \code{\link{select_optimal_K_r_mpca_twostage}}
#'
#' @export
select_K_r <- function(
    list_of_data,
    method = c("MFA", "MPCA"),
    strategy = c("grid", "twostage"),
    Kvec = 1:5,
    rvec = 1:5,
    r_fixed = 2,
    criterion = c("ICL", "BIC"),
    icl_max_iter = NULL,
    icl_temperature = 2.0,
    max_iter = 50,
    nIterFA = 5,
    nIterPCA = 5,
    tol = 1e-3,
    mpca_method = "EM",
    n_init = 1,
    use_kmeans_init = TRUE,
    subject_rdim_for_kmeans = 2,
    mc_cores_grid = 1,
    mc_cores = 1,
    n_threads = 1,
    verbose = TRUE
) {
  # Input validation
  method <- match.arg(method)
  strategy <- match.arg(strategy)
  criterion <- match.arg(criterion)

  stopifnot(
    "list_of_data must be a list" = is.list(list_of_data),
    "list_of_data cannot be empty" = length(list_of_data) > 0,
    "Kvec must be a numeric vector" = is.numeric(Kvec),
    "rvec must be a numeric vector" = is.numeric(rvec),
    "mc_cores_grid must be positive" = mc_cores_grid > 0,
    "mc_cores must be positive" = mc_cores > 0,
    "n_threads must be non-negative" = n_threads >= 0
  )
  
  # Parallel safety: prevent nested parallelization
  if (mc_cores_grid > 1 && n_threads > 1) {
    if (verbose) {
      message("Nested parallelism detected: forcing n_threads = 1 for stability.")
    }
    n_threads <- 1
  }
  
  # Windows fallback
  if (.Platform$OS.type == "windows" && mc_cores_grid > 1) {
    if (verbose) {
      warning("Windows detected: forcing sequential grid execution (mc_cores_grid = 1)")
    }
    mc_cores_grid <- 1
  }
  
  # Initialize result structure
  result <- list(
    best_K = NULL,
    best_r = NULL,
    best_model = NULL,
    summary_df = NULL,
    meta = list(
      method = method,
      strategy = strategy,
      criterion = criterion,
      Kvec = Kvec,
      rvec = rvec,
      r_fixed = if (strategy == "twostage") r_fixed else NULL
    )
  )
  
  # Dispatch based on method and strategy
  if (method == "MFA") {
    if (strategy == "grid") {
      # MFA grid search
      sel_result <- select_optimal_K_r_mfa(
        list_of_data = list_of_data,
        Kvec = Kvec,
        rvec = rvec,
        max_iter = max_iter,
        nIterFA = nIterFA,
        tol = tol,
        n_init = n_init,
        use_kmeans_init = use_kmeans_init,
        subject_rdim_for_kmeans = subject_rdim_for_kmeans,
        mc_cores_grid = mc_cores_grid,
        mc_cores = mc_cores,
        n_threads = n_threads
      )
      
      result$best_K <- sel_result$best_model_info$K
      result$best_r <- sel_result$best_model_info$r
      result$best_model <- sel_result$best_model_info$model
      result$summary_df <- sel_result$summary
      
    } else if (strategy == "twostage") {
      # MFA two-stage selection
      sel_result <- select_optimal_K_r_mfa_twostage(
        list_of_data = list_of_data,
        Kvec = Kvec,
        rvec = rvec,
        r_fixed = r_fixed,
        criterion = criterion,
        icl_max_iter = icl_max_iter,
        icl_temperature = icl_temperature,
        max_iter = max_iter,
        nIterFA = nIterFA,
        tol = tol,
        n_init = n_init,
        use_kmeans_init = use_kmeans_init,
        subject_rdim_for_kmeans = subject_rdim_for_kmeans,
        mc_cores = mc_cores,
        n_threads = n_threads
      )
      
      result$best_K <- sel_result$best_K
      result$best_r <- sel_result$best_r
      result$best_model <- sel_result$best_model
      result$summary_df <- list(
        stage1 = sel_result$stage1_summary,
        stage2 = sel_result$stage2_summary
      )
    }
    
  } else if (method == "MPCA") {
    if (strategy == "grid") {
      # MPCA grid search
      sel_result <- select_optimal_K_r_mpca(
        list_of_data = list_of_data,
        Kvec = Kvec,
        rvec = rvec,
        max_iter = max_iter,
        nIterPCA = nIterPCA,
        tol = tol,
        method = mpca_method,
        n_init = n_init,
        use_kmeans_init = use_kmeans_init,
        subject_rdim_for_kmeans = subject_rdim_for_kmeans,
        mc_cores_grid = mc_cores_grid,
        mc_cores = mc_cores,
        n_threads = n_threads
      )
      
      result$best_K <- sel_result$best_model_info$K
      result$best_r <- sel_result$best_model_info$r
      result$best_model <- sel_result$best_model_info$model
      result$summary_df <- sel_result$summary
      
    } else if (strategy == "twostage") {
      # MPCA two-stage selection
      sel_result <- select_optimal_K_r_mpca_twostage(
        list_of_data = list_of_data,
        Kvec = Kvec,
        rvec = rvec,
        r_fixed = r_fixed,
        criterion = criterion,
        icl_max_iter = icl_max_iter,
        icl_temperature = icl_temperature,
        max_iter = max_iter,
        nIterPCA = nIterPCA,
        tol = tol,
        method = mpca_method,
        n_init = n_init,
        use_kmeans_init = use_kmeans_init,
        subject_rdim_for_kmeans = subject_rdim_for_kmeans,
        mc_cores = mc_cores,
        n_threads = n_threads
      )
      
      result$best_K <- sel_result$best_K
      result$best_r <- sel_result$best_r
      result$best_model <- sel_result$best_model
      result$summary_df <- list(
        stage1 = sel_result$stage1_summary,
        stage2 = sel_result$stage2_summary
      )
    }
  }
  
  return(result)
}


#' Extract Selection Summary as Data Frame
#'
#' Converts the summary from \code{\link{select_K_r}} to a flat data frame,
#' useful for saving and plotting.
#'
#' @param select_result Result from \code{\link{select_K_r}}.
#'
#' @return A data frame with columns K, r, logLik, BIC, and optionally stage
#'   (for twostage strategy).
#'
#' @export
extract_selection_summary <- function(select_result) {
  summary_df <- select_result$summary_df
  strategy <- select_result$meta$strategy
  
  if (strategy == "grid") {
    # Grid search returns a single data frame
    return(summary_df)
  } else if (strategy == "twostage") {
    # Two-stage returns a list with stage1 and stage2
    stage1 <- summary_df$stage1
    stage2 <- summary_df$stage2
    
    stage1$stage <- "stage1_K_search"
    stage2$stage <- "stage2_r_search"
    
    combined <- rbind(stage1, stage2)
    return(combined)
  }
  
  return(summary_df)
}
