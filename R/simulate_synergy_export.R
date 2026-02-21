#' Generate v2 Dynamic Synergy Data for Model Selection Experiments
#'
#' This function provides a standardized interface for generating v2 dynamic
#' synergy data suitable for model selection experiments. It wraps
#' \code{\link{simulate_dynamic_synergy_data}} and returns data in a format
#' optimized for parameter sweep experiments with model selection.
#'
#' @param N Number of subjects (time series). Default is 30.
#' @param K Number of clusters (subgroups). Default is 3.
#' @param r Number of synergies per cluster. Default is 3.
#' @param M Number of observed channels (muscles). Default is 8.
#' @param T_each Time-series length per subject. Default is 200.
#' @param cluster_sep_spatial Numeric scalar (0-2) controlling spatial structure
#'   differentiation between clusters. Default is 1.0.
#' @param cluster_sep_temporal Numeric scalar (0-2) controlling temporal coordination
#'   differentiation between clusters. Default is 1.0.
#' @param cluster_sep_stability Numeric scalar (0-2) controlling stability-adaptability
#'   differentiation between clusters. Default is 1.0.
#' @param freq_range Numeric vector of length 2 specifying the range of base
#'   frequencies (in Hz) for the periodic activation signals. Default is c(0.5, 3.0).
#' @param amp_range Numeric vector of length 2 specifying the range of amplitudes
#'   for the periodic activation signals. Default is c(0.5, 1.5).
#' @param base_level Numeric scalar specifying the baseline activation level.
#'   Default is 0.05.
#' @param rectify_output Logical; if TRUE, apply half-wave rectification.
#'   Default is TRUE.
#' @param seed Random seed for reproducibility. Required for deterministic results.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{list_of_data}}{List of length N; each element is a (T_each x M) matrix
#'     of generated EMG signals. This is the standard input format for MFA/MPCA fitting.}
#'   \item{\code{true_cluster}}{Integer vector of length N indicating true cluster
#'     assignments (1 to K). Named as \code{true_cluster} for consistency.}
#'   \item{\code{true_Lambda_list}}{List of length K; each element is a (M x r) matrix
#'     of true cluster-specific synergy loadings. Essential for subspace recovery
#'     evaluation. This is the canonical name for ground truth basis matrices.}
#'   \item{\code{params}}{List containing all generation parameters for reproducibility.
#'     Includes both canonical short names (T, sep_spatial, sep_temporal, sep_stability)
#'     and full names (T_each, cluster_sep_spatial, etc.) for compatibility.
#'     Also includes r_mode = "global" for future extension.}
#'   \item{\code{Lambda_list}}{List of length K; each element is a (M x r) matrix
#'     of cluster-specific synergy loadings. Kept for backward compatibility;
#'     identical to \code{true_Lambda_list}.}
#'   \item{\code{C_list}}{List of length N; each element is a (T_each x r) matrix
#'     of synergy activation signals (for visualization/validation).}
#'   \item{\code{spatial_params}}{List containing spatial structure parameters
#'     for each cluster.}
#'   \item{\code{temporal_params}}{List containing temporal coordination parameters
#'     for each cluster.}
#'   \item{\code{stability_params}}{List containing stability-adaptability parameters
#'     for each cluster.}
#' }
#'
#' @details
#' This function is designed for use in parameter sweep experiments where model
#' selection (K, r) is performed on generated data. The key differences from
#' \code{\link{simulate_dynamic_synergy_data}} are:
#' \itemize{
#'   \item Returns \code{list_of_data} (renamed from \code{X_list}) for consistency
#'     with fitting functions
#'   \item Returns \code{true_cluster} (renamed from \code{z_true}) for consistency
#'     with evaluation functions
#'   \item Includes a \code{params} list for tracking generation parameters
#'   \item Seed is required (not optional) to ensure reproducibility
#' }
#'
#' @examples
#' \dontrun{
#' # Generate data for model selection experiment
#' data <- generate_synergy_data(
#'   N = 30, K = 3, r = 3, M = 8, T_each = 200,
#'   cluster_sep_spatial = 1.0,
#'   cluster_sep_temporal = 1.0,
#'   cluster_sep_stability = 1.0,
#'   seed = 123
#' )
#'
#' # Access the data for fitting
#' list_of_data <- data$list_of_data
#' true_cluster <- data$true_cluster
#'
#' # Fit MFA model
#' fit <- mfa_em_fit(list_of_data, K = 3, r = 3)
#'
#' # Evaluate clustering accuracy
#' ari <- mclust::adjustedRandIndex(true_cluster, fit$z)
#' }
#'
#' @seealso \code{\link{simulate_dynamic_synergy_data}} for the underlying
#'   simulation function, \code{\link{run_param_sweep_modelsel}} for parameter
#'   sweep experiments.
#'
#' @export
generate_synergy_data <- function(
    N = 30,
    K = 3,
    r = 3,
    M = 8,
    T_each = 200,
    cluster_sep_spatial = 1.0,
    cluster_sep_temporal = 1.0,
    cluster_sep_stability = 1.0,
    freq_range = c(0.5, 3.0),
    amp_range = c(0.5, 1.5),
    base_level = 0.05,
    rectify_output = TRUE,
    seed
) {
  # Input validation
  if (missing(seed)) {
    stop("seed is required for reproducibility in model selection experiments")
  }
  
  stopifnot(
    "N must be a positive integer" = is.numeric(N) && N > 0 && N == floor(N),
    "K must be a positive integer" = is.numeric(K) && K > 0 && K == floor(K),
    "r must be a positive integer" = is.numeric(r) && r > 0 && r == floor(r),
    "M must be a positive integer" = is.numeric(M) && M > 0 && M == floor(M),
    "T_each must be a positive integer" = is.numeric(T_each) && T_each > 0 && T_each == floor(T_each),
    "cluster_sep_spatial must be non-negative" = is.numeric(cluster_sep_spatial) && cluster_sep_spatial >= 0,
    "cluster_sep_temporal must be non-negative" = is.numeric(cluster_sep_temporal) && cluster_sep_temporal >= 0,
    "cluster_sep_stability must be non-negative" = is.numeric(cluster_sep_stability) && cluster_sep_stability >= 0,
    "seed must be numeric" = is.numeric(seed)
  )


  # Call the underlying simulation function

sim <- simulate_dynamic_synergy_data(
    N = N,
    K = K,
    r = r,
    M = M,
    T_each = T_each,
    cluster_sep_spatial = cluster_sep_spatial,
    cluster_sep_temporal = cluster_sep_temporal,
    cluster_sep_stability = cluster_sep_stability,
    freq_range = freq_range,
    amp_range = amp_range,
    base_level = base_level,
    rectify_output = rectify_output,
    seed = seed
  )

  # Build standardized output
  list(
    # Primary outputs for model fitting and evaluation
    list_of_data = sim$X_list,
    true_cluster = sim$z_true,
    
    # True basis matrices for subspace recovery evaluation (key for Step 3)
    true_Lambda_list = sim$Lambda_list,
    
    # Parameters for tracking/reproducibility
    # Includes both canonical short names and full names for compatibility
    params = list(
      N = N,
      K = K,
      r = r,
      M = M,
      T = T_each,                           # Canonical short name
      T_each = T_each,                      # Full name for compatibility
      sep_spatial = cluster_sep_spatial,    # Canonical short name
      sep_temporal = cluster_sep_temporal,  # Canonical short name
      sep_stability = cluster_sep_stability, # Canonical short name
      cluster_sep_spatial = cluster_sep_spatial,    # Full name for compatibility
      cluster_sep_temporal = cluster_sep_temporal,  # Full name for compatibility
      cluster_sep_stability = cluster_sep_stability, # Full name for compatibility
      seed = seed,
      r_mode = "global",                    # Fixed for future extension
      freq_range = freq_range,
      amp_range = amp_range,
      base_level = base_level,
      rectify_output = rectify_output
    ),
    
    # Additional outputs for visualization/validation
    # Lambda_list kept for backward compatibility with existing code
    Lambda_list = sim$Lambda_list,
    C_list = sim$C_list,
    spatial_params = sim$spatial_params,
    temporal_params = sim$temporal_params,
    stability_params = sim$stability_params
  )
}
