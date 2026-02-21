# Internal Helper: Fit All Six Methods on Data
#
# This helper function applies six analysis methods (SingleFA, SinglePCA,
# TwoStep_FA, TwoStep_PCA, MixtureFA, MixturePCA) to a given dataset.
# It is used by both compare_static_vs_dynamic() and compare_baselines().
#
# @param list_of_data List of N matrices, each (T_each x M)
# @param z_true Integer vector of length N with true cluster assignments
# @param N Number of subjects
# @param K Number of clusters
# @param r Number of synergies
# @param M Number of muscles
# @param seed Random seed for reproducibility
# @param n_init Number of random initializations for mixture methods
# @param use_kmeans_init Logical; if TRUE, also try k-means initialization
# @param mc_cores Number of cores for parallel computation
# @param method_pca Method for PCA fitting (default "EM")
#
# @return Data frame with columns: Method, BIC, ARI, SSE, VAF
fit_all_methods_on_data <- function(list_of_data, z_true, N, K, r, M,
                                    seed, n_init = 1, use_kmeans_init = FALSE,
                                    mc_cores = 1, method_pca = "EM", align_basis = TRUE,
                                    refine_assignments = TRUE,
                                    true_Lambda_list = NULL) {
  # Compute total sum of squares (SST) for VAF
  total_SST <- calc_total_SST(list_of_data)
  
  # 1) Single FA
  singleFA_res <- fit_single_factor_analysis(list_of_data, r = r)
  BIC_singleFA <- singleFA_res$BIC
  SSE_singleFA <- calc_reconstruction_error_singleFA(list_of_data, singleFA_res)
  VAF_singleFA <- 1 - (SSE_singleFA / total_SST)
  
  # 2) Single PCA
  singlePCA_res <- fit_single_pca(list_of_data, r = r)
  BIC_singlePCA <- singlePCA_res$BIC
  SSE_singlePCA <- calc_reconstruction_error_singlePCA(list_of_data, singlePCA_res)
  VAF_singlePCA <- 1 - (SSE_singlePCA / total_SST)
  
  # 3) Two-step FA
  seed_twostep <- seed + 2
  twoFA_res <- fit_two_step_baseline(
    list_of_data, 
    r = r, 
    K = K, 
    method = "FA",
    seed = seed_twostep,
    z_true = z_true,
    refine_assignments = refine_assignments,
    align_basis = align_basis
  )
  BIC_twoFA <- twoFA_res$BIC
  ARI_twoFA <- twoFA_res$ARI
  SSE_twoFA <- calc_reconstruction_error_twostep(list_of_data, twoFA_res)
  VAF_twoFA <- 1 - (SSE_twoFA / total_SST)
  
  # 4) Two-step PCA
  twoPCA_res <- fit_two_step_baseline(
    list_of_data, 
    r = r, 
    K = K, 
    method = "PCA",
    seed = seed_twostep,
    z_true = z_true,
    refine_assignments = refine_assignments,
    align_basis = align_basis
  )
  BIC_twoPCA <- twoPCA_res$BIC
  ARI_twoPCA <- twoPCA_res$ARI
  SSE_twoPCA <- calc_reconstruction_error_twostep(list_of_data, twoPCA_res)
  VAF_twoPCA <- 1 - (SSE_twoPCA / total_SST)
  
  # 5) Mixture FA
  seed_fit <- seed + 1
  mixtureFA_fit <- mfa_em_fit(
    list_of_data,
    K = K,
    r = r,
    max_iter = 50,
    nIterFA  = 20,
    tol      = 1e-4,
    n_init   = n_init,
    use_kmeans_init = use_kmeans_init,
    subject_rdim_for_kmeans = r,
    mc_cores = mc_cores,
    seed = seed_fit
  )
  ll_mixFA  <- compute_logLik_mfa(list_of_data, mixtureFA_fit)
  BIC_mixFA <- compute_BIC_mfa(ll_mixFA, K = K, r = r, M = M, N_obs = N)
  SSE_mixFA <- calc_reconstruction_error_mixtureFA(list_of_data, mixtureFA_fit)
  VAF_mixFA <- 1 - (SSE_mixFA / total_SST)
  
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("Package 'mclust' is required for ARI computation. Please install it.")
  }
  ari_mixFA <- mclust::adjustedRandIndex(mixtureFA_fit$z, z_true)
  
  # 6) Mixture PCA
  mixturePCA_fit <- mixture_pca_em_fit(
    list_of_data,
    K = K,
    r = r,
    max_iter  = 50,
    nIterPCA  = 20,
    tol       = 1e-4,
    method    = method_pca,
    n_init    = n_init,
    use_kmeans_init = use_kmeans_init,
    subject_rdim_for_kmeans = r,
    mc_cores = mc_cores,
    seed = seed_fit
  )
  ll_mixPCA  <- compute_logLik_mpca(list_of_data, mixturePCA_fit)
  BIC_mixPCA <- compute_BIC_mpca(ll_mixPCA, K = K, r = r, M = M, N_obs = N)
  SSE_mixPCA <- calc_reconstruction_error_mixturePCA(list_of_data, mixturePCA_fit)
  VAF_mixPCA <- 1 - (SSE_mixPCA / total_SST)
  ari_mixPCA <- mclust::adjustedRandIndex(mixturePCA_fit$z, z_true)
  
  # Compute subspace similarity for each method if true_Lambda_list is provided
  subsim_singleFA <- NA_real_
  subsim_singlePCA <- NA_real_
  subsim_twoFA <- NA_real_
  subsim_twoPCA <- NA_real_
  subsim_mixFA <- NA_real_
  subsim_mixPCA <- NA_real_
  
  if (!is.null(true_Lambda_list) && length(true_Lambda_list) > 0) {
    # For Single methods: compare single basis to each true cluster basis, weighted average
    # Single FA
    if (!is.null(singleFA_res$loadings)) {
      single_Lambda <- singleFA_res$loadings
      cluster_sizes <- table(z_true)
      weights <- as.numeric(cluster_sizes) / sum(cluster_sizes)
      sims <- sapply(seq_along(true_Lambda_list), function(k) {
        compute_subspace_similarity(true_Lambda_list[[k]], single_Lambda)
      })
      subsim_singleFA <- sum(weights * sims, na.rm = TRUE)
    }
    
    # Single PCA (uses P instead of loadings)
    if (!is.null(singlePCA_res$P)) {
      single_Lambda <- singlePCA_res$P
      cluster_sizes <- table(z_true)
      weights <- as.numeric(cluster_sizes) / sum(cluster_sizes)
      sims <- sapply(seq_along(true_Lambda_list), function(k) {
        compute_subspace_similarity(true_Lambda_list[[k]], single_Lambda)
      })
      subsim_singlePCA <- sum(weights * sims, na.rm = TRUE)
    }
    
    # TwoStep FA: use cluster-specific bases from cluster_models
    if (!is.null(twoFA_res$cluster_models) && !is.null(twoFA_res$z_est)) {
      # Extract W matrices from cluster_models (each element is list(W=..., n_members=...))
      twoFA_bases <- lapply(twoFA_res$cluster_models, function(cm) {
        if (!is.null(cm) && !is.null(cm$W)) cm$W else NULL
      })
      # Check if we have valid bases
      if (any(!sapply(twoFA_bases, is.null))) {
        alignment <- align_clusters_to_truth(z_true, twoFA_res$z_est)
        if (!is.na(alignment$accuracy)) {
          subspace_result <- compute_clusterwise_subspace_score(
            true_Lambda_list,
            twoFA_bases,
            alignment$mapping
          )
          subsim_twoFA <- subspace_result$mean_score
        }
      }
    }

    # TwoStep PCA: use cluster-specific bases from cluster_models
    if (!is.null(twoPCA_res$cluster_models) && !is.null(twoPCA_res$z_est)) {
      # Extract W matrices from cluster_models
      twoPCA_bases <- lapply(twoPCA_res$cluster_models, function(cm) {
        if (!is.null(cm) && !is.null(cm$W)) cm$W else NULL
      })
      # Check if we have valid bases
      if (any(!sapply(twoPCA_bases, is.null))) {
        alignment <- align_clusters_to_truth(z_true, twoPCA_res$z_est)
        if (!is.na(alignment$accuracy)) {
          subspace_result <- compute_clusterwise_subspace_score(
            true_Lambda_list,
            twoPCA_bases,
            alignment$mapping
          )
          subsim_twoPCA <- subspace_result$mean_score
        }
      }
    }
    
    # Mixture FA: use Lambda list
    if (!is.null(mixtureFA_fit$Lambda)) {
      alignment <- align_clusters_to_truth(z_true, mixtureFA_fit$z)
      if (!is.na(alignment$accuracy)) {
        subspace_result <- compute_clusterwise_subspace_score(
          true_Lambda_list,
          mixtureFA_fit$Lambda,
          alignment$mapping
        )
        subsim_mixFA <- subspace_result$mean_score
      }
    }
    
    # Mixture PCA: use W list or P list as fallback
    # (W may not be returned in all code paths, but P is always available)
    mixPCA_bases <- if (!is.null(mixturePCA_fit$W)) {
      mixturePCA_fit$W
    } else if (!is.null(mixturePCA_fit$P)) {
      mixturePCA_fit$P
    } else {
      NULL
    }
    if (!is.null(mixPCA_bases)) {
      alignment <- align_clusters_to_truth(z_true, mixturePCA_fit$z)
      if (!is.na(alignment$accuracy)) {
        subspace_result <- compute_clusterwise_subspace_score(
          true_Lambda_list,
          mixPCA_bases,
          alignment$mapping
        )
        subsim_mixPCA <- subspace_result$mean_score
      }
    }
  }
  
  # Compile results
  results_df <- data.frame(
    Method = c("SingleFA", "SinglePCA", "TwoStep_FA", "TwoStep_PCA", "MixtureFA", "MixturePCA"),
    BIC    = c(BIC_singleFA, BIC_singlePCA, BIC_twoFA, BIC_twoPCA, BIC_mixFA, BIC_mixPCA),
    ARI    = c(NA, NA, ARI_twoFA, ARI_twoPCA, ari_mixFA, ari_mixPCA),
    SSE    = c(SSE_singleFA, SSE_singlePCA, SSE_twoFA, SSE_twoPCA, SSE_mixFA, SSE_mixPCA),
    VAF    = c(VAF_singleFA, VAF_singlePCA, VAF_twoFA, VAF_twoPCA, VAF_mixFA, VAF_mixPCA),
    subspace_similarity = c(subsim_singleFA, subsim_singlePCA, subsim_twoFA, subsim_twoPCA, subsim_mixFA, subsim_mixPCA)
  )
  
  # Store fitted model objects as attributes for optional diagnostics
  attr(results_df, "mixtureFA_fit") <- mixtureFA_fit
  attr(results_df, "mixturePCA_fit") <- mixturePCA_fit
  attr(results_df, "singleFA_fit") <- singleFA_res
  attr(results_df, "singlePCA_fit") <- singlePCA_res
  attr(results_df, "twoFA_fit") <- twoFA_res
  attr(results_df, "twoPCA_fit") <- twoPCA_res
  
  # Return data frame (backward compatible)
  results_df
}

#' Compare Analysis Methods on Dynamic Synergy Data (v2)
#'
#' This function uses the v2 dynamic synergy simulator
#' (\code{\link{simulate_dynamic_synergy_data}}), which provides three
#' orthogonal control axes: spatial structure, temporal coordination, and
#' stability-adaptability. This allows for systematic evaluation of how
#' different motor control strategies affect muscle synergy structures.
#'
#' The function generates data using the Dynamic_v2 model and applies six
#' analysis methods (SingleFA, SinglePCA, TwoStep_FA, TwoStep_PCA, MixtureFA,
#' MixturePCA) to evaluate their performance.
#'
#' @param model_types Character vector specifying which model types to use.
#'   Default is \code{c("Dynamic_v2")}. Only "Dynamic_v2" is supported.
#' @param N Number of subjects (time series). Default is 30.
#' @param K True number of clusters. Default is 3.
#' @param r True factor dimension (number of synergies). Default is 3.
#' @param M Number of observed channels (muscles). Default is 8.
#' @param T_each Time-series length per subject. Default is 200.
#' @param cluster_sep_spatial Numeric scalar (0-2) controlling spatial structure
#'   differentiation between clusters. Higher values create more distinct spatial
#'   modules. Default is 1.0.
#' @param cluster_sep_temporal Numeric scalar (0-2) controlling temporal coordination
#'   differentiation between clusters. Higher values create more distinct temporal
#'   patterns. Default is 1.0.
#' @param cluster_sep_stability Numeric scalar (0-2) controlling stability-adaptability
#'   differentiation between clusters. Higher values create more distinct noise
#'   profiles. Default is 1.0.
#' @param seed Random seed for data generation. Default is 123.
#' @param output_dir Directory path where results will be saved. Default is
#'   \code{"results/compare_dynamic_v2"}. The directory will be created if it
#'   doesn't exist.
#' @param n_init Number of random initializations for mixture methods. Default is 1.
#' @param use_kmeans_init Logical; if TRUE, also try k-means initialization.
#'   Default is FALSE.
#' @param mc_cores Number of cores for parallel computation. Default is 1.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{\code{model_type}}{Character string: "Dynamic_v2"}
#'   \item{\code{Method}}{Analysis method name (e.g., "MixtureFA", "SinglePCA")}
#'   \item{\code{BIC}}{Bayesian Information Criterion}
#'   \item{\code{ARI}}{Adjusted Rand Index (clustering accuracy)}
#'   \item{\code{SSE}}{Sum of Squared Errors (reconstruction error)}
#'   \item{\code{VAF}}{Variance Accounted For (1 - SSE/totalSST)}
#'   \item{\code{seed}}{Random seed used for data generation}
#'   \item{\code{cluster_sep_spatial}}{Spatial separation parameter}
#'   \item{\code{cluster_sep_temporal}}{Temporal separation parameter}
#'   \item{\code{cluster_sep_stability}}{Stability separation parameter}
#'   \item{\code{N}}{Number of subjects}
#'   \item{\code{K}}{Number of clusters}
#'   \item{\code{r}}{Number of synergies}
#'   \item{\code{M}}{Number of muscles}
#'   \item{\code{T_each}}{Time series length}
#' }
#'
#' @details
#' The function saves results to an RDS file in the specified output directory
#' with filename format \code{results_YYYY-MM-DD.rds}.
#'
#' Data is generated using \code{\link{simulate_dynamic_synergy_data}},
#' which creates time-varying synergy activation patterns with three orthogonal
#' control axes: (1) Spatial structure, (2) Temporal coordination, and
#' (3) Stability-Adaptability.
#'
#' The generated data is analyzed using six methods, allowing comparison of
#' method performance on realistic dynamic synergy data.
#'
#' @examples
#' \dontrun{
#' # Basic comparison with v2 simulator
#' results <- compare_static_vs_dynamic(
#'   N = 30, K = 3, r = 3, M = 8, T_each = 200,
#'   cluster_sep_spatial = 1.0,
#'   cluster_sep_temporal = 1.0,
#'   cluster_sep_stability = 1.0,
#'   seed = 123
#' )
#' print(results)
#'
#' # Compare with high spatial separation
#' results_spatial <- compare_static_vs_dynamic(
#'   N = 30, K = 3, r = 3, M = 8, T_each = 200,
#'   cluster_sep_spatial = 2.0,
#'   cluster_sep_temporal = 0.5,
#'   cluster_sep_stability = 1.0,
#'   seed = 456
#' )
#'
#' # Visualize results
#' library(ggplot2)
#' ggplot(results, aes(x = Method, y = ARI)) +
#'   geom_bar(stat = "identity", fill = "steelblue") +
#'   labs(title = "Method Performance on Dynamic v2 Data")
#' }
#'
#' @export
compare_static_vs_dynamic <- function(
    model_types = c("Dynamic_v2"),
    N = 30,
    K = 3,
    r = 3,
    M = 8,
    T_each = 200,
    cluster_sep_spatial = 1.0,
    cluster_sep_temporal = 1.0,
    cluster_sep_stability = 1.0,
    seed = 123,
    output_dir = "results/compare_dynamic_v2",
    n_init = 1,
    use_kmeans_init = FALSE,
    mc_cores = 1
) {
  # Create output directory if it doesn't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Initialize results list
  results <- list()
  
  # Loop over each model type
  for (model_type in model_types) {
    message("Running simulation for: ", model_type)
    
    # Generate data using dynamic synergy model v2 with three orthogonal axes
    if (model_type == "Dynamic_v2") {
      sim_data <- simulate_dynamic_synergy_data(
        N = N,
        K = K,
        r = r,
        M = M,
        T_each = T_each,
        cluster_sep_spatial = cluster_sep_spatial,
        cluster_sep_temporal = cluster_sep_temporal,
        cluster_sep_stability = cluster_sep_stability,
        seed = seed
      )
      # simulate_dynamic_synergy_data returns X_list instead of list_of_data
      list_of_data <- sim_data$X_list
      z_true <- sim_data$z_true
      
    } else {
      stop("Unknown model_type: ", model_type, 
           ". Must be 'Dynamic_v2'.")
    }
    
    # Now apply the 6 methods to the generated data using the helper function
    comp <- fit_all_methods_on_data(
      list_of_data = list_of_data,
      z_true = z_true,
      N = N,
      K = K,
      r = r,
      M = M,
      seed = seed,
      n_init = n_init,
      use_kmeans_init = use_kmeans_init,
      mc_cores = mc_cores,
      method_pca = "EM"
    )
    
    # Add model_type, seed, and parameter columns
    comp$model_type <- model_type
    comp$seed <- seed
    comp$cluster_sep_spatial <- cluster_sep_spatial
    comp$cluster_sep_temporal <- cluster_sep_temporal
    comp$cluster_sep_stability <- cluster_sep_stability
    comp$N <- N
    comp$K <- K
    comp$r <- r
    comp$M <- M
    comp$T_each <- T_each
    
    # Store results
    results[[model_type]] <- comp
  }
  
  # Combine results into a single data frame
  df <- do.call(rbind, results)
  rownames(df) <- NULL
  
  # Save results to file
  output_file <- file.path(output_dir, paste0("results_", Sys.Date(), ".rds"))
  saveRDS(df, output_file)
  message("Results saved to: ", output_file)
  
  # Return the combined data frame
  return(df)
}

#' Compare Analysis Methods on v2 Dynamic Synergy Data
#'
#' This function enables pure method-level benchmarking by generating data using
#' only the v2 dynamic synergy simulator and applying all six analysis methods
#' (SingleFA, SinglePCA, TwoStep_FA, TwoStep_PCA, MixtureFA, MixturePCA).
#' This complements \code{\link{compare_static_vs_dynamic}} by focusing on
#' method performance rather than model comparison.
#'
#' @param N Number of subjects (time series). Default is 50.
#' @param K Number of clusters. Default is 2.
#' @param r Number of synergies (factor dimension). Default is 2.
#' @param M Number of observed channels (muscles). Default is 6.
#' @param T_each Time-series length per subject. Default is 100.
#' @param cluster_sep_spatial Numeric scalar (0-2) controlling spatial structure
#'   differentiation between clusters. Higher values create more distinct spatial
#'   modules. Default is 1.0.
#' @param cluster_sep_temporal Numeric scalar (0-2) controlling temporal coordination
#'   differentiation between clusters. Higher values create more distinct temporal
#'   patterns. Default is 1.0.
#' @param cluster_sep_stability Numeric scalar (0-2) controlling stability-adaptability
#'   differentiation between clusters. Higher values create more distinct noise
#'   profiles. Default is 1.0.
#' @param seed Random seed for data generation. Default is 123.
#' @param n_init Number of random initializations for mixture methods. Default is 1.
#' @param use_kmeans_init Logical; if TRUE, also try k-means initialization.
#'   Default is FALSE.
#' @param mc_cores Number of cores for parallel computation. Default is 1.
#' @param method_pca Method for PCA fitting. Default is "EM".
#' @param align_basis Logical. If \code{TRUE} (default), applies Hungarian algorithm
#'   alignment in Step 1.5 of TwoStep methods to ensure corresponding components across
#'   subjects are in the same order before clustering. Set to \code{FALSE} to disable
#'   alignment and compare performance with/without this step.
#' @param refine_assignments Logical. If \code{TRUE} (default), performs a single
#'   reassignment pass after fitting cluster-specific models in Step 4 of TwoStep methods.
#'   Each subject is reassigned to the cluster that minimizes reconstruction error. Set to
#'   \code{FALSE} to disable reassignment and use only the initial k-means cluster assignments.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{\code{Method}}{Analysis method name (e.g., "MixtureFA", "SinglePCA")}
#'   \item{\code{BIC}}{Bayesian Information Criterion}
#'   \item{\code{ARI}}{Adjusted Rand Index (clustering accuracy)}
#'   \item{\code{SSE}}{Sum of Squared Errors (reconstruction error)}
#'   \item{\code{VAF}}{Variance Accounted For (1 - SSE/totalSST)}
#'   \item{\code{model_type}}{Always "Dynamic_v2" for this function}
#'   \item{\code{seed}}{Random seed used for data generation}
#'   \item{\code{cluster_sep_spatial}}{Spatial separation parameter}
#'   \item{\code{cluster_sep_temporal}}{Temporal separation parameter}
#'   \item{\code{cluster_sep_stability}}{Stability separation parameter}
#'   \item{\code{N}}{Number of subjects}
#'   \item{\code{K}}{Number of clusters}
#'   \item{\code{r}}{Number of synergies}
#'   \item{\code{M}}{Number of muscles}
#'   \item{\code{T_each}}{Time series length}
#' }
#'
#' @details
#' This function generates data using \code{\link{simulate_dynamic_synergy_data}}
#' with three orthogonal control axes (spatial, temporal, stability), then applies
#' all six analysis methods to evaluate their performance. This is useful for:
#' \itemize{
#'   \item Method benchmarking on realistic dynamic synergy data
#'   \item Parameter sensitivity analysis for method performance
#'   \item Integration with \code{\link{run_param_sweep}} for systematic evaluation
#' }
#'
#' The function complements \code{\link{compare_static_vs_dynamic}} by focusing
#' on method performance rather than model comparison. Together, these functions
#' provide comprehensive evaluation of both models and methods.
#'
#' @examples
#' \dontrun{
#' # Basic method comparison on v2 data
#' results <- compare_baselines(
#'   N = 50, K = 2, r = 2, M = 6, T_each = 100,
#'   cluster_sep_spatial = 1.0,
#'   cluster_sep_temporal = 1.0,
#'   cluster_sep_stability = 1.0,
#'   seed = 123
#' )
#' print(results)
#'
#' # Compare methods with high spatial separation
#' results_spatial <- compare_baselines(
#'   N = 50, K = 2, r = 2, M = 6, T_each = 100,
#'   cluster_sep_spatial = 2.0,
#'   cluster_sep_temporal = 0.5,
#'   cluster_sep_stability = 1.0,
#'   seed = 456
#' )
#'
#' # Visualize method performance
#' library(ggplot2)
#' ggplot(results, aes(x = Method, y = ARI)) +
#'   geom_bar(stat = "identity", fill = "steelblue") +
#'   labs(title = "Method Performance on v2 Dynamic Data")
#' }
#'
#' @seealso
#' \code{\link{compare_static_vs_dynamic}} for model comparison,
#' \code{\link{run_param_sweep}} for systematic parameter sweep,
#' \code{\link{simulate_dynamic_synergy_data}} for data generation
#'
#' @export
compare_baselines <- function(
    N = 50,
    K = 2,
    r = 2,
    M = 6,
    T_each = 100,
    cluster_sep_spatial = 1.0,
    cluster_sep_temporal = 1.0,
    cluster_sep_stability = 1.0,
    seed = 123,
    n_init = 1,
    use_kmeans_init = FALSE,
    mc_cores = 1,
    method_pca = "EM",
    align_basis = TRUE,
    refine_assignments = TRUE
) {
  message("Generating v2 dynamic synergy data...")
  
  # Generate data using v2 simulator
  sim_data <- simulate_dynamic_synergy_data(
    N = N,
    K = K,
    r = r,
    M = M,
    T_each = T_each,
    cluster_sep_spatial = cluster_sep_spatial,
    cluster_sep_temporal = cluster_sep_temporal,
    cluster_sep_stability = cluster_sep_stability,
    seed = seed
  )
  
  # Extract data and ground truth
  list_of_data <- sim_data$X_list
  z_true <- sim_data$z_true
  true_Lambda_list <- sim_data$Lambda_list
  
  message("Fitting all methods...")
  
  # Apply all six methods using the helper function
  comp <- fit_all_methods_on_data(
    list_of_data = list_of_data,
    z_true = z_true,
    N = N,
    K = K,
    r = r,
    M = M,
    seed = seed,
    n_init = n_init,
    use_kmeans_init = use_kmeans_init,
    mc_cores = mc_cores,
    method_pca = method_pca,
    align_basis = align_basis,
    refine_assignments = refine_assignments,
    true_Lambda_list = true_Lambda_list
  )
  
  # Add metadata columns
  comp$model_type <- "Dynamic_v2"
  comp$seed <- seed
  comp$cluster_sep_spatial <- cluster_sep_spatial
  comp$cluster_sep_temporal <- cluster_sep_temporal
  comp$cluster_sep_stability <- cluster_sep_stability
  comp$N <- N
  comp$K <- K
  comp$r <- r
  comp$M <- M
  comp$T_each <- T_each
  
  message("Method comparison complete!")
  
  # Return the results
  return(comp)
}


#' Check SSE Diagnostics for Mixture Models
#'
#' This function performs diagnostic checks to investigate SSE (Sum of Squared Errors)
#' behavior in mixture model fits, particularly when SSE increases monotonically with
#' latent dimension r. It examines three key aspects: EM convergence quality, noise
#' variance (Psi) distribution, and VAF-SSE correlation.
#'
#' @param fit_results A data frame returned by \code{fit_all_methods_on_data},
#'   containing columns: Method, BIC, ARI, SSE, VAF. The fitted model objects may be
#'   stored as attributes (mixtureFA_fit, mixturePCA_fit) if available.
#' @param list_of_data The original data used for fitting (list of N matrices, each T_i x M).
#'   Required for VAF-SSE correlation theoretical line calculation.
#' @param mixtureFA_fit Optional fitted MFA model object. If NULL, will attempt to extract
#'   from fit_results attributes. If not available, EM convergence and Psi diagnostics
#'   for MFA will be skipped.
#' @param mixturePCA_fit Optional fitted MPCA model object. If NULL, will attempt to extract
#'   from fit_results attributes. If not available, EM convergence and Psi diagnostics
#'   for MPCA will be skipped.
#' @param output_dir Optional directory path where diagnostic plots will be saved.
#'   If NULL (default), plots are displayed on screen only.
#' @param psi_small_thresh Threshold for detecting overfitting (Psi values below this).
#'   Default is 0.01.
#' @param psi_large_thresh Threshold for detecting poor fit (Psi values above this).
#'   Default is 100.
#'
#' @return A list with diagnostic results:
#' \describe{
#'   \item{\code{convergence_ok}}{Logical; TRUE if convergence appears stable}
#'   \item{\code{psi_ok}}{Logical; TRUE if Psi values are in reasonable range}
#'   \item{\code{correlation_ok}}{Logical; TRUE if VAF-SSE correlation is strongly negative}
#'   \item{\code{warnings}}{Character vector of warning messages}
#'   \item{\code{psi_summary_mfa}}{Summary statistics for MFA Psi values}
#'   \item{\code{psi_summary_mpca}}{Summary statistics for MPCA Psi values}
#'   \item{\code{vaf_sse_correlation}}{Correlation coefficient between VAF and SSE}
#' }
#'
#' @details
#' The function performs three main diagnostic checks:
#'
#' \strong{1. EM Convergence Quality Check}
#'
#' Since the C++ EM implementation doesn't return a convergence trace, this check
#' examines the final model quality by comparing log-likelihoods and checking for
#' degenerate solutions (e.g., empty clusters, extreme parameter values).
#'
#' \strong{2. Psi (Noise Variance) Distribution Check}
#'
#' Examines the distribution of noise variance parameters in both MFA and MPCA models.
#' Extreme values (very small < 0.01 or very large > 100) may indicate:
#' \itemize{
#'   \item Overfitting (Psi approaches 0)
#'   \item Poor model fit (Psi approaches infinity)
#'   \item Numerical instability
#' }
#'
#' Creates histograms showing the distribution of Psi values across all clusters
#' and muscles.
#'
#' \strong{3. VAF-SSE Correlation Verification}
#'
#' Since VAF = 1 - SSE/SST, there should be a strong negative correlation (approximately -1)
#' between VAF and SSE across methods. Deviations from this expected relationship
#' may indicate:
#' \itemize{
#'   \item Scaling inconsistencies in SSE calculation
#'   \item Different SST values being used across methods
#'   \item Numerical precision issues
#' }
#'
#' Creates a scatter plot of VAF vs SSE with the expected theoretical relationship.
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic usage with fit_all_methods_on_data()
#' # The fitted models are stored as attributes and automatically extracted
#' sim_data <- generate_synergy_data(N = 30, K = 3, r = 3, M = 8, T_each = 200, seed = 123)
#' fit_results <- fit_all_methods_on_data(
#'   list_of_data = sim_data$list_of_data,
#'   z_true = sim_data$true_cluster,
#'   N = 30, K = 3, r = 3, M = 8,
#'   seed = 123
#' )
#'
#' # Run diagnostics (fitted models extracted from attributes automatically)
#' diagnostics <- check_SSE_diagnostics(
#'   fit_results = fit_results,
#'   list_of_data = sim_data$list_of_data,
#'   output_dir = "diagnostics"
#' )
#'
#' # Check results
#' print(diagnostics$warnings)
#' print(diagnostics$psi_summary_mfa)
#' print(diagnostics$vaf_sse_correlation)
#'
#' # Example 2: Pass fitted models explicitly
#' # Useful when you have fitted models from other sources
#' mfa_fit <- mfa_em_fit(sim_data$list_of_data, K = 3, r = 3)
#' mpca_fit <- mixture_pca_em_fit(sim_data$list_of_data, K = 3, r = 3)
#'
#' # Create a minimal results data frame
#' results_df <- data.frame(
#'   Method = c("MixtureFA", "MixturePCA"),
#'   SSE = c(100, 120),
#'   VAF = c(0.85, 0.80)
#' )
#'
#' # Run diagnostics with explicit fitted models
#' diagnostics2 <- check_SSE_diagnostics(
#'   fit_results = results_df,
#'   list_of_data = sim_data$list_of_data,
#'   mixtureFA_fit = mfa_fit,
#'   mixturePCA_fit = mpca_fit,
#'   output_dir = "diagnostics2",
#'   psi_small_thresh = 0.005,  # Custom threshold
#'   psi_large_thresh = 200     # Custom threshold
#' )
#'
#' # Example 3: VAF-SSE correlation only (no fitted models)
#' # Useful for quick checks without EM convergence/Psi diagnostics
#' diagnostics3 <- check_SSE_diagnostics(
#'   fit_results = fit_results,
#'   list_of_data = sim_data$list_of_data
#' )
#' # Note: EM convergence and Psi diagnostics will be skipped with a warning
#' }
#'
#' @export
check_SSE_diagnostics <- function(fit_results, 
                                   list_of_data, 
                                   mixtureFA_fit = NULL,
                                   mixturePCA_fit = NULL,
                                   output_dir = NULL,
                                   psi_small_thresh = 0.01,
                                   psi_large_thresh = 100) {
  # Initialize warnings vector
  warnings <- character(0)
  
  # Create output directory if specified
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # ============================================================================
  # INPUT VALIDATION AND EXTRACTION
  # ============================================================================
  
  # Validate fit_results input
  if (!is.data.frame(fit_results) && !is.list(fit_results)) {
    stop("fit_results must be either a data frame or a list")
  }
  
  # Extract VAF and SSE values based on input type
  if (is.data.frame(fit_results)) {
    # Standard case: fit_results is a data frame from fit_all_methods_on_data()
    if (!all(c("Method", "SSE", "VAF") %in% names(fit_results))) {
      stop("fit_results data frame must contain columns: Method, SSE, VAF")
    }
    vaf_values <- fit_results$VAF
    sse_values <- fit_results$SSE
    
    # Try to extract fitted models from attributes if not provided explicitly
    if (is.null(mixtureFA_fit)) {
      mixtureFA_fit <- attr(fit_results, "mixtureFA_fit")
    }
    if (is.null(mixturePCA_fit)) {
      mixturePCA_fit <- attr(fit_results, "mixturePCA_fit")
    }
  } else if (is.list(fit_results)) {
    # Legacy support: fit_results is a list with $results data frame
    if (!is.null(fit_results$results) && is.data.frame(fit_results$results)) {
      vaf_values <- fit_results$results$VAF
      sse_values <- fit_results$results$SSE
      
      # Extract fitted models from list if not provided explicitly
      if (is.null(mixtureFA_fit) && !is.null(fit_results$mixtureFA_fit)) {
        mixtureFA_fit <- fit_results$mixtureFA_fit
      }
      if (is.null(mixturePCA_fit) && !is.null(fit_results$mixturePCA_fit)) {
        mixturePCA_fit <- fit_results$mixturePCA_fit
      }
    } else if (!is.null(fit_results$VAF) && !is.null(fit_results$SSE)) {
      # Direct VAF/SSE vectors in list
      vaf_values <- fit_results$VAF
      sse_values <- fit_results$SSE
    } else {
      stop("fit_results list must contain either $results data frame or $VAF and $SSE vectors")
    }
  }
  
  # ============================================================================
  # MODEL VALIDATION
  # ============================================================================
  
  # Validate mixtureFA_fit structure if provided
  if (!is.null(mixtureFA_fit)) {
    if (!is.list(mixtureFA_fit)) {
      warnings <- c(warnings, "mixtureFA_fit is not a list - skipping MFA diagnostics")
      mixtureFA_fit <- NULL
    } else {
      # Check for essential fields
      required_fields_mfa <- c("z", "Psi", "logLik")
      missing_fields <- setdiff(required_fields_mfa, names(mixtureFA_fit))
      if (length(missing_fields) > 0) {
        warnings <- c(warnings, sprintf(
          "mixtureFA_fit missing essential fields (%s) - skipping MFA diagnostics",
          paste(missing_fields, collapse = ", ")
        ))
        mixtureFA_fit <- NULL
      }
    }
  }
  
  # Validate mixturePCA_fit structure if provided
  if (!is.null(mixturePCA_fit)) {
    if (!is.list(mixturePCA_fit)) {
      warnings <- c(warnings, "mixturePCA_fit is not a list - skipping MPCA diagnostics")
      mixturePCA_fit <- NULL
    } else {
      # Check for essential fields
      required_fields_mpca <- c("z", "sigma2", "logLik")
      missing_fields <- setdiff(required_fields_mpca, names(mixturePCA_fit))
      if (length(missing_fields) > 0) {
        warnings <- c(warnings, sprintf(
          "mixturePCA_fit missing essential fields (%s) - skipping MPCA diagnostics",
          paste(missing_fields, collapse = ", ")
        ))
        mixturePCA_fit <- NULL
      }
    }
  }
  
  # Warn if fitted models are not available after validation
  if (is.null(mixtureFA_fit) && is.null(mixturePCA_fit)) {
    warnings <- c(warnings, 
                  "No valid fitted mixture models available - EM convergence and Psi diagnostics will be skipped. ",
                  "To enable full diagnostics, pass valid mixtureFA_fit and/or mixturePCA_fit parameters.")
  }
  
  # ============================================================================
  # 1. EM CONVERGENCE QUALITY CHECK
  # ============================================================================
  
  convergence_ok <- TRUE
  
  # Check MFA convergence quality
  if (!is.null(mixtureFA_fit)) {
    # Check for empty clusters
    cluster_sizes <- table(mixtureFA_fit$z)
    if (any(cluster_sizes == 0)) {
      warnings <- c(warnings, "MFA: Empty clusters detected - EM may not have converged properly")
      convergence_ok <- FALSE
    }
    
    # Check for very small clusters (< 5% of data)
    N <- length(mixtureFA_fit$z)
    if (any(cluster_sizes < 0.05 * N)) {
      warnings <- c(warnings, sprintf("MFA: Very small clusters detected (min size: %d / %d subjects)", 
                                     min(cluster_sizes), N))
    }
    
    # Check if log-likelihood is finite
    if (!is.finite(mixtureFA_fit$logLik)) {
      warnings <- c(warnings, "MFA: Non-finite log-likelihood - numerical issues detected")
      convergence_ok <- FALSE
    }
  }
  
  # Check MPCA convergence quality
  if (!is.null(mixturePCA_fit)) {
    # Check for empty clusters
    cluster_sizes <- table(mixturePCA_fit$z)
    if (any(cluster_sizes == 0)) {
      warnings <- c(warnings, "MPCA: Empty clusters detected - EM may not have converged properly")
      convergence_ok <- FALSE
    }
    
    # Check for very small clusters
    N <- length(mixturePCA_fit$z)
    if (any(cluster_sizes < 0.05 * N)) {
      warnings <- c(warnings, sprintf("MPCA: Very small clusters detected (min size: %d / %d subjects)", 
                                     min(cluster_sizes), N))
    }
    
    # Check if log-likelihood is finite
    if (!is.finite(mixturePCA_fit$logLik)) {
      warnings <- c(warnings, "MPCA: Non-finite log-likelihood - numerical issues detected")
      convergence_ok <- FALSE
    }
  }
  
  # ============================================================================
  # 2. PSI (NOISE VARIANCE) DISTRIBUTION CHECK
  # ============================================================================
  
  psi_ok <- TRUE
  psi_summary_mfa <- NULL
  psi_summary_mpca <- NULL
  
  # Check MFA Psi distribution
  if (!is.null(mixtureFA_fit) && !is.null(mixtureFA_fit$Psi)) {
    # Extract all Psi diagonal values across all clusters
    K <- length(mixtureFA_fit$Psi)
    all_psi_values <- numeric(0)
    
    for (k in 1:K) {
      psi_k <- diag(mixtureFA_fit$Psi[[k]])
      all_psi_values <- c(all_psi_values, psi_k)
    }
    
    # Compute summary statistics
    psi_summary_mfa <- summary(all_psi_values)
    
    # Check for extreme values
    if (any(all_psi_values < psi_small_thresh)) {
      n_small <- sum(all_psi_values < psi_small_thresh)
      warnings <- c(warnings, sprintf("MFA: %d Psi values < %.3f detected (possible overfitting)", 
                                     n_small, psi_small_thresh))
      psi_ok <- FALSE
    }
    
    if (any(all_psi_values > psi_large_thresh)) {
      n_large <- sum(all_psi_values > psi_large_thresh)
      warnings <- c(warnings, sprintf("MFA: %d Psi values > %.1f detected (possible poor fit)", 
                                     n_large, psi_large_thresh))
      psi_ok <- FALSE
    }
    
    # Create histogram
    if (!is.null(output_dir)) {
      pdf(file.path(output_dir, "psi_distribution_mfa.pdf"), width = 8, height = 6)
    }
    
    hist(all_psi_values, breaks = 30, 
         main = "MFA: Distribution of Psi (Noise Variance)",
         xlab = "Psi values",
         col = "steelblue",
         border = "white")
    abline(v = psi_small_thresh, col = "red", lty = 2, lwd = 2)
    abline(v = psi_large_thresh, col = "red", lty = 2, lwd = 2)
    legend("topright", 
           legend = sprintf("Warning thresholds (%.3f, %.1f)", psi_small_thresh, psi_large_thresh), 
           col = "red", lty = 2, lwd = 2)
    
    if (!is.null(output_dir)) {
      dev.off()
    }
  }
  
  # Check MPCA Psi distribution (sigma2 values)
  if (!is.null(mixturePCA_fit) && !is.null(mixturePCA_fit$sigma2)) {
    # For MPCA, Psi = sigma2 * I, so we check sigma2 values
    all_sigma2_values <- mixturePCA_fit$sigma2
    
    # Compute summary statistics
    psi_summary_mpca <- summary(all_sigma2_values)
    
    # Check for extreme values
    if (any(all_sigma2_values < psi_small_thresh)) {
      n_small <- sum(all_sigma2_values < psi_small_thresh)
      warnings <- c(warnings, sprintf("MPCA: %d sigma2 values < %.3f detected (possible overfitting)", 
                                     n_small, psi_small_thresh))
      psi_ok <- FALSE
    }
    
    if (any(all_sigma2_values > psi_large_thresh)) {
      n_large <- sum(all_sigma2_values > psi_large_thresh)
      warnings <- c(warnings, sprintf("MPCA: %d sigma2 values > %.1f detected (possible poor fit)", 
                                     n_large, psi_large_thresh))
      psi_ok <- FALSE
    }
    
    # Create histogram
    if (!is.null(output_dir)) {
      pdf(file.path(output_dir, "psi_distribution_mpca.pdf"), width = 8, height = 6)
    }
    
    hist(all_sigma2_values, breaks = 10, 
         main = "MPCA: Distribution of sigma2 (Noise Variance)",
         xlab = "sigma2 values",
         col = "coral",
         border = "white")
    abline(v = psi_small_thresh, col = "red", lty = 2, lwd = 2)
    abline(v = psi_large_thresh, col = "red", lty = 2, lwd = 2)
    legend("topright", 
           legend = sprintf("Warning thresholds (%.3f, %.1f)", psi_small_thresh, psi_large_thresh), 
           col = "red", lty = 2, lwd = 2)
    
    if (!is.null(output_dir)) {
      dev.off()
    }
  }
  
  # ============================================================================
  # 3. VAF-SSE CORRELATION VERIFICATION
  # ============================================================================
  
  correlation_ok <- TRUE
  vaf_sse_correlation <- NA
  
  # Remove NA values (vaf_values and sse_values already extracted in input validation)
  valid_idx <- !is.na(vaf_values) & !is.na(sse_values)
  vaf_values_clean <- vaf_values[valid_idx]
  sse_values_clean <- sse_values[valid_idx]
  
  if (length(vaf_values_clean) > 1 && length(sse_values_clean) > 1) {
    # Compute correlation
    vaf_sse_correlation <- cor(vaf_values_clean, sse_values_clean, use = "complete.obs")
    
    # Check if correlation is strongly negative (should be approximately -1)
    if (is.finite(vaf_sse_correlation)) {
      if (vaf_sse_correlation > -0.9) {
        warnings <- c(warnings, sprintf(
          "VAF-SSE correlation = %.3f (expected approximately -1.0) - possible scaling inconsistency",
          vaf_sse_correlation
        ))
        correlation_ok <- FALSE
      }
      
      # Check if correlation is positive (definitely wrong)
      if (vaf_sse_correlation > 0) {
        warnings <- c(warnings, sprintf(
          "VAF-SSE correlation = %.3f is POSITIVE - serious scaling issue detected",
          vaf_sse_correlation
        ))
        correlation_ok <- FALSE
      }
    } else {
      warnings <- c(warnings, "VAF-SSE correlation is non-finite - numerical issues detected")
      correlation_ok <- FALSE
    }
    
    # Create scatter plot
    if (!is.null(output_dir)) {
      pdf(file.path(output_dir, "vaf_sse_correlation.pdf"), width = 8, height = 6)
    }
    
    plot(sse_values_clean, vaf_values_clean, 
         pch = 19, cex = 1.5, col = "steelblue",
         xlab = "SSE (Sum of Squared Errors)",
         ylab = "VAF (Variance Accounted For)",
         main = sprintf("VAF vs SSE (correlation = %.3f)", vaf_sse_correlation))
    
    # Add theoretical relationship line (VAF = 1 - SSE/SST)
    # Compute SST from the data
    if (!is.null(list_of_data)) {
      total_SST <- calc_total_SST(list_of_data)
      sse_range <- seq(min(sse_values_clean), max(sse_values_clean), length.out = 100)
      vaf_theoretical <- 1 - sse_range / total_SST
      lines(sse_range, vaf_theoretical, col = "red", lwd = 2, lty = 2)
      legend("topright", 
             legend = c("Observed", "Theoretical (VAF = 1 - SSE/SST)"),
             col = c("steelblue", "red"),
             pch = c(19, NA),
             lty = c(NA, 2),
             lwd = c(NA, 2))
    }
    
    grid()
    
    if (!is.null(output_dir)) {
      dev.off()
    }
  } else {
    warnings <- c(warnings, "Insufficient data points to compute VAF-SSE correlation")
    correlation_ok <- FALSE
  }
  
  # ============================================================================
  # SUMMARY
  # ============================================================================
  
  if (length(warnings) == 0) {
    message("All diagnostic checks passed")
  } else {
    message("Diagnostic warnings detected:")
    for (w in warnings) {
      message("  - ", w)
    }
  }
  
  # Return diagnostic results
  return(list(
    convergence_ok = convergence_ok,
    psi_ok = psi_ok,
    correlation_ok = correlation_ok,
    warnings = warnings,
    psi_summary_mfa = psi_summary_mfa,
    psi_summary_mpca = psi_summary_mpca,
    vaf_sse_correlation = vaf_sse_correlation
  ))
}
