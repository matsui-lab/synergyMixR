# Roxygen Parameter Documentation Templates
#
# This file contains reusable parameter documentation that can be
# inherited by other functions using @inheritParams.

#' Package-level imports for base R functions
#'
#' This dummy function provides package-level imports for base R functions
#' used throughout the package. These imports are needed to pass R CMD check.
#'
#' @importFrom grDevices colorRampPalette dev.off pdf rainbow
#' @importFrom graphics abline axis barplot grid hist image legend lines par points text
#' @importFrom stats as.formula median quantile time var
#' @importFrom utils flush.console write.csv modifyList
#' @keywords internal
#' @noRd
NULL

#' Parameter Documentation Templates
#'
#' Standard parameter documentation for synergyMixR functions.
#' Use `@inheritParams param-templates` to inherit these definitions.
#'
#' @param list_of_data A list of N matrices, each of dimension (T_i x M),
#'   where T_i is the number of time points for subject i and M is the
#'   number of muscles/variables.
#' @param K Integer. Number of clusters (mixture components). Must be >= 1.
#' @param r Integer. Number of synergies (latent factors). Must be >= 1.
#' @param max_iter Integer. Maximum number of EM iterations. Default: 100.
#' @param iter_fa Integer. Number of sub-iterations for factor analyzer update.
#'   Default: 20.
#' @param iter_pca Integer. Number of sub-iterations for PCA/PPCA update.
#'   Default: 20.
#' @param tol Numeric. Convergence tolerance for EM algorithm. Default: 1e-3.
#' @param n_init Integer. Number of random initializations. Default: 1.
#' @param use_kmeans_init Logical. If TRUE, also run one initialization using
#'   k-means clustering on subject-level features. Default: TRUE.
#' @param kmeans_rdim Integer. PCA dimension for extracting subject-level
#'   features when use_kmeans_init is TRUE. Default: r.
#' @param mc_cores Integer. Number of CPU cores for parallel processing.
#'   Default: 1 (sequential).
#' @param n_threads Integer. Number of OpenMP threads for C++ parallel
#'   computation. 0 = auto-detect. Default: 0.
#' @param model_type Character. Either "MFA" (Mixture Factor Analysis) or
#'   "MPCA" (Mixture PCA).
#' @param verbose Logical. If TRUE, print progress messages. Default: TRUE.
#' @param seed Integer. Random seed for R-level reproducibility. Default: NULL.
#' @param cpp_seed Integer. Random seed for C++ RNG. 0 = no seed. Default: 0.
#' @param progress_callback Function. Optional callback for progress reporting.
#'   Should accept two arguments: current iteration and total iterations.
#' @param fit A fitted model object from mfa_em_fit() or mixture_pca_em_fit().
#' @param rvec Integer vector. Candidate r values to evaluate in model selection.
#'   Default: 1:6.
#' @param method Character. Model selection method (varies by function).
#' @param scope Character. Scope of selection: "global", "clusterwise", or "hybrid".
#'
#' @name param-templates
#' @keywords internal
NULL


#' Return Value Templates
#'
#' Standard return value documentation for synergyMixR functions.
#'
#' @name return-templates
#' @keywords internal
#'
#' @section MFA Fit Return Value:
#' For MFA models (mfa_em_fit), returns a list with:
#' \describe{
#'   \item{z}{Integer vector of length N. Hard cluster assignments (1 to K).}
#'   \item{pi}{Numeric vector of length K. Cluster mixing proportions.}
#'   \item{mu}{List of length K. Each element is an M-dimensional mean vector.}
#'   \item{Lambda}{List of length K. Each element is an M x r factor loading matrix.}
#'   \item{Psi}{List of length K. Each element is an M x M diagonal noise covariance.}
#'   \item{logLik}{Numeric. Final log-likelihood value.}
#'   \item{resp}{N x K matrix of cluster responsibilities (soft assignments).}
#' }
#'
#' @section MPCA Fit Return Value:
#' For MPCA models (mixture_pca_em_fit), returns a list with:
#' \describe{
#'   \item{z}{Integer vector of length N. Hard cluster assignments (1 to K).}
#'   \item{pi}{Numeric vector of length K. Cluster mixing proportions.}
#'   \item{mu}{List of length K. Each element is an M-dimensional mean vector.}
#'   \item{W}{List of length K. Each element is an M x r projection matrix.}
#'   \item{P}{List of length K. Each element is an M x r principal directions matrix.}
#'   \item{D}{List of length K. Each element is an r x r diagonal scaling matrix.}
#'   \item{Psi}{List of length K. Each element is an M x M diagonal noise covariance.}
#'   \item{sigma2}{Numeric vector of length K. Noise variances per cluster.}
#'   \item{logLik}{Numeric. Final log-likelihood value.}
#'   \item{resp}{N x K matrix of cluster responsibilities (soft assignments).}
#' }
NULL
