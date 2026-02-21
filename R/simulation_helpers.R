#' Internal Helper Functions for Simulation
#'
#' This file contains internal utility functions used by the simulation subsystem.
#' These functions are not exported and are intended for internal package use only.
#'
#' @name simulation_helpers
#' @keywords internal
NULL

# Internal BIC computation helper
compute_bic <- function(logLik, K, r, M, N_total_rows) {
  num_params <- K * (M * r + 2 * M) + (K - 1)
  -2 * logLik + num_params * log(N_total_rows)
}

# Internal helper for reconstruction (not exported)
# Different implementation from visualization_mfa.R version
reconstruct_singleFA_sim <- function(X, mu, Lambda, Psi_diag) {
  Xc <- sweep(X, 2, mu, "-")
  M <- ncol(X)
  r <- ncol(Lambda)

  Sigma <- Lambda %*% t(Lambda)
  diag(Sigma) <- diag(Sigma) + Psi_diag
  invS <- solve(Sigma)

  W <- t(Lambda) %*% invS  # (r x M)
  A <- diag(r) + W %*% Lambda
  Ainv <- solve(A)
  EZt <- Ainv %*% W %*% t(Xc)
  EZ  <- t(EZt)
  
  # Reconstruct: Xhat = mu + Lambda * EZ^T
  Xhat <- sweep(EZ %*% t(Lambda), 2, mu, "+")
  Xhat
}

#' Compute Factor Scores for Single FA Model
#'
#' Internal helper to compute factor scores for a single cluster's FA parameters.
#'
#' @param X Numeric matrix (T x M) of data for one subject.
#' @param mu Mean vector (length M) for the cluster.
#' @param Lambda Loading matrix (M x r) for the cluster.
#' @param diagPsi Diagonal of noise covariance (length M).
#' @return A matrix of factor scores (T x r).
#' @keywords internal
compute_factor_scores_singleFA <- function(X, mu, Lambda, diagPsi) {
  Xc <- sweep(X, 2, mu, "-")
  M <- ncol(X)
  r <- ncol(Lambda)

  Sigma <- Lambda %*% t(Lambda)
  diag(Sigma) <- diag(Sigma) + diagPsi
  invS <- solve(Sigma)

  W <- t(Lambda) %*% invS  # (r x M)
  A <- diag(r) + W %*% Lambda
  Ainv <- solve(A)
  EZt <- Ainv %*% W %*% t(Xc)
  EZ  <- t(EZt)
  EZ
}

#' Reconstruct Data Using Single PCA Model
#'
#' Internal helper to reconstruct data matrix using PCA parameters.
#'
#' @param X Numeric matrix (T x M) of data for one subject.
#' @param pca_res PCA result containing P (loadings) and mu (mean).
#' @return A matrix of reconstructed data (T x M).
#' @keywords internal
reconstruct_singlePCA <- function(X, pca_res) {
  P  <- pca_res$P
  mu <- pca_res$mu
  Xc <- sweep(X, 2, mu, "-")
  Z <- Xc %*% P
  Xhat <- Z %*% t(P)
  Xhat <- sweep(Xhat, 2, mu, "+")
  Xhat
}

#' Reconstruct Data for One Cluster in MixtureFA Model
#'
#' Internal helper to reconstruct data for a single subject using one cluster's MFA parameters.
#'
#' @param X Numeric matrix (T x M) of data for one subject.
#' @param k Cluster index.
#' @param mfa_fit MFA model fit containing cluster parameters.
#' @return A matrix of reconstructed data (T x M).
#' @keywords internal
reconstruct_mixtureFA_oneCluster <- function(X, k, mfa_fit) {
  mu_k      <- mfa_fit$mu[[k]]
  Lambda_k  <- mfa_fit$Lambda[[k]]
  diagPsi_k <- diag(mfa_fit$Psi[[k]])
  Xhat <- reconstruct_singleFA(X, mu_k, Lambda_k, diagPsi_k)
  Xhat
}

#' Reconstruct Data for One Cluster in MixturePCA Model
#'
#' Internal helper to reconstruct data for a single subject using one cluster's MPCA parameters.
#'
#' @param X Numeric matrix (T x M) of data for one subject.
#' @param k Cluster index.
#' @param mpca_fit MPCA model fit containing cluster parameters.
#' @return A matrix of reconstructed data (T x M).
#' @keywords internal
reconstruct_mixturePCA_oneCluster <- function(X, k, mpca_fit) {
  P_k  <- mpca_fit$P[[k]]
  mu_k <- mpca_fit$mu[[k]]
  # Simple PCA reconstruction
  Xc <- sweep(X, 2, mu_k, "-")
  Z  <- Xc %*% P_k
  Xhat <- Z %*% t(P_k)
  Xhat <- sweep(Xhat, 2, mu_k, "+")
  Xhat
}

#' Execute Code with Controlled Internal Threading
#'
#' This helper function sets environment variables to control the number of threads
#' used by BLAS/LAPACK and OpenMP operations, executes the provided expression,
#' and then restores the original environment variables.
#'
#' @param threads Integer number of threads to use (default: 1). Must be at least 1.
#' @param expr Expression to evaluate with controlled threading.
#' @return The result of evaluating \code{expr}.
#'
#' @details
#' This function sets the following environment variables to control threading:
#' \itemize{
#'   \item \code{OMP_NUM_THREADS} - OpenMP thread count
#'   \item \code{OPENBLAS_NUM_THREADS} - OpenBLAS thread count
#'   \item \code{MKL_NUM_THREADS} - Intel MKL thread count
#'   \item \code{BLIS_NUM_THREADS} - BLIS thread count
#'   \item \code{OMP_DYNAMIC} - Disable dynamic thread adjustment
#'   \item \code{MKL_DYNAMIC} - Disable MKL dynamic thread adjustment
#' }
#'
#' The original values are restored via \code{on.exit()}, ensuring cleanup
#' even if \code{expr} throws an error.
#'
#' This is particularly useful when using \code{parallel::mclapply()} or similar
#' parallelization at a higher level, to prevent CPU oversubscription from nested
#' parallelism.
#'
#' @keywords internal
with_internal_threads <- function(threads = 1, expr) {
  # Validate threads parameter
  threads <- as.integer(threads)
  if (threads < 1) {
    warning("threads must be at least 1; setting to 1")
    threads <- 1L
  }
  
  threads_str <- as.character(threads)
  
  # Save current environment variables
  old_omp <- Sys.getenv("OMP_NUM_THREADS", unset = NA)
  old_openblas <- Sys.getenv("OPENBLAS_NUM_THREADS", unset = NA)
  old_mkl <- Sys.getenv("MKL_NUM_THREADS", unset = NA)
  old_blis <- Sys.getenv("BLIS_NUM_THREADS", unset = NA)
  old_omp_dynamic <- Sys.getenv("OMP_DYNAMIC", unset = NA)
  old_mkl_dynamic <- Sys.getenv("MKL_DYNAMIC", unset = NA)
  
  # Set thread control environment variables
  Sys.setenv(
    OMP_NUM_THREADS = threads_str,
    OPENBLAS_NUM_THREADS = threads_str,
    MKL_NUM_THREADS = threads_str,
    BLIS_NUM_THREADS = threads_str,
    OMP_DYNAMIC = "FALSE",
    MKL_DYNAMIC = "FALSE"
  )
  
  # Restore original values on exit
  on.exit({
    if (is.na(old_omp)) {
      Sys.unsetenv("OMP_NUM_THREADS")
    } else {
      Sys.setenv(OMP_NUM_THREADS = old_omp)
    }
    
    if (is.na(old_openblas)) {
      Sys.unsetenv("OPENBLAS_NUM_THREADS")
    } else {
      Sys.setenv(OPENBLAS_NUM_THREADS = old_openblas)
    }
    
    if (is.na(old_mkl)) {
      Sys.unsetenv("MKL_NUM_THREADS")
    } else {
      Sys.setenv(MKL_NUM_THREADS = old_mkl)
    }
    
    if (is.na(old_blis)) {
      Sys.unsetenv("BLIS_NUM_THREADS")
    } else {
      Sys.setenv(BLIS_NUM_THREADS = old_blis)
    }
    
    if (is.na(old_omp_dynamic)) {
      Sys.unsetenv("OMP_DYNAMIC")
    } else {
      Sys.setenv(OMP_DYNAMIC = old_omp_dynamic)
    }
    
    if (is.na(old_mkl_dynamic)) {
      Sys.unsetenv("MKL_DYNAMIC")
    } else {
      Sys.setenv(MKL_DYNAMIC = old_mkl_dynamic)
    }
  }, add = TRUE)
  
  # Capture current thread settings for debugging/testing
  thread_info <- list(
    omp = Sys.getenv("OMP_NUM_THREADS"),
    openblas = Sys.getenv("OPENBLAS_NUM_THREADS"),
    mkl = Sys.getenv("MKL_NUM_THREADS"),
    blis = Sys.getenv("BLIS_NUM_THREADS")
  )
  
  # Evaluate the expression
  result <- force(expr)
  
  # Attach thread info to result if it's a list
  if (is.list(result)) {
    result$thread_info <- thread_info
  }
  
  result
}
