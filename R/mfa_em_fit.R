#' Single-Run Mixture Factor Analysis with a Preset Initial Assignment
#'
#' This helper function calls the C++ routine \code{mfaTimeseriesCpp()} with a user-specified
#' cluster assignment (\code{z_init}) and then computes the final log-likelihood and
#' responsibilities in R.
#'
#' @param list_of_data A list of \code{(T_i x M)} matrices, one for each subject.
#' @param K Number of clusters.
#' @param r Factor dimension in each cluster.
#' @param z_init An integer vector of length \code{N} (where \code{N=length(list_of_data)}),
#'   giving the initial cluster label (1..K) for each subject.
#' @param max_iter Maximum EM iterations in C++.
#' @param nIterFA Number of sub-iterations for the factor analyzer update in C++.
#' @param tol Convergence tolerance used in C++.
#' @param verbose Logical; if \code{TRUE}, print progress messages. Default is \code{TRUE}.
#' @param n_threads Integer; number of OpenMP threads to use (0 = auto). Default is \code{0}.
#' @param seed Integer; random seed for C++ RNG (0 = default/no seed). Default is \code{0}.
#' @param progress_callback Optional function to report progress from C++. Should accept two arguments: current iteration and total iterations.
#' @param em_type Character; type of EM algorithm to use. Either "hard" (classification EM with
#'   hard assignments) or "soft" (standard EM with responsibility-weighted updates). Default is "hard".
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{z}}{Hard cluster assignments, length \code{N}.}
#'   \item{\code{pi}}{Cluster mixing proportions, length \code{K}.}
#'   \item{\code{mu}}{List of length \code{K}, each a mean vector (\code{M x 1}).}
#'   \item{\code{Lambda}}{List of length \code{K}, each an \code{M x r} factor loading matrix.}
#'   \item{\code{Psi}}{List of length \code{K}, each an \code{M x M} diagonal noise matrix.}
#'   \item{\code{logLik}}{Final log-likelihood computed in R.}
#'   \item{\code{resp}}{An \code{N x K} matrix of responsibilities.}
#' }
#'
#' @examples
#' \dontrun{
#' # Suppose we have list_of_data, K=3, r=2:
#' z_init <- sample.int(3, size=length(list_of_data), replace=TRUE)
#' fit_one <- mfa_em_fit_cpp_singleInit(list_of_data, K=3, r=2,
#'                                      z_init=z_init)
#' fit_one$logLik
#' }
#'
#' @export
mfa_em_fit_cpp_singleInit <- function(list_of_data,
                                      K,
                                      r,
                                      z_init,
                                      max_iter = 50,
                                      nIterFA  = 20,
                                      tol      = 1e-3,
                                      verbose  = TRUE,
                                      n_threads = 0,
                                      seed = 0,
                                      progress_callback = NULL,
                                      em_type = c("hard", "soft"))
{
  em_type <- match.arg(em_type, c("hard", "soft"))

  # 1) Call the appropriate C++ function based on em_type
  if(em_type == "soft"){
    fit_cpp <- mfaTimeseriesSoftCpp(
      list_of_data = list_of_data,
      K = K,
      r = r,
      max_iter = max_iter,
      nIterFA  = nIterFA,
      tol      = tol,
      z_init   = z_init,
      verbose  = verbose,
      n_threads = n_threads,
      seed     = seed
    )
  } else {
    fit_cpp <- mfaTimeseriesCpp(
      list_of_data = list_of_data,
      K = K,
      r = r,
      max_iter = max_iter,
      nIterFA  = nIterFA,
      tol      = tol,
      z_init   = z_init,
      verbose  = verbose,
      n_threads = n_threads,
      seed     = seed,
      progress_callback = progress_callback
    )
  }

  # 2) Compute final log-likelihood & responsibilities in R
  N <- length(list_of_data)
  Kc <- length(fit_cpp$Lambda)
  if(Kc != K){
    stop("C++ output mismatch: length(Lambda) != K.")
  }

  # Build Sigma_k
  Sigma_list <- vector("list", Kc)
  for(k2 in seq_len(Kc)){
    Lambda_k <- fit_cpp$Lambda[[k2]]
    Psi_k    <- fit_cpp$Psi[[k2]]
    Sig_k    <- Lambda_k %*% t(Lambda_k) + Psi_k
    Sigma_list[[k2]] <- Sig_k
  }
  pi_vec <- fit_cpp$pi

  # Accumulate log-likelihood
  logLik_val <- 0
  resp <- matrix(0, nrow=N, ncol=Kc)

  for(i in seq_len(N)){
    Xi <- list_of_data[[i]]
    T_i <- nrow(Xi)
    logvals <- numeric(Kc)

    for(k2 in seq_len(Kc)){
      mu_k  <- fit_cpp$mu[[k2]]
      Sig_k <- Sigma_list[[k2]]
      dens_t <- mvtnorm::dmvnorm(Xi, mean=mu_k, sigma=Sig_k, log=TRUE)
      sumLog <- sum(dens_t)
      logvals[k2] <- log(pi_vec[k2] + 1e-16) + sumLog
    }

    # log-sum-exp
    m0 <- max(logvals)
    li <- m0 + log(sum(exp(logvals - m0)))
    logLik_val <- logLik_val + li

    # responsibilities
    for(k2 in seq_len(Kc)){
      resp[i,k2] <- exp(logvals[k2] - li)
    }
  }

  # For soft EM, use gamma from C++ as responsibilities if available
  if(em_type == "soft" && !is.null(fit_cpp$gamma)){
    resp <- fit_cpp$gamma
  }

  # Final output
  out <- list(
    z      = fit_cpp$z,
    pi     = fit_cpp$pi,
    mu     = fit_cpp$mu,
    Lambda = fit_cpp$Lambda,
    Psi    = fit_cpp$Psi,
    logLik = logLik_val,
    resp   = resp
  )
  out
}


#' Fit a Mixture Factor Analysis Model (with optional multi-initialization)
#'
#' This function calls the C++ function \code{mfaTimeseriesCpp()} to perform a Mixture
#' Factor Analysis EM algorithm for fixed \code{K} and \code{r}. By default, it runs
#' a single pass with the internal initialization from C++. However, if you specify
#' multiple initial attempts (via \code{n_init>1}) and/or \code{use_kmeans_init=TRUE},
#' this function will try several different initial cluster assignments (in parallel
#' using \code{mclapply}), then return the best solution (maximizing the final log-likelihood).
#'
#' @param list_of_data A list of matrices (each \code{(T_i x M)}) to be modeled.
#' @param K Number of clusters.
#' @param r Factor dimension.
#' @param max_iter Maximum EM iterations for \code{mfaTimeseriesCpp}.
#' @param iter_fa Number of sub-iterations for the factor analyzer update in C++.
#'   (Replaces deprecated \code{nIterFA})
#' @param tol Convergence tolerance for \code{mfaTimeseriesCpp}.
#' @param n_init Number of random initial assignments to try (in addition to the default
#'   single-run or k-means if requested). Defaults to \code{1}.
#' @param use_kmeans_init If \code{TRUE}, we also run one initialization where we
#'   extract subject features (via PCA) and apply \code{kmeans} to get a cluster assignment.
#'   Defaults to \code{FALSE}.
#' @param kmeans_rdim PCA dimension for extracting subject-level features,
#'   used only if \code{use_kmeans_init=TRUE}. Default is \code{r}.
#'   (Replaces deprecated \code{subject_rdim_for_kmeans})
#' @param mc_cores Number of cores for parallel execution via \code{mclapply}.
#'   Defaults to \code{1} (no parallel).
#' @param seed Random seed for reproducibility (default: NULL).
#' @param verbose Logical; if \code{TRUE}, print progress messages during EM iterations. Default is \code{TRUE}.
#' @param n_threads Integer; number of OpenMP threads to use for parallel computation (0 = auto). Default is \code{0}.
#' @param cpp_seed Integer; random seed for C++ RNG (0 = default/no seed). Default is \code{0}.
#' @param progress_callback Optional function to report progress from C++. Should accept two arguments: current iteration and total iterations.
#' @param ... Additional arguments for backwards compatibility. Deprecated arguments like
#'   \code{nIterFA} and \code{subject_rdim_for_kmeans} are accepted with warnings.
#'
#' @return A list with the same elements as a single run: \code{z, pi, mu, Lambda, Psi, logLik, resp}.
#'   It corresponds to the best solution (i.e. highest log-likelihood) among all tried inits.
#'
#' @details
#' The single-run approach in C++ already does a built-in initialization if \code{z_init}
#' is not provided. Therefore, if \code{n_init=1} and \code{use_kmeans_init=FALSE}, we just call
#' \code{mfaTimeseriesCpp} once (like the original design).
#' Otherwise:
#' \enumerate{
#'   \item If \code{use_kmeans_init=TRUE}, we do one run where we assign clusters by k-means
#'         on some subject-level features (extracted by PCA).
#'   \item We also sample \code{n_init} random assignments (each subject assigned to a random cluster).
#'   \item For each assignment, we call \code{\link{mfa_em_fit_cpp_singleInit}}.
#'   \item We compare all results' \code{logLik} and keep the best one.
#' }
#'
#' This provides a more robust initialization scheme to avoid poor local maxima.
#'
#' @examples
#' \dontrun{
#' # Suppose we have a list_of_data for N=100 subjects, each (T_i x M=8).
#' # We want to fit K=3, r=2, and try 5 random inits plus one k-means init:
#' fit <- mfa_em_fit(
#'   list_of_data, K=3, r=2,
#'   max_iter=50, iter_fa=20, tol=1e-3,
#'   n_init=5,
#'   use_kmeans_init=TRUE,
#'   kmeans_rdim=2,
#'   mc_cores=2,
#'   seed=123
#' )
#' print(fit$logLik)
#' head(fit$z)
#' }
#'
#' @seealso \code{\link{mixture_pca_em_fit}} for MPCA fitting,
#'   \code{\link{select_r}} for selecting optimal r,
#'   \code{\link{plot_cluster_synergy_loadings_mfa}} for visualization
#'
#' @export
mfa_em_fit <- function(
    list_of_data,
    K,
    r,
    max_iter = 50,
    iter_fa  = 20,
    tol      = 1e-3,
    n_init   = 1,
    use_kmeans_init = TRUE,
    kmeans_rdim = r,
    mc_cores = 1,
    seed = NULL,
    verbose = TRUE,
    n_threads = 0,
    cpp_seed = 0,
    progress_callback = NULL,
    ...
)
{
  # Handle deprecated argument names
  dots <- list(...)
  if ("nIterFA" %in% names(dots)) {
    warning("Argument 'nIterFA' is deprecated. Use 'iter_fa' instead.", call. = FALSE)
    if (missing(iter_fa)) iter_fa <- dots$nIterFA
  }
  if ("subject_rdim_for_kmeans" %in% names(dots)) {
    warning("Argument 'subject_rdim_for_kmeans' is deprecated. Use 'kmeans_rdim' instead.", call. = FALSE)
    if (missing(kmeans_rdim)) kmeans_rdim <- dots$subject_rdim_for_kmeans
  }

  # Input validation using standardized validation functions
  validate_list_of_data(list_of_data)
  validate_model_params(K = K, r = r, max_iter = max_iter, n_init = n_init, tol = tol)

  N <- length(list_of_data)
  M <- ncol(list_of_data[[1]])

  validate_K_vs_N(K, N)
  validate_r_vs_M(r, M)

  # Additional parameter validation
  stopifnot(
    "iter_fa must be positive" = is.numeric(iter_fa) && iter_fa > 0,
    "use_kmeans_init must be logical" = is.logical(use_kmeans_init),
    "mc_cores must be a positive integer" = is.numeric(mc_cores) && mc_cores > 0 && mc_cores == floor(mc_cores)
  )
  
  # Set seed for reproducibility if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Windows parallel fallback
  if (.Platform$OS.type == "windows" && mc_cores > 1) {
    warning("Windows detected: forcing sequential execution (mc_cores = 1)")
    mc_cores <- 1
  }
  
  # Set n_threads to 1 if not specified (default behavior for thread control)
  if (n_threads == 0) {
    n_threads <- 1
  }
  
  # Wrap the entire fitting process with thread control to prevent oversubscription
  # when this function is called from parallel::mclapply() at a higher level
  result <- with_internal_threads(threads = n_threads, {
  
  # If n_init=1 and use_kmeans_init=FALSE => single-run original approach
  # i.e., just call the C++ code once, relying on its built-in (i%K) init
  if(n_init == 1 && !use_kmeans_init){
    # (original style) => call mfaTimeseriesCpp with no z_init
    fit_cpp <- mfaTimeseriesCpp(
      list_of_data = list_of_data,
      K = K,
      r = r,
      max_iter = max_iter,
      nIterFA  = iter_fa,
      tol      = tol,
      # no z_init => it uses the default (i%K)+1
      verbose  = verbose,
      n_threads = n_threads,
      seed     = cpp_seed,
      progress_callback = progress_callback
    )
    # Now compute final logLik & resp in R (similar to your original post-fitting code)
    N <- length(list_of_data)
    Sigma_list <- vector("list", K)
    for(k2 in seq_len(K)){
      Lambda_k <- fit_cpp$Lambda[[k2]]
      Psi_k    <- fit_cpp$Psi[[k2]]
      Sigma_list[[k2]] <- Lambda_k %*% t(Lambda_k) + Psi_k
    }
    pi_vec <- fit_cpp$pi
    logLik_val <- 0
    resp <- matrix(0, nrow=N, ncol=K)

    for(i in seq_len(N)){
      Xi <- list_of_data[[i]]
      logvals <- numeric(K)
      for(k2 in seq_len(K)){
        mu_k  <- fit_cpp$mu[[k2]]
        sig_k <- Sigma_list[[k2]]
        dens_t <- mvtnorm::dmvnorm(Xi, mean=mu_k, sigma=sig_k, log=TRUE)
        sumLog <- sum(dens_t)
        logvals[k2] <- log(pi_vec[k2] + 1e-16) + sumLog
      }
      # log-sum-exp
      m0 <- max(logvals)
      li <- m0 + log(sum(exp(logvals - m0)))
      logLik_val <- logLik_val + li

      # responsibilities
      for(k2 in seq_len(K)){
        resp[i,k2] <- exp(logvals[k2] - li)
      }
    }
    out <- list(
      z      = fit_cpp$z,
      pi     = fit_cpp$pi,
      mu     = fit_cpp$mu,
      Lambda = fit_cpp$Lambda,
      Psi    = fit_cpp$Psi,
      logLik = logLik_val,
      resp   = resp,
      thread_info = list(
        n_threads = n_threads,
        mc_cores = mc_cores,
        n_init = n_init,
        use_kmeans_init = use_kmeans_init,
        cpp_seed = cpp_seed,
        used_parallel = FALSE
      )
    )
    return(out)
  }

  # Otherwise => multi-init approach

  # We will store results in a list and pick the best by logLik
  best_fit    <- NULL
  best_logLik <- -Inf

  # 1) Possibly do a k-means-based init
  if(use_kmeans_init){
    # If user wants k-means, we do it once
    features <- extract_subject_features_by_singlePCA(
      list_of_data, r_dim=kmeans_rdim
    )
    z_init_km <- assign_by_kmeans(features, K=K)

    fit_km <- mfa_em_fit_cpp_singleInit(
      list_of_data = list_of_data,
      K = K,
      r = r,
      z_init = z_init_km,
      max_iter = max_iter,
      nIterFA  = iter_fa,
      tol      = tol,
      verbose  = verbose,
      n_threads = n_threads,
      seed     = cpp_seed,
      progress_callback = progress_callback
    )
    if(fit_km$logLik > best_logLik){
      best_fit    <- fit_km
      best_logLik <- fit_km$logLik
    }
  }

  # 2) Multiple random initial assignments => parallel
  if(!requireNamespace("parallel", quietly=TRUE)){
    stop("Package 'parallel' is required for mclapply. Please install it.")
  }

  N <- length(list_of_data)
  random_inits <- lapply(seq_len(n_init), function(x){
    sample.int(K, size=N, replace=TRUE)
  })

  fit_list <- parallel::mclapply(
    random_inits,
    FUN = function(z_init_rand){
      # single-run with that random init
      fit_rand <- mfa_em_fit_cpp_singleInit(
        list_of_data = list_of_data,
        K = K,
        r = r,
        z_init = z_init_rand,
        max_iter = max_iter,
        nIterFA  = iter_fa,
        tol      = tol,
        verbose  = verbose,
        n_threads = n_threads,
        seed     = cpp_seed,
        progress_callback = progress_callback
      )
      fit_rand
    },
    mc.cores = mc_cores
  )

  # Check all random results
  for(fit_rand in fit_list){
    if(fit_rand$logLik > best_logLik){
      best_fit    <- fit_rand
      best_logLik <- fit_rand$logLik
    }
  }

  # Add thread_info to the best fit result
  best_fit$thread_info <- list(
    n_threads = n_threads,
    mc_cores = mc_cores,
    n_init = n_init,
    use_kmeans_init = use_kmeans_init,
    cpp_seed = cpp_seed,
    used_parallel = (mc_cores > 1)
  )
  
  # Return the best
  best_fit
  
  })  # End of with_internal_threads()
  
  return(result)
}
