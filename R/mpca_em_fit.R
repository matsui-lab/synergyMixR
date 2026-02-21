#' Single-Run Mixture PCA with a Preset Initial Assignment
#'
#' This helper function calls the updated C++ routine \code{mpcaTimeseriesCpp()}
#' with a user-specified cluster assignment (\code{z_init}), then converts
#' the returned \code{W} into \code{(P, D)} and constructs \code{Psi}, and finally
#' computes the final log-likelihood and responsibilities in R.
#'
#' @param list_of_data A list of matrices (each \code{T_i x M}).
#' @param K Number of clusters.
#' @param r Number of principal components (dimension).
#' @param max_iter Maximum EM iterations (C++).
#' @param nIterPCA Sub-iterations for updating each cluster's PPCA parameters in C++.
#' @param tol Convergence tolerance (C++).
#' @param method Either \code{"EM"} or \code{"closed_form"}, passed to \code{mpcaTimeseriesCpp}.
#' @param z_init An integer vector of length \code{N} giving initial cluster labels (1..K).
#' @param verbose Logical; if \code{TRUE}, print progress messages. Default is \code{TRUE}.
#' @param n_threads Integer; number of OpenMP threads to use (0 = auto). Default is \code{0}.
#' @param seed Integer; random seed for C++ RNG (0 = default/no seed). Default is \code{0}.
#' @param progress_callback Optional function to report progress from C++. Should accept two arguments: current iteration and total iterations.
#' @param em_type Character; type of EM algorithm to use. Either "hard" (classification EM with
#'   hard assignments) or "soft" (standard EM with responsibility-weighted updates). Default is "hard".
#'
#' @return A list with elements:
#'   \item{z}{Hard cluster assignments (1..K).}
#'   \item{pi}{Mixing proportions (length K).}
#'   \item{mu}{List of length K, each a mean vector (length M).}
#'   \item{W}{List of length K, each an \code{(M x r)} loading matrix (optional).}
#'   \item{P}{List of length K, each an \code{(M x r)} matrix of principal directions.}
#'   \item{D}{List of length K, each an \code{(r x r)} diagonal matrix of column norms.}
#'   \item{Psi}{List of length K, each an \code{(M x M)} diagonal matrix \code{sigma2 * I}.}
#'   \item{logLik}{Final log-likelihood computed in R.}
#'   \item{sigma2}{Numeric vector of length K, the final \code{sigma2} for each cluster.}
#'   \item{resp}{An \code{(N x K)} matrix of responsibilities.}
#'
#' @examples
#' \dontrun{
#' # Suppose we have list_of_data, K=3, r=2, method="EM", plus some z_init:
#' z_init <- sample.int(3, size=length(list_of_data), replace=TRUE)
#' fit_one <- mixture_pca_em_fit_cpp_singleInit(
#'   list_of_data, K=3, r=2, max_iter=50, nIterPCA=20, tol=1e-3,
#'   method="EM", z_init=z_init
#' )
#' print(fit_one$logLik)
#' # Now we can call compute_factor_scores_mpca(fit_one, list_of_data) because fit_one$W exists
#' }
#'
#' @export
mixture_pca_em_fit_cpp_singleInit <- function(list_of_data,
                                              K,
                                              r,
                                              max_iter = 50,
                                              nIterPCA = 20,
                                              tol      = 1e-3,
                                              method   = "EM",
                                              z_init,
                                              verbose  = TRUE,
                                              n_threads = 0,
                                              seed = 0,
                                              progress_callback = NULL,
                                              em_type = c("hard", "soft"))
{
  em_type <- match.arg(em_type, c("hard", "soft"))

  # 1) Call the appropriate C++ function based on em_type
  if(em_type == "soft"){
    fit_cpp <- mpcaTimeseriesSoftCpp(
      list_of_data = list_of_data,
      K = K,
      r = r,
      max_iter  = max_iter,
      nIterPCA  = nIterPCA,
      tol       = tol,
      method    = method,
      z_init    = z_init,
      verbose   = verbose,
      n_threads = n_threads,
      seed      = seed
    )
  } else {
    fit_cpp <- mpcaTimeseriesCpp(
      list_of_data = list_of_data,
      K = K,
      r = r,
      max_iter  = max_iter,
      nIterPCA  = nIterPCA,
      tol       = tol,
      method    = method,
      z_init    = z_init,
      verbose   = verbose,
      n_threads = n_threads,
      seed      = seed,
      progress_callback = progress_callback
    )
  }

  # 2) Extract results from C++
  #    W_list_cpp, mu_list_cpp, sigma2_vec, pi_vec, z_cpp
  z_cpp      <- fit_cpp$z
  W_list_cpp <- fit_cpp$W       # (K elements, each (M x r))
  mu_list_cpp<- fit_cpp$mu      # (K elements)
  sigma2_vec <- fit_cpp$sigma2  # length K
  pi_vec     <- fit_cpp$pi      # length K

  # 3) Convert W -> (P, D) and build Psi
  K_check <- length(W_list_cpp)
  if(K_check != K){
    stop("C++ output mismatch: length(W_list_cpp) != K.")
  }

  P_list   <- vector("list", K)
  D_list   <- vector("list", K)
  Psi_list <- vector("list", K)

  for(k2 in seq_len(K)){
    W_k <- W_list_cpp[[k2]]    # (M x r)
    M_k <- nrow(W_k)
    r_k <- ncol(W_k)

    # column norms -> diagonal of D
    col_scales <- numeric(r_k)
    for(j in seq_len(r_k)){
      cs_j <- sqrt(sum(W_k[,j]^2))
      if(cs_j < 1e-12) cs_j <- 1e-12
      col_scales[j] <- cs_j
    }
    # P_k = W_k with normalized columns
    P_k <- W_k
    for(j in seq_len(r_k)){
      P_k[, j] <- W_k[, j] / col_scales[j]
    }
    D_k <- diag(col_scales, r_k, r_k)

    sig2_k <- sigma2_vec[k2]
    Psi_k  <- diag(sig2_k, M_k)

    P_list[[k2]]  <- P_k
    D_list[[k2]]  <- D_k
    Psi_list[[k2]]<- Psi_k
  }

  # 4) Compute final log-likelihood and responsibilities in R
  N <- length(list_of_data)
  Sigma_list <- vector("list", K)
  for(k2 in seq_len(K)){
    # Sigma_k = W_k W_k^T + sigma2 * I
    # or equivalently P_k D_k^2 P_k^T + sigma2 * I
    W_k  <- W_list_cpp[[k2]]  # or use P_k * D_k
    sig2 <- sigma2_vec[k2]
    Sig_k <- W_k %*% t(W_k)
    diag(Sig_k) <- diag(Sig_k) + sig2
    Sigma_list[[k2]] <- Sig_k
  }

  logLik_val <- 0
  resp <- matrix(0, nrow=N, ncol=K)
  for(i in seq_len(N)){
    Xi <- list_of_data[[i]]
    T_i <- nrow(Xi)
    logvals <- numeric(K)
    for(k2 in seq_len(K)){
      mu_k  <- mu_list_cpp[[k2]]
      Sig_k <- Sigma_list[[k2]]
      dens_t <- mvtnorm::dmvnorm(Xi, mean=mu_k, sigma=Sig_k, log=TRUE)
      sumLog <- sum(dens_t)
      logvals[k2] <- log(pi_vec[k2] + 1e-16) + sumLog
    }
    # log-sum-exp
    m0 <- max(logvals)
    li <- m0 + log(sum(exp(logvals - m0)))
    logLik_val <- logLik_val + li

    for(k2 in seq_len(K)){
      resp[i,k2] <- exp(logvals[k2] - li)
    }
  }

  # For soft EM, use gamma from C++ as responsibilities if available
  if(em_type == "soft" && !is.null(fit_cpp$gamma)){
    resp <- fit_cpp$gamma
  }

  # 5) Return final structure
  #    **Include $W in the output** so compute_factor_scores_mpca() won't complain
  out <- list(
    z      = z_cpp,
    pi     = pi_vec,
    mu     = mu_list_cpp,
    W      = W_list_cpp,  # <--- keep the original loadings
    P      = P_list,
    D      = D_list,
    Psi    = Psi_list,
    logLik = logLik_val,
    sigma2 = sigma2_vec,
    resp   = resp
  )
  return(out)
}


#' Fit a Mixture PCA Model (Using mpcaTimeseriesCpp) with optional multi-initialization
#'
#' This function calls the C++ function \code{mpcaTimeseriesCpp()} to perform a
#' Mixture PCA (PPCA) EM algorithm for fixed \code{K} and \code{r}, then converts
#' the \code{W} matrices into \code{(P, D)}, creates \code{Psi}, and computes the final
#' log-likelihood/responsibilities in R. By default, it runs a single pass with
#' the built-in initialization in C++. However, if \code{n_init > 1} or
#' \code{use_kmeans_init=TRUE}, it performs multiple initializations in parallel
#' and returns the best-fitting result (highest log-likelihood).
#'
#' @param list_of_data A list of matrices (each \code{T_i x M}).
#' @param K Number of clusters.
#' @param r Number of principal components.
#' @param max_iter Maximum EM iterations for the C++ routine.
#' @param iter_pca Sub-iterations for updating each cluster's PPCA parameters.
#'   (Replaces deprecated \code{nIterPCA})
#' @param tol Convergence tolerance for \code{mpcaTimeseriesCpp}.
#' @param method Either \code{"EM"} or \code{"closed_form"}, passed down to the C++ routine.
#' @param n_init Integer; how many random initial assignments to try (in addition to
#'   the default single-run or k-means if requested). Defaults to \code{1}.
#' @param use_kmeans_init Logical; if \code{TRUE}, we also run one initialization
#'   where we assign clusters by k-means on subject-level PCA features. Defaults to \code{FALSE}.
#' @param kmeans_rdim The PCA dimension for the subject-level feature extraction,
#'   used only if \code{use_kmeans_init=TRUE}. Defaults to \code{r}.
#'   (Replaces deprecated \code{subject_rdim_for_kmeans})
#' @param mc_cores Number of cores for parallel execution via \code{mclapply}.
#'   Defaults to \code{1} (no parallel).
#' @param seed Random seed for reproducibility (default: NULL).
#' @param verbose Logical; if \code{TRUE}, print progress messages during EM iterations. Default is \code{TRUE}.
#' @param n_threads Integer; number of OpenMP threads to use for parallel computation (0 = auto). Default is \code{0}.
#' @param cpp_seed Integer; random seed for C++ RNG (0 = default/no seed). Default is \code{0}.
#' @param progress_callback Optional function to report progress from C++. Should accept two arguments: current iteration and total iterations.
#' @param ... Additional arguments for backwards compatibility. Deprecated arguments like
#'   \code{nIterPCA} and \code{subject_rdim_for_kmeans} are accepted with warnings.
#'
#' @return A list with elements:
#'   \item{z}{Hard cluster assignments (1..K).}
#'   \item{pi}{Mixing proportions.}
#'   \item{mu}{List of length K, each a mean vector.}
#'   \item{P}{List of length K, each \code{(M x r)} principal directions.}
#'   \item{D}{List of length K, each \code{(r x r)} diagonal.}
#'   \item{Psi}{List of length K, each \code{(M x M)} diagonal.}
#'   \item{logLik}{Final log-likelihood over all data.}
#'   \item{sigma2}{Numeric vector of length K for \code{sigma2}.}
#'   \item{resp}{\code{(N x K)} matrix of responsibilities.}
#'
#' @details
#' If \code{n_init=1} and \code{use_kmeans_init=FALSE}, this function runs exactly one pass
#' with the default initialization in C++ (i.e. subject \code{i} is assigned to cluster
#' \code{(i \% K) + 1}). Otherwise:
#' \enumerate{
#'   \item If \code{use_kmeans_init=TRUE}, we do a run where we apply k-means to some
#'         subject-level PCA features to get \code{z_init}, then run EM once.
#'   \item We generate \code{n_init} random initial assignments, run EM for each in parallel,
#'         and collect the results.
#'   \item We compare all solutions by final log-likelihood and pick the best.
#' }
#'
#' @examples
#' \dontrun{
#' # Suppose we have a list_of_data, K=3, r=2, and we want to try 5 random inits + k-means:
#' fit <- mixture_pca_em_fit(
#'   list_of_data, K=3, r=2, max_iter=50, iter_pca=20, tol=1e-3, method="EM",
#'   n_init=5, use_kmeans_init=TRUE, kmeans_rdim=2, mc_cores=2
#' )
#' print(fit$logLik)
#' head(fit$z)
#' }
#'
#' @seealso \code{\link{mfa_em_fit}} for MFA fitting,
#'   \code{\link{select_r}} for selecting optimal r,
#'   \code{\link{plot_cluster_synergy_loadings_mpca}} for visualization
#'
#' @export
mixture_pca_em_fit <- function(list_of_data,
                               K,
                               r,
                               max_iter = 50,
                               iter_pca = 20,
                               tol      = 1e-3,
                               method   = "EM",
                               n_init   = 1,
                               use_kmeans_init = TRUE,
                               kmeans_rdim = r,
                               mc_cores = 1,
                               seed = NULL,
                               verbose = TRUE,
                               n_threads = 0,
                               cpp_seed = 0,
                               progress_callback = NULL,
                               ...)
{
  # Handle deprecated argument names
  dots <- list(...)
  if ("nIterPCA" %in% names(dots)) {
    warning("Argument 'nIterPCA' is deprecated. Use 'iter_pca' instead.", call. = FALSE)
    if (missing(iter_pca)) iter_pca <- dots$nIterPCA
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
    "iter_pca must be positive" = is.numeric(iter_pca) && iter_pca > 0,
    "method must be 'EM' or 'closed_form'" = method %in% c("EM", "closed_form"),
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
  
  # If n_init=1 and use_kmeans_init=FALSE => original single-run approach
  if(n_init == 1 && !use_kmeans_init){
    # Just call mpcaTimeseriesCpp once with no z_init => same as original
    fit_cpp <- mpcaTimeseriesCpp(
      list_of_data = list_of_data,
      K = K,
      r = r,
      max_iter  = max_iter,
      nIterPCA  = iter_pca,
      tol       = tol,
      method    = method,
      verbose   = verbose,
      n_threads = n_threads,
      seed      = cpp_seed,
      progress_callback = progress_callback
    )

    # Then do the original interpretation (W->(P,D), etc.) + logLik calculation
    # (Identical to your original code in mixture_pca_em_fit)
    N <- length(list_of_data)
    if(N < 1){
      stop("No data in list_of_data.")
    }
    W_list_cpp <- fit_cpp$W
    mu_list_cpp<- fit_cpp$mu
    sigma2_vec <- fit_cpp$sigma2
    pi_vec     <- fit_cpp$pi
    z_cpp      <- fit_cpp$z

    if(is.null(sigma2_vec)){
      stop("mpcaTimeseriesCpp did not return 'sigma2'.")
    }

    # Convert W->(P,D), build Psi
    K_check <- length(W_list_cpp)
    if(K_check != K){
      stop("Mismatch in K: length(W_list_cpp) != K.")
    }
    P_list  <- vector("list", K)
    D_list  <- vector("list", K)
    Psi_list<- vector("list", K)

    for(k2 in seq_len(K)){
      W_k <- W_list_cpp[[k2]]
      M_k <- nrow(W_k)
      r_k <- ncol(W_k)
      col_scales <- numeric(r_k)
      for(j in seq_len(r_k)){
        cs_j <- sqrt(sum(W_k[,j]^2))
        if(cs_j < 1e-12) cs_j <- 1e-12
        col_scales[j] <- cs_j
      }
      P_k <- W_k
      for(j in seq_len(r_k)){
        P_k[, j] <- W_k[, j] / col_scales[j]
      }
      D_k <- diag(col_scales, r_k, r_k)
      sig2_k <- sigma2_vec[k2]
      Psi_k  <- diag(sig2_k, M_k)

      P_list[[k2]]  <- P_k
      D_list[[k2]]  <- D_k
      Psi_list[[k2]]<- Psi_k
    }

    # Compute final logLik & resp
    Sigma_list <- vector("list", K)
    for(k2 in seq_len(K)){
      P_k   <- P_list[[k2]]
      D_k   <- D_list[[k2]]
      sig2_k<- sigma2_vec[k2]
      Sig_k <- P_k %*% (D_k^2) %*% t(P_k)
      diag(Sig_k) <- diag(Sig_k) + sig2_k
      Sigma_list[[k2]] <- Sig_k
    }
    logLik_val <- 0
    resp <- matrix(0, nrow=N, ncol=K)
    for(i in seq_len(N)){
      Xi <- list_of_data[[i]]
      logvals <- numeric(K)
      for(k2 in seq_len(K)){
        mu_k  <- mu_list_cpp[[k2]]
        Sig_k <- Sigma_list[[k2]]
        dens_t <- mvtnorm::dmvnorm(Xi, mean=mu_k, sigma=Sig_k, log=TRUE)
        sumLog <- sum(dens_t)
        logvals[k2] <- log(pi_vec[k2] + 1e-16) + sumLog
      }
      m0 <- max(logvals)
      li <- m0 + log(sum(exp(logvals - m0)))
      logLik_val <- logLik_val + li
      for(k2 in seq_len(K)){
        resp[i,k2] <- exp(logvals[k2] - li)
      }
    }

    out <- list(
      z      = z_cpp,
      pi     = pi_vec,
      mu     = mu_list_cpp,
      P      = P_list,
      D      = D_list,
      Psi    = Psi_list,
      logLik = logLik_val,
      sigma2 = sigma2_vec,
      resp   = resp
    )
    return(out)
  }

  # Otherwise => multi-init approach
  best_fit    <- NULL
  best_logLik <- -Inf

  # 1) If k-means init is requested => run once
  if(use_kmeans_init){
    if(!requireNamespace("parallel", quietly=TRUE)){
      stop("Package 'parallel' is required for mclapply. Please install it.")
    }
    features <- extract_subject_features_by_singlePCA(
      list_of_data, r_dim=kmeans_rdim
    )
    z_init_km <- assign_by_kmeans(features, K=K)

    fit_km <- mixture_pca_em_fit_cpp_singleInit(
      list_of_data = list_of_data,
      K = K,
      r = r,
      max_iter  = max_iter,
      nIterPCA  = iter_pca,
      tol       = tol,
      method    = method,
      z_init    = z_init_km,
      verbose   = verbose,
      n_threads = n_threads,
      seed      = cpp_seed,
      progress_callback = progress_callback
    )
    if(fit_km$logLik > best_logLik){
      best_fit    <- fit_km
      best_logLik <- fit_km$logLik
    }
  }

  # 2) Random initializations => parallel
  N <- length(list_of_data)
  random_inits <- lapply(seq_len(n_init), function(i) {
    sample.int(K, size=N, replace=TRUE)
  })

  if(!requireNamespace("parallel", quietly=TRUE)){
    stop("Package 'parallel' is required for mclapply. Please install it.")
  }
  fit_list <- parallel::mclapply(
    random_inits,
    FUN = function(z_init_rand){
      fit_rand <- mixture_pca_em_fit_cpp_singleInit(
        list_of_data = list_of_data,
        K = K,
        r = r,
        max_iter  = max_iter,
        nIterPCA  = iter_pca,
        tol       = tol,
        method    = method,
        z_init    = z_init_rand,
        verbose   = verbose,
        n_threads = n_threads,
        seed      = cpp_seed,
        progress_callback = progress_callback
      )
      fit_rand
    },
    mc.cores = mc_cores
  )

  # Compare random solutions
  for(fit_rand in fit_list){
    if(fit_rand$logLik > best_logLik){
      best_fit    <- fit_rand
      best_logLik <- fit_rand$logLik
    }
  }

  best_fit
  
  })  # End of with_internal_threads()
  
  return(result)
}

