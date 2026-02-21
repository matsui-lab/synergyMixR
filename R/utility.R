#' Compute Global VAF (1 - SSE/SST) for an MFA Model
#'
#' This function takes a fitted MFA model (with elements \code{z, mu, Lambda, Psi})
#' and reconstructs all subjects' data from the factor scores, computing the
#' ratio of explained variance (\eqn{1 - \mathrm{SSE}/\mathrm{SST}}).
#'
#' @param list_of_data A list of matrices, each \code{(T_i x M)}, containing the original data.
#' @param res A fitted MFA model result. It must contain:
#'   \itemize{
#'     \item{\code{z}: A length-N vector of cluster assignments (1..K).}
#'     \item{\code{mu}: A list of length K, each a mean vector of length M.}
#'     \item{\code{Lambda}: A list of length K, each an \code{M x r} loading matrix.}
#'     \item{\code{Psi}: A list of length K, each an \code{M x M} diagonal matrix.}
#'   }
#' @return A numeric scalar for the global variance accounted for, \eqn{1 - \mathrm{SSE}/\mathrm{SST}}.
#'
#' @details
#' For each subject \code{i}, we identify its cluster assignment \code{k_i}, compute factor scores
#' via \code{\link{compute_factor_scores}}, reconstruct \code{Xhat} from the factor model,
#' accumulate the sum of squares of residuals (SSE), and the sum of squares total (SST).
#'
#' @seealso \code{\link{compute_factor_scores}}
#'
#' @examples
#' \dontrun{
#' gvaf <- compute_global_vaf_mfa(my_data_list, fit_mfa)
#' cat("Global VAF:", gvaf, "\n")
#' }
#'
#' @export
compute_global_vaf_mfa <- function(list_of_data, res){
  # res: assumed to contain z, mu, Lambda, Psi (lists)
  z_vec <- res$z
  N <- length(list_of_data)
  M <- ncol(list_of_data[[1]])

  # Pre-allocate matrices for memory efficiency (avoid O(N^2) rbind)
  total_rows <- sum(vapply(list_of_data, nrow, integer(1)))
  X_all <- matrix(0, nrow = total_rows, ncol = M)
  Xhat_all <- matrix(0, nrow = total_rows, ncol = M)

  idx <- 1
  for(i in seq_len(N)){
    k_i <- z_vec[i]  # cluster ID for subject i
    mu_k <- res$mu[[k_i]]
    Lambda_k <- res$Lambda[[k_i]]
    Psi_k <- res$Psi[[k_i]]

    # observed data
    X_i <- list_of_data[[i]]  # (T_i x M)
    T_i <- nrow(X_i)

    # factor scores
    EZ_i <- compute_factor_scores(X_i, mu_k, Lambda_k, Psi_k)

    # reconstruct
    Xhat_i <- EZ_i %*% t(Lambda_k)
    Xhat_i <- sweep(Xhat_i, 2, mu_k, FUN="+")

    # Assign to pre-allocated matrices
    rows <- idx:(idx + T_i - 1)
    X_all[rows, ] <- X_i
    Xhat_all[rows, ] <- Xhat_i
    idx <- idx + T_i
  }

  SSE_mat <- (X_all - Xhat_all)^2
  SSE <- sum(SSE_mat)
  SST <- sum(X_all^2)

  vaf_global <- 1 - SSE/SST
  vaf_global
}


#' Compute Factor Scores for a Single MFA Cluster
#'
#' Given a single cluster's parameters (\code{mu, Lambda, Psi}) and data \code{X},
#' this function returns the estimated factor scores.
#'
#' @param X A numeric matrix of size \code{(T x M)}.
#' @param mu A length-M mean vector.
#' @param Lambda An \code{M x r} loading matrix.
#' @param Psi An \code{M x M} diagonal matrix for unique variances.
#'
#' @return A \code{(T x r)} matrix of factor scores.
#'
#' @details
#' This follows the standard posterior mean formula in factor analysis:
#' \eqn{\mathrm{E}[Z|X] = (I + \Lambda^T \Psi^{-1} \Lambda)^{-1} \Lambda^T \Psi^{-1} (X - \mu)}.
#'
#' @seealso \code{\link{compute_global_vaf_mfa}}
#'
#' @examples
#' \dontrun{
#' # Suppose we have X (T x M), mu (M), Lambda (M x r), Psi (M x M).
#' scores <- compute_factor_scores(X, mu, Lambda, Psi)
#' head(scores)
#' }
#'
#' @export
compute_factor_scores <- function(X, mu, Lambda, Psi){
  # Wrapper function that calls the optimized implementation
  compute_factor_scores_optimized(X, mu, Lambda, Psi)
}


#' Compute Cluster Sizes from a Fitted MFA (or Mixture PCA) Model
#'
#' Given a result with a vector \code{z} of length N (the cluster assignments),
#' this function returns a frequency table, i.e., how many subjects are assigned
#' to each cluster.
#'
#' @param res A fitted mixture model object that contains \code{z}.
#'
#' @return A table of cluster frequencies.
#'
#' @examples
#' \dontrun{
#' tab <- compute_cluster_sizes(fit_mfa)
#' print(tab)
#' }
#' @export
compute_cluster_sizes <- function(res){
  z_vec <- res$z
  tab <- table(z_vec)
  tab
}


#' Post-hoc Evaluation of MFA Grid Search
#'
#' Given the output of \code{\link{select_optimal_K_r_mfa}}, this function computes
#' additional metrics for each fitted model in \code{all_models}, such as
#' global VAF (\code{\link{compute_global_vaf_mfa}}) and cluster sizes
#' (\code{\link{compute_cluster_sizes}}).
#'
#' @param list_of_data The same list of data used to fit the models.
#' @param selection_obj The object returned by \code{\link{select_optimal_K_r_mfa}}.
#'
#' @return A data frame with columns:
#' \item{K}{Number of clusters.}
#' \item{r}{Factor dimension.}
#' \item{logLik}{Log-likelihood from the fit.}
#' \item{BIC}{BIC from the fit.}
#' \item{GlobalVAF}{Global variance accounted for, \eqn{1 - \mathrm{SSE}/\mathrm{SST}}.}
#' \item{minClusterSize}{Minimum cluster size.}
#' \item{maxClusterSize}{Maximum cluster size.}
#'
#' @details
#' This function loops over each entry in \code{selection_obj$all_models},
#' extracts the model, calls \code{\link{compute_global_vaf_mfa}} and
#' \code{\link{compute_cluster_sizes}}, and records the results in a combined table.
#'
#' @seealso \code{\link{select_optimal_K_r_mfa}}, \code{\link{compute_factor_scores}}
#'
#' @examples
#' \dontrun{
#' sel_obj <- select_optimal_K_r_mfa(my_data, Kvec=2:4, rvec=1:3)
#' df_posthoc <- posthoc_mfa_evaluation(my_data, sel_obj)
#' head(df_posthoc)
#' }
#'
#' @export
posthoc_mfa_evaluation <- function(list_of_data, selection_obj){
  all_mods <- selection_obj$all_models
  if(is.null(all_mods) || length(all_mods)==0){
    stop("No models found in selection_obj$all_models")
  }
  df_list <- list()

  for(i in seq_along(all_mods)){
    obj_i <- all_mods[[i]]
    K_val <- obj_i$K
    r_val <- obj_i$r
    fit   <- obj_i$model

    # 1) Global VAF
    gvaf_i <- compute_global_vaf_mfa(list_of_data, fit)
    # 2) cluster sizes
    size_tab <- compute_cluster_sizes(fit)
    min_size <- min(size_tab)
    max_size <- max(size_tab)

    df_list[[i]] <- data.frame(
      K    = K_val,
      r    = r_val,
      logLik = obj_i$logLik,
      BIC    = obj_i$BIC,
      GlobalVAF = gvaf_i,
      minClusterSize = min_size,
      maxClusterSize = max_size
    )
  }

  df_posthoc <- dplyr::bind_rows(df_list)
  df_posthoc <- df_posthoc[order(df_posthoc$BIC), ]
  rownames(df_posthoc) <- NULL
  df_posthoc
}


#' Retrieve a Specific Model by (K, r) from a Grid Search
#'
#' After calling \code{\link{select_optimal_K_r_mfa}} (or a similar function),
#' you have \code{$all_models} containing fits for various \code{(K, r)}. This
#' function extracts the model for a user-specified pair \code{(K_target, r_target)}.
#'
#' @param selection_obj The object returned by \code{\link{select_optimal_K_r_mfa}}.
#' @param K_target Desired number of clusters.
#' @param r_target Desired factor dimension.
#'
#' @return The fitted model (the same structure returned by \code{mfa_em_fit}),
#'   or \code{NULL} if no match found.
#'
#' @examples
#' \dontrun{
#' sel_obj <- select_optimal_K_r_mfa(my_data, Kvec=2:4, rvec=1:3)
#' best_model <- get_model_by_K_r(sel_obj, K_target=3, r_target=2)
#' }
#'
#' @export
get_model_by_K_r <- function(selection_obj, K_target, r_target){
  all_mods <- selection_obj$all_models
  if(is.null(all_mods) || length(all_mods)==0) return(NULL)

  idx <- which(
    sapply(all_mods, function(x) x$K) == K_target &
      sapply(all_mods, function(x) x$r) == r_target
  )
  if(length(idx)<1){
    message(sprintf("No model found for K=%d, r=%d.", K_target, r_target))
    return(NULL)
  }
  all_mods[[ idx[1] ]]$model
}



#' Compute Global VAF (1 - SSE/SST) for Mixture PCA
#'
#' Similar to \code{\link{compute_global_vaf_mfa}}, but uses a Mixture PCA model structure.
#' Given a model with \code{z, mu, P, D, Psi}, we reconstruct each subject's data from
#' the principal component scores, accumulate SSE, and compare to total sum of squares (SST).
#'
#' @param list_of_data A list of \code{(T_i x M)} data matrices.
#' @param res The fitted Mixture PCA model, typically from \code{\link{mixture_pca_em_fit}}.
#'   Must contain:
#'   \itemize{
#'     \item{\code{z}: cluster assignments (length N).}
#'     \item{\code{mu}: list of length K, each \code{(M)} vector.}
#'     \item{\code{P}: list of length K, each \code{(M x r)}.}
#'     \item{\code{D}: list of length K, each \code{(r x r)} diagonal.}
#'     \item{\code{Psi}: list of length K, each \code{(M x M)} diagonal.}
#'   }
#'
#' @return A numeric scalar, \eqn{1 - \mathrm{SSE}/\mathrm{SST}} across all subjects and time steps.
#'
#' @details
#' The user must define how to reconstruct data from PCA parameters. Typically,
#' for each cluster \code{k}, we have \code{W = P * D}, so we either compute
#' posterior scores or do a direct projection.
#'
#' @export
compute_global_vaf_mpca <- function(list_of_data, res){
  z_vec <- res$z
  N <- length(list_of_data)
  M <- ncol(list_of_data[[1]])

  # Pre-allocate matrices for memory efficiency (avoid O(N^2) rbind)
  total_rows <- sum(vapply(list_of_data, nrow, integer(1)))
  X_all <- matrix(0, nrow = total_rows, ncol = M)
  Xhat_all <- matrix(0, nrow = total_rows, ncol = M)

  idx <- 1
  for(i in seq_len(N)){
    k_i <- z_vec[i]  # cluster ID for subject i
    mu_k <- res$mu[[k_i]]
    W_k <- res$W[[k_i]]  # or res$P[[k_i]] depending on structure

    # Handle different possible structures
    if (is.null(W_k) && !is.null(res$P)) {
      P_k <- res$P[[k_i]]
      D_k <- res$D[[k_i]]
      W_k <- P_k %*% D_k  # Reconstruct W from P and D
    }

    if (is.null(W_k)) {
      stop("Neither W nor P matrices found in MPCA model structure")
    }

    # observed data
    X_i <- list_of_data[[i]]  # (T_i x M)
    T_i <- nrow(X_i)

    # For MPCA, we need to compute factor scores differently
    # Using the posterior mean for PPCA: z = (W^T W + sigma2 I)^{-1} W^T (x - mu)
    sigma2_k <- if (!is.null(res$sigma2)) res$sigma2[k_i] else 1.0

    r <- ncol(W_k)

    # Center the data
    Xc_i <- sweep(X_i, 2, mu_k, FUN = "-")

    # Compute posterior factor scores for PPCA (vectorized)
    A <- t(W_k) %*% W_k + sigma2_k * diag(r)
    A_inv <- solve(A)
    W_t <- t(W_k)

    # Factor scores: Z_i = Xc_i %*% W %*% A_inv^T (vectorized over time)
    Z_i <- Xc_i %*% W_k %*% t(A_inv)

    # reconstruct: Xhat = Z * W^T + mu
    Xhat_i <- Z_i %*% t(W_k)
    Xhat_i <- sweep(Xhat_i, 2, mu_k, FUN = "+")

    # Assign to pre-allocated matrices
    rows <- idx:(idx + T_i - 1)
    X_all[rows, ] <- X_i
    Xhat_all[rows, ] <- Xhat_i
    idx <- idx + T_i
  }

  SSE_mat <- (X_all - Xhat_all)^2
  SSE <- sum(SSE_mat)
  SST <- sum(X_all^2)

  vaf_global <- 1 - SSE/SST
  vaf_global
}


#' Compute Cluster Sizes for Mixture PCA
#'
#' Like \code{\link{compute_cluster_sizes}}, but for a Mixture PCA model result.
#'
#' @param res A fitted mixture PCA model result, containing \code{$z}.
#'
#' @return A table of frequencies of cluster assignments.
#'
#' @export
compute_cluster_sizes_mpca <- function(res){
  # same as MFA, basically
  z_vec <- res$z
  table(z_vec)
}


#' Post-hoc Evaluation of Mixture PCA Grid Search
#'
#' Takes the output of \code{\link{select_optimal_K_r_mpca}} and computes additional
#' metrics (global VAF, cluster sizes, etc.) for each entry in \code{all_models}.
#'
#' @param list_of_data The same list of data used for the grid search.
#' @param selection_obj The object returned by \code{\link{select_optimal_K_r_mpca}}.
#'
#' @return A data frame with columns \code{K, r, logLik, BIC, GlobalVAF, minClusterSize, maxClusterSize}.
#'
#' @export
posthoc_mpca_evaluation <- function(list_of_data, selection_obj){
  all_mods <- selection_obj$all_models
  if(is.null(all_mods) || length(all_mods)==0){
    stop("No models found in selection_obj$all_models")
  }

  df_list <- list()

  for(i in seq_along(all_mods)){
    obj_i <- all_mods[[i]]
    K_val <- obj_i$K
    r_val <- obj_i$r
    fit   <- obj_i$model  # mixture PCA fit

    # compute global VAF
    gvaf_i <- compute_global_vaf_mpca(list_of_data, fit)

    # cluster sizes
    size_tab <- compute_cluster_sizes_mpca(fit)
    min_size <- min(size_tab)
    max_size <- max(size_tab)

    df_list[[i]] <- data.frame(
      K = K_val,
      r = r_val,
      logLik = obj_i$logLik,
      BIC = obj_i$BIC,
      GlobalVAF = gvaf_i,
      minClusterSize = min_size,
      maxClusterSize = max_size
    )
  }

  df_posthoc <- dplyr::bind_rows(df_list)
  df_posthoc <- df_posthoc[order(df_posthoc$BIC), ]
  rownames(df_posthoc) <- NULL
  df_posthoc
}


#' Retrieve a Mixture PCA Model by (K, r)
#'
#' Given a \code{selection_obj} from \code{\link{select_optimal_K_r_mpca}},
#' extracts the model for the specified \code{(K_target, r_target)} if it exists.
#'
#' @param selection_obj The output of \code{\link{select_optimal_K_r_mpca}}.
#' @param K_target Desired K.
#' @param r_target Desired r.
#'
#' @return The fitted model (same structure as from \code{mixture_pca_em_fit}), or \code{NULL}.
#'
#' @export
get_model_by_K_r_mpca <- function(selection_obj, K_target, r_target){
  all_mods <- selection_obj$all_models
  if(is.null(all_mods) || length(all_mods)==0) return(NULL)

  idx <- which(
    sapply(all_mods, function(x) x$K) == K_target &
      sapply(all_mods, function(x) x$r) == r_target
  )
  if(length(idx)<1){
    message(sprintf("No PCA model found for K=%d, r=%d.", K_target, r_target))
    return(NULL)
  }
  all_mods[[ idx[1] ]]$model
}

# ============================================================
# compute_logLik_mfa
# ============================================================

#' Compute the Mixture-Model Log-Likelihood for MFA
#'
#' Given a fitted MFA model (with \code{K} clusters, each containing
#' \code{Lambda[[k]]}, \code{mu[[k]]}, \code{Psi[[k]]}, and mixing proportion \code{pi[k]}),
#' this function computes the standard mixture log-likelihood:
#' \deqn{
#'   \sum_{i=1}^{N} \log\left(
#'     \sum_{k=1}^{K} \pi_k \prod_{t=1}^{T_i} \mathcal{N}(x_{i,t} \mid \mu_k, \Sigma_k)
#'   \right),
#' }
#' where Sigma_k = Lambda_k * Lambda_k^T + Psi_k.
#'
#' @param list_of_data A list of length \code{N}, each an \code{(T_i x M)} matrix.
#' @param mfa_fit A fitted MFA model, containing:
#'   \itemize{
#'     \item \code{K} clusters (implicit from \code{length(mfa_fit$Lambda)}),
#'     \item \code{Lambda[[k]]}, \code{mu[[k]]}, \code{Psi[[k]]},
#'     \item \code{pi} as a length-K vector of mixing proportions.
#'   }
#'
#' @return A numeric scalar, the total log-likelihood.
#'
#' @details
#' Internally uses rowwise multivariate normal densities. You will need
#' \code{mvtnorm::dmvnorm} or an equivalent. Ensure each Sigma_k is well-defined
#' and invertible if needed.
#'
#' @export
compute_logLik_mfa <- function(list_of_data, mfa_fit) {
  if(!requireNamespace("mvtnorm", quietly=TRUE)) {
    stop("Package 'mvtnorm' is required for compute_logLik_mfa. Please install it.")
  }

  Lambda_list <- mfa_fit$Lambda
  mu_list     <- mfa_fit$mu
  Psi_list    <- mfa_fit$Psi
  pi_vec      <- mfa_fit$pi   # length K

  K <- length(Lambda_list)
  N <- length(list_of_data)

  # Precompute Sigma_k for each cluster
  Sigma_list <- vector("list", K)
  for(k in seq_len(K)) {
    Lam_k <- Lambda_list[[k]]
    Psi_k <- Psi_list[[k]]
    Sig_k <- Lam_k %*% t(Lam_k) + Psi_k
    Sigma_list[[k]] <- Sig_k
  }

  total_loglik <- 0
  for(i in seq_len(N)) {
    X_i <- list_of_data[[i]]
    T_i <- nrow(X_i)
    logvals <- numeric(K)

    # For each cluster k, compute log(pi_k) + sum of rowwise log density
    for(k in seq_len(K)) {
      mu_k  <- mu_list[[k]]
      Sig_k <- Sigma_list[[k]]
      # rowwise log densities
      dens_vec <- mvtnorm::dmvnorm(X_i, mean=mu_k, sigma=Sig_k, log=TRUE)
      sum_logdens <- sum(dens_vec)
      logvals[k] <- log(pi_vec[k] + 1e-16) + sum_logdens
    }
    # mixture log-likelihood for subject i
    # log( sum_k exp(logvals[k]) )
    # use a stable log-sum-exp
    maxv <- max(logvals)
    logLi <- maxv + log(sum(exp(logvals - maxv)))
    total_loglik <- total_loglik + logLi
  }

  total_loglik
}

#' Compute BIC for an MFA model
#'
#' Given the log-likelihood from \code{\link{compute_logLik_mfa}}, and known
#' \code{K}, \code{r}, \code{M}, \code{N}, this function computes
#' BIC = -2 * log(L) + nu * log(N),
#' where nu is the parameter count for the MFA model. The corrected formula is
#' nu = K * (M * r + r * (r-1)/2 + M + M) + (K - 1).
#'
#' @param loglik The total log-likelihood from \code{compute_logLik_mfa}.
#' @param K Number of clusters.
#' @param r Factor dimension.
#' @param M Observed dimension (channels).
#' @param N_obs Number of observations (typically sum of all T_i across subjects).
#'
#' @return A numeric BIC value.
#'
#' @details
#' The parameter count includes:
#' \itemize{
#'   \item \code{M*r}: factor loadings (Lambda)
#'   \item \code{r*(r-1)/2}: factor covariances (off-diagonal elements)
#'   \item \code{M}: unique variances (diagonal Psi)
#'   \item \code{M}: cluster means (mu)
#'   \item \code{K-1}: mixing proportions
#' }
#' This corrected formula properly penalizes larger r values, improving r estimation stability.
#' 
#' Note: N_obs should typically be the total number of time points (sum of T_i) rather than
#' the number of subjects, as BIC penalizes based on sample size.
#'
#' @export
compute_BIC_mfa <- function(loglik, K, r, M, N_obs) {
  # Corrected parameter count including factor covariances
  # For each cluster k:
  #   - M*r: factor loadings (Lambda_k)
  #   - r*(r-1)/2: factor covariances (off-diagonal elements of factor covariance)
  #   - M: unique variances (diagonal of Psi_k)
  #   - M: cluster mean (mu_k)
  # Plus (K-1) for mixing proportions
  param_count <- K * (M*r + r*(r-1)/2 + M + M) + (K - 1)

  BIC_val <- -2 * loglik + param_count * log(N_obs)
  BIC_val
}

# ============================================================
# compute_logLik_mpca
# ============================================================

#' Compute the Mixture-Model Log-Likelihood for Mixture PCA
#'
#' Given a fitted Mixture PCA model with \code{$P[[k]]}, \code{$mu[[k]]},
#' \code{$sigma2[k]}, and \code{$pi[k]}, we compute
#' \deqn{
#'   \sum_{i=1}^{N} \log\Bigl(\sum_{k=1}^{K} \pi_k
#'     \prod_{t=1}^{T_i}\mathcal{N}\bigl(x_{i,t}\mid \mu_k,\Sigma_k\bigr)\Bigr)
#' }
#' where Sigma_k = W_k * W_k^T + sigma2_k * I. A typical approach is to
#' reconstruct Sigma_k from P[[k]] if needed (or W[[k]]).
#'
#' @param list_of_data A list of length \code{N}, each \code{(T_i x M)}.
#' @param mpca_fit A fitted Mixture PCA model, containing:
#'   \itemize{
#'     \item \code{P[[k]]} or \code{W[[k]]} for cluster \code{k},
#'     \item \code{mu[[k]]}, \code{sigma2[k]}, \code{pi[k]}.
#'   }
#'
#' @return A numeric scalar, the total mixture log-likelihood.
#'
#' @details
#' We assume each cluster k has Sigma_k = W_k * W_k^T + sigma2_k * I.
#' You must assemble that to call \code{mvtnorm::dmvnorm} row by row.
#'
#' @export
compute_logLik_mpca <- function(list_of_data, mpca_fit) {
  if(!requireNamespace("mvtnorm", quietly=TRUE)) {
    stop("Package 'mvtnorm' is required for compute_logLik_mpca. Please install it.")
  }
  P_list    <- mpca_fit$P  # or mpca_fit$W
  mu_list   <- mpca_fit$mu
  sigma2vec <- mpca_fit$sigma2
  pi_vec    <- mpca_fit$pi


  K <- length(P_list)
  N <- length(list_of_data)

  # Precompute Sigma_k
  Sigma_list <- vector("list", K)
  for(k in seq_len(K)) {
    P_k  <- P_list[[k]]  # (M x r)
    M    <- nrow(P_k)
    sig2 <- sigma2vec[k]
    # W_k W_k^T + sig2 I
    # here W_k = P_k * some scale? or if P_k is the final W? adapt if needed
    # We'll treat P_k as W_k directly, you might adapt code if you have D_k
    Sig_k <- P_k %*% t(P_k)
    diag(Sig_k) <- diag(Sig_k) + sig2
    Sigma_list[[k]] <- Sig_k
  }

  total_loglik <- 0
  for(i in seq_len(N)) {
    X_i <- list_of_data[[i]]
    T_i <- nrow(X_i)
    logvals <- numeric(K)

    for(k in seq_len(K)) {
      mu_k  <- mu_list[[k]]
      Sig_k <- Sigma_list[[k]]
      dens_vec <- mvtnorm::dmvnorm(X_i, mean=mu_k, sigma=Sig_k, log=TRUE)
      sum_logdens <- sum(dens_vec)
      logvals[k] <- log(pi_vec[k] + 1e-16) + sum_logdens
    }
    maxv <- max(logvals)
    logLi <- maxv + log(sum(exp(logvals - maxv)))
    total_loglik <- total_loglik + logLi
  }

  total_loglik
}

#' Compute BIC for a Mixture PCA Model
#'
#' Given the log-likelihood from \code{\link{compute_logLik_mpca}}, and the
#' number of clusters \code{K}, dimension \code{r}, etc., this function
#' computes the corrected BIC:
#' BIC = -2 * log(L) + nu * log(N),
#' where nu is the parameter count for MPCA. The corrected formula is
#' \code{K*(M*r + r*(r-1)/2 + M + 1) + (K-1)}, i.e. for each cluster
#' we have \code{(M*r)} for the directions, \code{r*(r-1)/2} for factor covariances,
#' \code{M} for the mean, and \code{1} for sigma^2. Then \code{(K-1)} for mixing proportions.
#'
#' @param loglik The total mixture log-likelihood from \code{compute_logLik_mpca}.
#' @param K Number of clusters.
#' @param r Dimension of principal components.
#' @param M Observed dimension.
#' @param N_obs Number of observations (typically sum of all T_i across subjects).
#'
#' @return A numeric BIC value.
#' 
#' @details
#' The parameter count includes:
#' \itemize{
#'   \item \code{M*r}: loading matrix (W or P)
#'   \item \code{r*(r-1)/2}: factor covariances (off-diagonal elements)
#'   \item \code{M}: cluster means (mu)
#'   \item \code{1}: isotropic noise variance (sigma2)
#'   \item \code{K-1}: mixing proportions
#' }
#' This corrected formula properly penalizes larger r values, improving r estimation stability.
#' 
#' Note: N_obs should typically be the total number of time points (sum of T_i) rather than
#' the number of subjects, as BIC penalizes based on sample size.
#' 
#' @export
compute_BIC_mpca <- function(loglik, K, r, M, N_obs) {
  # Corrected parameter count including factor covariances
  # For each cluster k:
  #   - M*r: loading matrix (W_k or P_k)
  #   - r*(r-1)/2: factor covariances (off-diagonal elements)
  #   - M: cluster mean (mu_k)
  #   - 1: isotropic noise variance (sigma2_k)
  # Plus (K-1) for mixing proportions
  param_count <- K * (M*r + r*(r-1)/2 + M + 1) + (K - 1)
  BIC_val <- -2 * loglik + param_count * log(N_obs)
  BIC_val
}


# ============================================================
# ICL (Integrated Complete-data Likelihood) Functions
# ============================================================

#' Compute Classification Entropy from Responsibilities
#'
#' Computes the entropy term used in ICL (Integrated Complete-data Likelihood).
#' Higher entropy indicates more uncertain cluster assignments.
#'
#' @param resp An (N x K) matrix of responsibilities (posterior probabilities),
#'   where resp[i,k] is the probability that subject i belongs to cluster k.
#'
#' @return A numeric scalar representing the entropy: -sum(resp * log(resp)).
#'
#' @details
#' The entropy is computed as:
#' \deqn{H = -\sum_{i=1}^{N} \sum_{k=1}^{K} r_{ik} \log(r_{ik})}
#' where \eqn{r_{ik}} is the responsibility of subject i for cluster k.
#' A small constant (1e-16) is added to avoid log(0).
#'
#' @examples
#' \dontrun{
#' # Responsibilities from MFA fit
#' entropy <- compute_entropy(fit_mfa$resp)
#' }
#'
#' @export
compute_entropy <- function(resp) {
  # Add small constant to avoid log(0)
  resp_safe <- resp + 1e-16
  entropy <- -sum(resp_safe * log(resp_safe))
  entropy
}

#' Compute ICL for a Mixture Factor Analysis Model
#'
#' ICL (Integrated Complete-data Likelihood) is a model selection criterion
#' that penalizes uncertain cluster assignments more heavily than BIC.
#' This makes it more conservative for selecting the number of clusters K.
#'
#' @param loglik The total log-likelihood from the MFA model.
#' @param K Number of clusters.
#' @param r Factor dimension.
#' @param M Observed dimension (channels).
#' @param N_obs Number of observations (typically sum of all T_i across subjects).
#' @param resp An (N x K) matrix of responsibilities (posterior probabilities).
#'
#' @return A numeric ICL value.
#'
#' @details
#' ICL is computed as:
#' \deqn{ICL = BIC + 2 \times Entropy}
#' where Entropy = -sum(resp * log(resp)).
#'
#' The entropy penalty discourages models where cluster assignments are uncertain,
#' which typically happens when K is overestimated (extra clusters with low membership).
#' This makes ICL more conservative than BIC for K selection.
#'
#' @references
#' Biernacki, C., Celeux, G., & Govaert, G. (2000). Assessing a mixture model
#' for clustering with the integrated completed likelihood. IEEE Transactions
#' on Pattern Analysis and Machine Intelligence, 22(7), 719-725.
#'
#' @seealso \code{\link{compute_BIC_mfa}}, \code{\link{compute_entropy}}
#'
#' @examples
#' \dontrun{
#' icl <- compute_ICL_mfa(fit$logLik, K = 2, r = 3, M = 8,
#'                        N_obs = 1000, resp = fit$resp)
#' }
#'
#' @export
compute_ICL_mfa <- function(loglik, K, r, M, N_obs, resp) {
  # Compute BIC first
  bic_val <- compute_BIC_mfa(loglik, K, r, M, N_obs)

  # Compute entropy penalty
  entropy <- compute_entropy(resp)

  # ICL = BIC + 2 * entropy (higher entropy = worse ICL)
  icl_val <- bic_val + 2 * entropy
  icl_val
}

#' Compute ICL for a Mixture PCA Model
#'
#' ICL (Integrated Complete-data Likelihood) is a model selection criterion
#' that penalizes uncertain cluster assignments more heavily than BIC.
#' This makes it more conservative for selecting the number of clusters K.
#'
#' @param loglik The total log-likelihood from the MPCA model.
#' @param K Number of clusters.
#' @param r Principal component dimension.
#' @param M Observed dimension (channels).
#' @param N_obs Number of observations (typically sum of all T_i across subjects).
#' @param resp An (N x K) matrix of responsibilities (posterior probabilities).
#'
#' @return A numeric ICL value.
#'
#' @details
#' ICL is computed as:
#' \deqn{ICL = BIC + 2 \times Entropy}
#' where Entropy = -sum(resp * log(resp)).
#'
#' @seealso \code{\link{compute_BIC_mpca}}, \code{\link{compute_entropy}}
#'
#' @export
compute_ICL_mpca <- function(loglik, K, r, M, N_obs, resp) {
  # Compute BIC first
  bic_val <- compute_BIC_mpca(loglik, K, r, M, N_obs)

  # Compute entropy penalty
  entropy <- compute_entropy(resp)

  # ICL = BIC + 2 * entropy (higher entropy = worse ICL)
  icl_val <- bic_val + 2 * entropy
  icl_val
}

#' Compute ICL for a Mixture Model (Unified Interface)
#'
#' Unified interface for computing ICL (Integrated Complete-data Likelihood)
#' for either MFA or MPCA models.
#'
#' @param loglik The total log-likelihood from the fitted model.
#' @param K Number of clusters.
#' @param r Factor/principal component dimension.
#' @param M Observed dimension (channels).
#' @param N_obs Number of observations.
#' @param resp An (N x K) matrix of responsibilities.
#' @param model_type Either "MFA" or "MPCA".
#'
#' @return A numeric ICL value.
#'
#' @seealso \code{\link{compute_ICL_mfa}}, \code{\link{compute_ICL_mpca}}
#'
#' @export
compute_ICL <- function(loglik, K, r, M, N_obs, resp, model_type) {
  model_type <- toupper(model_type)
  if (model_type == "MFA") {
    compute_ICL_mfa(loglik, K, r, M, N_obs, resp)
  } else if (model_type == "MPCA") {
    compute_ICL_mpca(loglik, K, r, M, N_obs, resp)
  } else {
    stop("model_type must be 'MFA' or 'MPCA'")
  }
}


#' Select K Using Bootstrap Stability Selection
#'
#' This function implements bootstrap stability selection for choosing the optimal
#' number of clusters K. It assesses the stability of cluster assignments across
#' bootstrap resamples and selects K that produces the most stable clustering.
#'
#' @param list_of_data A list of data matrices, each (T_i x M).
#' @param Kvec Vector of candidate K values to test.
#' @param r Number of factors (fixed).
#' @param B Number of bootstrap iterations (default: 50).
#' @param subsample_ratio Proportion of subjects to sample (default: 0.8).
#' @param max_iter Maximum EM iterations per fit.
#' @param n_init Number of random initializations per fit.
#' @param model_type Either "MFA" or "MPCA" (default: "MFA").
#' @param verbose Print progress messages (default: TRUE).
#'
#' @return A list containing:
#' \describe{
#'   \item{best_K}{The selected optimal K based on stability}
#'   \item{stability_scores}{Stability score for each K}
#'   \item{stability_df}{Data frame with K and stability scores}
#' }
#'
#' @details
#' The stability selection procedure:
#' 1. For each candidate K, perform B bootstrap iterations
#' 2. In each iteration, subsample subjects and fit the model
#' 3. Compute pairwise co-clustering matrix (1 if subjects i,j in same cluster)
#' 4. Average co-clustering matrices across bootstrap samples
#' 5. Compute stability score as mean of (max probability per subject - 1/K)
#' 6. Select K with highest stability score
#'
#' Higher stability indicates more consistent cluster assignments across resamples.
#'
#' @references
#' Ben-Hur, A., Elisseeff, A., & Guyon, I. (2001). A stability based method for
#' discovering structure in clustered data. Pacific Symposium on Biocomputing.
#'
#' @export
select_K_bootstrap_stability <- function(
    list_of_data,
    Kvec = 2:5,
    r = 2,
    B = 50,
    subsample_ratio = 0.8,
    max_iter = 50,
    n_init = 1,
    model_type = c("MFA", "MPCA"),
    verbose = TRUE
) {
  model_type <- match.arg(model_type)
  N <- length(list_of_data)
  n_subsample <- max(2, floor(N * subsample_ratio))

  if (verbose) {
    cat("=== Bootstrap Stability Selection for K ===\n")
    cat(sprintf("N subjects: %d, Subsample size: %d, B iterations: %d\n",
                N, n_subsample, B))
    cat(sprintf("Model type: %s, r fixed at: %d\n\n", model_type, r))
  }

  stability_scores <- numeric(length(Kvec))
  names(stability_scores) <- as.character(Kvec)

  for (ki in seq_along(Kvec)) {
    K <- Kvec[ki]
    if (verbose) cat(sprintf("Testing K=%d...\n", K))

    # Store cluster assignments for each bootstrap sample
    # Each row is a bootstrap sample, each column is a subject
    # Value is cluster assignment (NA if not sampled)
    assignments_matrix <- matrix(NA, nrow = B, ncol = N)

    for (b in seq_len(B)) {
      # Subsample subjects
      sampled_idx <- sort(sample(N, n_subsample, replace = FALSE))
      subsampled_data <- list_of_data[sampled_idx]

      # Fit model
      fit <- tryCatch({
        if (model_type == "MFA") {
          suppressMessages({
            mfa_em_fit(
              list_of_data = subsampled_data,
              K = K,
              r = r,
              max_iter = max_iter,
              n_init = n_init,
              use_kmeans_init = TRUE
            )
          })
        } else {
          suppressMessages({
            mixture_pca_em_fit(
              list_of_data = subsampled_data,
              K = K,
              r = r,
              max_iter = max_iter,
              n_init = n_init,
              use_kmeans_init = TRUE
            )
          })
        }
      }, error = function(e) NULL)

      if (!is.null(fit)) {
        # Get hard cluster assignments
        cluster_assignments <- apply(fit$resp, 1, which.max)
        assignments_matrix[b, sampled_idx] <- cluster_assignments
      }
    }

    # Compute instability as variance of co-clustering probabilities
    # For perfect clustering, all pairwise probabilities should be 0 or 1
    # Higher variance of intermediate values = more instability
    cocluster_probs <- numeric(0)

    for (i in 1:(N-1)) {
      for (j in (i+1):N) {
        both_sampled <- !is.na(assignments_matrix[, i]) & !is.na(assignments_matrix[, j])
        n_both <- sum(both_sampled)

        if (n_both >= 5) {  # Require minimum samples for reliable estimate
          prob_same <- mean(assignments_matrix[both_sampled, i] == assignments_matrix[both_sampled, j])
          cocluster_probs <- c(cocluster_probs, prob_same)
        }
      }
    }

    if (length(cocluster_probs) > 0) {
      # Stability = 1 - 4 * mean(p * (1-p))
      # This equals 1 when all probs are 0 or 1 (perfect stability)
      # This equals 0 when all probs are 0.5 (maximum instability)
      instability <- mean(cocluster_probs * (1 - cocluster_probs))
      stability_scores[ki] <- 1 - 4 * instability
    } else {
      stability_scores[ki] <- 0
    }

    if (verbose) {
      cat(sprintf("  K=%d: stability=%.4f (n_pairs=%d)\n", K, stability_scores[ki], length(cocluster_probs)))
    }
  }

  # Select K with highest stability
  best_idx <- which.max(stability_scores)
  best_K <- Kvec[best_idx]

  if (verbose) {
    cat(sprintf("\nSelected K=%d (stability=%.4f)\n", best_K, stability_scores[best_idx]))
  }

  stability_df <- data.frame(
    K = Kvec,
    stability = stability_scores
  )

  list(
    best_K = best_K,
    stability_scores = stability_scores,
    stability_df = stability_df
  )
}


#' Select K Using Cross-Validation
#'
#' This function implements K-fold cross-validation for choosing the optimal
#' number of clusters K. It evaluates model generalization by computing
#' log-likelihood on held-out subjects.
#'
#' @param list_of_data A list of data matrices, each (T_i x M).
#' @param Kvec Vector of candidate K values to test.
#' @param r Number of factors (fixed).
#' @param nfolds Number of cross-validation folds (default: 5).
#' @param max_iter Maximum EM iterations per fit.
#' @param n_init Number of random initializations per fit.
#' @param model_type Either "MFA" or "MPCA" (default: "MFA").
#' @param verbose Print progress messages (default: TRUE).
#'
#' @return A list containing:
#' \describe{
#'   \item{best_K}{The selected optimal K based on CV log-likelihood}
#'   \item{cv_scores}{Mean CV log-likelihood for each K}
#'   \item{cv_se}{Standard error of CV log-likelihood for each K}
#'   \item{cv_df}{Data frame with K, mean CV score, and SE}
#' }
#'
#' @details
#' The cross-validation procedure:
#' 1. Split subjects into nfolds groups
#' 2. For each K and each fold:
#'    - Train model on (nfolds-1) folds
#'    - Compute log-likelihood on held-out fold
#' 3. Average log-likelihoods across folds for each K
#' 4. Select K with highest average CV log-likelihood
#'
#' Higher CV log-likelihood indicates better generalization.
#'
#' @export
select_K_crossvalidation <- function(
    list_of_data,
    Kvec = 2:5,
    r = 2,
    nfolds = 5,
    max_iter = 50,
    n_init = 1,
    model_type = c("MFA", "MPCA"),
    verbose = TRUE
) {
  model_type <- match.arg(model_type)
  N <- length(list_of_data)
  M <- ncol(list_of_data[[1]])

  if (N < nfolds) {
    stop("Number of subjects (N) must be >= nfolds")
  }

  # Create fold assignments
  fold_ids <- rep(1:nfolds, length.out = N)
  fold_ids <- sample(fold_ids)  # Randomize

  if (verbose) {
    cat("=== Cross-Validation for K Selection ===\n")
    cat(sprintf("N subjects: %d, %d-fold CV\n", N, nfolds))
    cat(sprintf("Model type: %s, r fixed at: %d\n\n", model_type, r))
  }

  cv_scores <- matrix(NA, nrow = length(Kvec), ncol = nfolds)
  rownames(cv_scores) <- as.character(Kvec)

  for (ki in seq_along(Kvec)) {
    K <- Kvec[ki]
    if (verbose) cat(sprintf("Testing K=%d...\n", K))

    for (fold in 1:nfolds) {
      # Split data
      test_idx <- which(fold_ids == fold)
      train_idx <- which(fold_ids != fold)

      train_data <- list_of_data[train_idx]
      test_data <- list_of_data[test_idx]

      # Fit model on training data
      fit <- tryCatch({
        if (model_type == "MFA") {
          suppressMessages({
            mfa_em_fit(
              list_of_data = train_data,
              K = K,
              r = r,
              max_iter = max_iter,
              n_init = n_init,
              use_kmeans_init = TRUE
            )
          })
        } else {
          suppressMessages({
            mixture_pca_em_fit(
              list_of_data = train_data,
              K = K,
              r = r,
              max_iter = max_iter,
              n_init = n_init,
              use_kmeans_init = TRUE
            )
          })
        }
      }, error = function(e) NULL)

      if (!is.null(fit)) {
        # Compute log-likelihood on test data
        test_loglik <- compute_held_out_loglik(
          test_data = test_data,
          fit = fit,
          model_type = model_type
        )
        # Normalize by number of observations in test data
        n_test_obs <- sum(sapply(test_data, nrow))
        cv_scores[ki, fold] <- test_loglik / n_test_obs
      }
    }

    if (verbose) {
      mean_cv <- mean(cv_scores[ki, ], na.rm = TRUE)
      se_cv <- sd(cv_scores[ki, ], na.rm = TRUE) / sqrt(sum(!is.na(cv_scores[ki, ])))
      cat(sprintf("  K=%d: CV loglik=%.2f (SE=%.2f)\n", K, mean_cv, se_cv))
    }
  }

  # Compute mean and SE for each K
  cv_mean <- apply(cv_scores, 1, mean, na.rm = TRUE)
  cv_se <- apply(cv_scores, 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

  # Find best K using "one standard error rule"
  # Select smallest K whose score is within 1 SE of the best score
  best_score_idx <- which.max(cv_mean)
  best_score <- cv_mean[best_score_idx]
  best_score_se <- cv_se[best_score_idx]

  # Threshold: best score - 1 SE
  threshold <- best_score - best_score_se

  # Find smallest K above threshold
  above_threshold <- which(cv_mean >= threshold)
  best_K <- Kvec[min(above_threshold)]

  if (verbose) {
    cat(sprintf("\nBest score: K=%d (CV=%.2f, SE=%.2f)\n",
                Kvec[best_score_idx], best_score, best_score_se))
    cat(sprintf("1-SE threshold: %.2f\n", threshold))
    cat(sprintf("Selected K=%d (CV loglik=%.2f) using 1-SE rule\n",
                best_K, cv_mean[which(Kvec == best_K)]))
  }

  cv_df <- data.frame(
    K = Kvec,
    cv_loglik = cv_mean,
    cv_se = cv_se
  )

  list(
    best_K = best_K,
    cv_scores = cv_mean,
    cv_se = cv_se,
    cv_df = cv_df
  )
}


#' Compute Log-Likelihood on Held-Out Data
#'
#' Internal function to compute log-likelihood of held-out subjects
#' given a fitted mixture model.
#'
#' @param test_data List of data matrices for held-out subjects.
#' @param fit Fitted model object from mfa_em_fit or mpca_em_fit.
#' @param model_type Either "MFA" or "MPCA".
#'
#' @return Total log-likelihood of test data.
#'
#' @keywords internal
compute_held_out_loglik <- function(test_data, fit, model_type) {
  # Get mixing proportions (handle both pi and pi_k naming)
  pi_vec <- if (!is.null(fit$pi)) fit$pi else fit$pi_k
  K <- length(pi_vec)
  M <- ncol(test_data[[1]])
  N_test <- length(test_data)

  total_loglik <- 0

  for (i in seq_len(N_test)) {
    X_i <- test_data[[i]]
    T_i <- nrow(X_i)

    # Compute log-likelihood for each cluster
    log_probs <- numeric(K)

    for (k in 1:K) {
      if (model_type == "MFA") {
        # MFA: X ~ N(mu_k, Lambda_k %*% t(Lambda_k) + Psi_k)
        mu_k <- as.vector(fit$mu[[k]])
        Lambda_k <- fit$Lambda[[k]]
        Psi_k <- fit$Psi[[k]]

        # Covariance matrix
        # Psi_k is already a diagonal matrix (MxM)
        Sigma_k <- Lambda_k %*% t(Lambda_k) + Psi_k
      } else {
        # MPCA: X ~ N(mu_k, W_k %*% t(W_k) + sigma2_k * I)
        mu_k <- as.vector(fit$mu[[k]])
        W_k <- fit$W[[k]]
        sigma2_k <- fit$sigma2[k]

        # Covariance matrix
        Sigma_k <- W_k %*% t(W_k) + sigma2_k * diag(M)
      }

      # Compute log-likelihood for all rows of subject i
      # Using Cholesky decomposition for numerical stability
      tryCatch({
        L <- chol(Sigma_k)
        log_det <- 2 * sum(log(diag(L)))

        ll_k <- 0
        for (t in 1:T_i) {
          x_centered <- X_i[t, ] - mu_k
          quad_form <- sum(backsolve(L, x_centered, transpose = TRUE)^2)
          ll_k <- ll_k + (-0.5 * (M * log(2 * pi) + log_det + quad_form))
        }

        log_probs[k] <- log(pi_vec[k]) + ll_k
      }, error = function(e) {
        log_probs[k] <<- -Inf
      })
    }

    # Log-sum-exp for numerical stability
    max_log_prob <- max(log_probs)
    if (is.finite(max_log_prob)) {
      total_loglik <- total_loglik + max_log_prob + log(sum(exp(log_probs - max_log_prob)))
    }
  }

  total_loglik
}


#' Select Minimum r Based on VAF Threshold for MFA
#'
#' This function implements a traditional EMG synergy approach: select the minimum
#' r (number of factors) that achieves a target VAF (Variance Accounted For) threshold.
#' This is commonly used in motor control research where VAF >= 90-95 percent is desired.
#'
#' @param list_of_data A list of data matrices, each (T_i x M).
#' @param K Number of clusters (must be pre-specified).
#' @param rvec Vector of candidate r values to test.
#' @param vaf_threshold Target VAF threshold (default: 0.90 for 90 percent).
#' @param max_iter Maximum EM iterations.
#' @param nIterFA Sub-iterations within the factor-analyzer update.
#' @param tol Convergence tolerance.
#' @param n_init Number of random initializations.
#' @param use_kmeans_init Whether to use k-means initialization.
#' @param subject_rdim_for_kmeans Dimension for PCA-based features.
#' @param mc_cores Number of CPU cores for parallelization.
#' @param n_threads Number of OpenMP threads.
#'
#' @return A list containing: summary (data frame with r, VAF, BIC for each r),
#'   best_r_vaf (minimum r achieving VAF threshold), best_r_bic (optimal r by BIC),
#'   best_model_vaf (model with minimum r achieving threshold), and
#'   best_model_bic (model with optimal BIC).
#'
#' @details
#' This function fits MFA models for each r in rvec with fixed K, computes both
#' VAF and BIC, and returns:
#' 1. The minimum r that achieves the VAF threshold (traditional approach)
#' 2. The optimal r by BIC (statistical approach)
#' 
#' Users can compare both criteria to make informed model selection decisions.
#'
#' @export
select_r_by_vaf_mfa <- function(
    list_of_data,
    K,
    rvec = 1:5,
    vaf_threshold = 0.90,
    max_iter = 50,
    nIterFA = 5,
    tol = 1e-3,
    n_init = 1,
    use_kmeans_init = TRUE,
    subject_rdim_for_kmeans = 2,
    mc_cores = 1,
    n_threads = 1
){
  N <- length(list_of_data)
  M <- ncol(list_of_data[[1]])
  N_total_rows <- sum(sapply(list_of_data, nrow))
  
  cat(sprintf("=== VAF-Based r Selection for MFA (K=%d) ===\n", K))
  cat(sprintf("Target VAF threshold: %.1f%%\n", vaf_threshold * 100))
  
  results <- list()
  for(i in seq_along(rvec)) {
    r <- rvec[i]
    cat(sprintf("  Fitting r=%d...\n", r))
    
    fit_mfa <- mfa_em_fit(
      list_of_data = list_of_data,
      K = K,
      r = r,
      max_iter = max_iter,
      nIterFA = nIterFA,
      tol = tol,
      n_init = n_init,
      use_kmeans_init = use_kmeans_init,
      subject_rdim_for_kmeans = subject_rdim_for_kmeans,
      mc_cores = mc_cores,
      n_threads = n_threads
    )
    
    loglik_val <- fit_mfa$logLik
    bic_val <- compute_BIC_mfa(loglik_val, K, r, M, N_total_rows)
    vaf_val <- compute_global_vaf_mfa(list_of_data, fit_mfa)
    
    results[[i]] <- list(
      r = r,
      logLik = loglik_val,
      BIC = bic_val,
      VAF = vaf_val,
      model = fit_mfa
    )
    
    cat(sprintf("    VAF = %.3f, BIC = %.2f\n", vaf_val, bic_val))
  }
  
  # Build summary
  df_summary <- do.call(rbind, lapply(results, function(res) {
    data.frame(r=res$r, VAF=res$VAF, BIC=res$BIC, logLik=res$logLik)
  }))
  
  # Find minimum r that achieves VAF threshold
  meets_threshold <- df_summary$VAF >= vaf_threshold
  if(any(meets_threshold)) {
    best_r_vaf <- min(df_summary$r[meets_threshold])
    cat(sprintf("\n* Minimum r achieving VAF >= %.1f%%: r = %d (VAF = %.3f)\n", 
                vaf_threshold * 100, best_r_vaf, 
                df_summary$VAF[df_summary$r == best_r_vaf]))
  } else {
    best_r_vaf <- max(rvec)
    cat(sprintf("\nX No r achieved VAF >= %.1f%%. Best VAF = %.3f at r = %d\n",
                vaf_threshold * 100, max(df_summary$VAF), best_r_vaf))
  }
  
  # Find best r by BIC
  df_summary_bic <- df_summary[order(df_summary$BIC), ]
  best_r_bic <- df_summary_bic$r[1]
  cat(sprintf("* Optimal r by BIC: r = %d (BIC = %.2f, VAF = %.3f)\n",
              best_r_bic, df_summary_bic$BIC[1], df_summary_bic$VAF[1]))
  
  # Get models
  best_model_vaf <- results[[which(sapply(results, function(x) x$r) == best_r_vaf)[1]]]$model
  best_model_bic <- results[[which(sapply(results, function(x) x$r) == best_r_bic)[1]]]$model
  
  list(
    summary = df_summary,
    best_r_vaf = best_r_vaf,
    best_r_bic = best_r_bic,
    best_model_vaf = best_model_vaf,
    best_model_bic = best_model_bic
  )
}

#' Select Minimum r Based on VAF Threshold for MPCA
#'
#' This function implements a traditional EMG synergy approach: select the minimum
#' r (number of principal components) that achieves a target VAF threshold.
#' This is commonly used in motor control research where VAF >= 90-95 percent is desired.
#'
#' @param list_of_data A list of data matrices, each (T_i x M).
#' @param K Number of clusters (must be pre-specified).
#' @param rvec Vector of candidate r values to test.
#' @param vaf_threshold Target VAF threshold (default: 0.90 for 90 percent).
#' @param max_iter Maximum EM iterations.
#' @param nIterPCA Sub-iterations within the PCA/PPCA update.
#' @param tol Convergence tolerance.
#' @param method "EM" or "closed_form".
#' @param n_init Number of random initializations.
#' @param use_kmeans_init Whether to use k-means initialization.
#' @param subject_rdim_for_kmeans Dimension for PCA-based features.
#' @param mc_cores Number of CPU cores for parallelization.
#' @param n_threads Number of OpenMP threads.
#'
#' @return A list containing: summary (data frame with r, VAF, BIC for each r),
#'   best_r_vaf (minimum r achieving VAF threshold), best_r_bic (optimal r by BIC),
#'   best_model_vaf (model with minimum r achieving threshold), and
#'   best_model_bic (model with optimal BIC).
#'
#' @details
#' This function fits MPCA models for each r in rvec with fixed K, computes both
#' VAF and BIC, and returns:
#' 1. The minimum r that achieves the VAF threshold (traditional approach)
#' 2. The optimal r by BIC (statistical approach)
#' 
#' Users can compare both criteria to make informed model selection decisions.
#'
#' @export
select_r_by_vaf_mpca <- function(
    list_of_data,
    K,
    rvec = 1:5,
    vaf_threshold = 0.90,
    max_iter = 50,
    nIterPCA = 5,
    tol = 1e-3,
    method = "EM",
    n_init = 1,
    use_kmeans_init = TRUE,
    subject_rdim_for_kmeans = 2,
    mc_cores = 1,
    n_threads = 1
){
  N <- length(list_of_data)
  M <- ncol(list_of_data[[1]])
  N_total_rows <- sum(sapply(list_of_data, nrow))
  
  cat(sprintf("=== VAF-Based r Selection for MPCA (K=%d) ===\n", K))
  cat(sprintf("Target VAF threshold: %.1f%%\n", vaf_threshold * 100))
  
  results <- list()
  for(i in seq_along(rvec)) {
    r <- rvec[i]
    cat(sprintf("  Fitting r=%d...\n", r))
    
    fit_pca <- mixture_pca_em_fit(
      list_of_data = list_of_data,
      K = K,
      r = r,
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
    
    loglik_val <- fit_pca$logLik
    bic_val <- compute_BIC_mpca(loglik_val, K, r, M, N_total_rows)
    vaf_val <- compute_global_vaf_mpca(list_of_data, fit_pca)
    
    results[[i]] <- list(
      r = r,
      logLik = loglik_val,
      BIC = bic_val,
      VAF = vaf_val,
      model = fit_pca
    )
    
    cat(sprintf("    VAF = %.3f, BIC = %.2f\n", vaf_val, bic_val))
  }
  
  # Build summary
  df_summary <- do.call(rbind, lapply(results, function(res) {
    data.frame(r=res$r, VAF=res$VAF, BIC=res$BIC, logLik=res$logLik)
  }))
  
  # Find minimum r that achieves VAF threshold
  meets_threshold <- df_summary$VAF >= vaf_threshold
  if(any(meets_threshold)) {
    best_r_vaf <- min(df_summary$r[meets_threshold])
    cat(sprintf("\n* Minimum r achieving VAF >= %.1f%%: r = %d (VAF = %.3f)\n", 
                vaf_threshold * 100, best_r_vaf, 
                df_summary$VAF[df_summary$r == best_r_vaf]))
  } else {
    best_r_vaf <- max(rvec)
    cat(sprintf("\nX No r achieved VAF >= %.1f%%. Best VAF = %.3f at r = %d\n",
                vaf_threshold * 100, max(df_summary$VAF), best_r_vaf))
  }
  
  # Find best r by BIC
  df_summary_bic <- df_summary[order(df_summary$BIC), ]
  best_r_bic <- df_summary_bic$r[1]
  cat(sprintf("* Optimal r by BIC: r = %d (BIC = %.2f, VAF = %.3f)\n",
              best_r_bic, df_summary_bic$BIC[1], df_summary_bic$VAF[1]))
  
  # Get models
  best_model_vaf <- results[[which(sapply(results, function(x) x$r) == best_r_vaf)[1]]]$model
  best_model_bic <- results[[which(sapply(results, function(x) x$r) == best_r_bic)[1]]]$model
  
  list(
    summary = df_summary,
    best_r_vaf = best_r_vaf,
    best_r_bic = best_r_bic,
    best_model_vaf = best_model_vaf,
    best_model_bic = best_model_bic
  )
}



# ============================================================
# Unified Metric Functions (Phase 1.2 API Unification)
# ============================================================

#' Compute Global Variance Accounted For (VAF)
#'
#' Unified interface for computing global VAF for both MFA and MPCA models.
#' Automatically detects model type from the fit object structure.
#'
#' @param list_of_data A list of N matrices, each of dimension (T_i x M)
#' @param fit A fitted model object (MFA or MPCA)
#' @param model_type Character, either "MFA" or "MPCA". If NULL, auto-detected
#'   from the fit object structure.
#'
#' @return A numeric scalar representing the global VAF (0 to 1)
#'
#' @details
#' VAF is computed as 1 - SSE/SST, where SSE is the sum of squared errors
#' between observed and reconstructed data, and SST is the total sum of squares.
#'
#' @seealso \code{\link{compute_global_vaf_mfa}}, \code{\link{compute_global_vaf_mpca}}
#'
#' @examples
#' \dontrun{
#' # Works with either MFA or MPCA fit objects
#' vaf <- compute_VAF(list_of_data, fit)
#'
#' # Or explicitly specify model type
#' vaf <- compute_VAF(list_of_data, fit, model_type = "MFA")
#' }
#'
#' @export
compute_VAF <- function(list_of_data, fit, model_type = NULL) {
  # Auto-detect model type if not specified
  if (is.null(model_type)) {
    model_type <- validate_fit_object(fit)
  } else {
    model_type <- validate_model_type(model_type)
  }

  if (model_type == "MFA") {
    compute_global_vaf_mfa(list_of_data, fit)
  } else {
    compute_global_vaf_mpca(list_of_data, fit)
  }
}


#' Compute Log-Likelihood for Mixture Models
#'
#' Unified interface for computing log-likelihood for both MFA and MPCA models.
#' Automatically detects model type from the fit object structure.
#'
#' @param list_of_data A list of N matrices, each of dimension (T_i x M)
#' @param fit A fitted model object (MFA or MPCA)
#' @param model_type Character, either "MFA" or "MPCA". If NULL, auto-detected
#'   from the fit object structure.
#'
#' @return A numeric scalar representing the total log-likelihood
#'
#' @seealso \code{\link{compute_logLik_mfa}}, \code{\link{compute_logLik_mpca}}
#'
#' @examples
#' \dontrun{
#' # Works with either MFA or MPCA fit objects
#' ll <- compute_logLik(list_of_data, fit)
#'
#' # Or explicitly specify model type
#' ll <- compute_logLik(list_of_data, fit, model_type = "MFA")
#' }
#'
#' @export
compute_logLik <- function(list_of_data, fit, model_type = NULL) {
  # Auto-detect model type if not specified
  if (is.null(model_type)) {
    model_type <- validate_fit_object(fit)
  } else {
    model_type <- validate_model_type(model_type)
  }

  if (model_type == "MFA") {
    compute_logLik_mfa(list_of_data, fit)
  } else {
    compute_logLik_mpca(list_of_data, fit)
  }
}


#' Compute Bayesian Information Criterion (BIC)
#'
#' Unified interface for computing BIC for both MFA and MPCA models.
#'
#' @param loglik Numeric scalar, the log-likelihood value
#' @param K Integer, number of clusters
#' @param r Integer, number of factors/synergies
#' @param M Integer, number of variables (muscles)
#' @param N_obs Integer, number of observations (total time points across subjects)
#' @param model_type Character, either "MFA" or "MPCA"
#'
#' @return A numeric scalar representing the BIC value (lower is better)
#'
#' @details
#' BIC = -2 * loglik + p * log(N_obs), where p is the number of parameters.
#' The parameter count differs between MFA and MPCA due to different model
#' structures.
#'
#' @seealso \code{\link{compute_BIC_mfa}}, \code{\link{compute_BIC_mpca}}
#'
#' @examples
#' \dontrun{
#' # Compute BIC for an MFA model
#' bic <- compute_BIC(loglik, K = 2, r = 3, M = 8, N_obs = 1000, model_type = "MFA")
#' }
#'
#' @export
compute_BIC <- function(loglik, K, r, M, N_obs, model_type) {
  model_type <- validate_model_type(model_type)

  if (model_type == "MFA") {
    compute_BIC_mfa(loglik, K, r, M, N_obs)
  } else {
    compute_BIC_mpca(loglik, K, r, M, N_obs)
  }
}


#' Compute Reconstruction Error
#'
#' Unified interface for computing reconstruction error for both MFA and MPCA models.
#' Automatically detects model type from the fit object structure.
#'
#' @param list_of_data A list of N matrices, each of dimension (T_i x M)
#' @param fit A fitted model object (MFA or MPCA)
#' @param model_type Character, either "MFA" or "MPCA". If NULL, auto-detected
#'   from the fit object structure.
#' @param metric Character, one of "sse" (sum of squared errors) or
#'   "rmse" (root mean squared error). Default: "sse"
#'
#' @return A numeric scalar representing the reconstruction error
#'
#' @examples
#' \dontrun{
#' # Compute SSE
#' sse <- compute_reconstruction_error(list_of_data, fit)
#'
#' # Compute RMSE
#' rmse <- compute_reconstruction_error(list_of_data, fit, metric = "rmse")
#' }
#'
#' @export
compute_reconstruction_error <- function(list_of_data, fit,
                                          model_type = NULL,
                                          metric = c("sse", "rmse")) {
  metric <- match.arg(metric)

  # Compute VAF first, then derive error
  vaf <- compute_VAF(list_of_data, fit, model_type)

  # Compute total sum of squares
  total_obs <- sum(vapply(list_of_data, function(x) sum(x^2), numeric(1)))
  sst <- total_obs

  # SSE = SST * (1 - VAF)
  sse <- sst * (1 - vaf)

  if (metric == "sse") {
    return(sse)
  } else {
    # RMSE = sqrt(SSE / n)
    n_total <- sum(vapply(list_of_data, function(x) length(x), integer(1)))
    return(sqrt(sse / n_total))
  }
}
