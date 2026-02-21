#' Reconstruct Data from Single Factor Analysis Model (Internal)
#'
#' Given a single subject's data and factor analysis parameters, this function
#' reconstructs the data using the factor scores computed via the posterior mean.
#' **Internal function**.
#'
#' @param X Original \code{(T x M)} data matrix (used for dimension reference)
#' @param mu A length-M mean vector
#' @param Lambda An \code{(M x r)} factor loading matrix
#' @param Psi_diag A length-M vector of diagonal elements from Psi matrix
#'
#' @return A \code{(T x M)} reconstructed matrix Xhat
#'
#' @details
#' This function first computes factor scores using \code{\link{compute_factor_scores}},
#' then reconstructs the data as Xhat = Z * Lambda^T + mu, where Z are the factor scores.
#'
#' @examples
#' \dontrun{
#' # Suppose we have X (T x M), mu (M), Lambda (M x r), Psi_diag (M)
#' Xhat <- reconstruct_singleFA(X, mu, Lambda, Psi_diag)
#' }
#'
#' @noRd
reconstruct_singleFA <- function(X, mu, Lambda, Psi_diag) {
  Psi <- diag(Psi_diag)
  scores <- compute_factor_scores(X, mu, Lambda, Psi)
  Xhat   <- scores %*% t(Lambda)
  Xhat   <- sweep(Xhat, 2, mu, "+")
  Xhat
}
