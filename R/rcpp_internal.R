#' Internal C++ Computational Functions
#'
#' Low-level computational functions implemented in C++ via Rcpp for internal package use.
#' These functions are automatically exported by Rcpp but are not intended for direct user access.
#'
#' @details 
#' **Note:** These functions are exported for internal use via Rcpp. Users should not call them directly.
#' 
#' \describe{
#'   \item{dmvnormFA_rowwiseCpp, dmvnormPCA_rowwiseCpp}{Multivariate normal density calculations for FA and PCA models}
#'   \item{faEMupdateCpp, pcaEMupdateCpp}{EM algorithm update steps for factor analysis and PCA}
#'   \item{mfaTimeseriesCpp, mpcaTimeseriesCpp}{Time series processing functions for mixture models}
#'   \item{pcaClosedFormCpp}{Closed-form PCA computation for efficiency}
#' }
#'
#' @importFrom Rcpp evalCpp
#' @useDynLib synergyMixR, .registration = TRUE
#' @importFrom stats kmeans prcomp sd cor rnorm runif
#' @importFrom utils txtProgressBar setTxtProgressBar head
#' @noRd
#' @name cpp_internal_functions
#' @aliases dmvnormFA_rowwiseCpp dmvnormPCA_rowwiseCpp faEMupdateCpp mfaTimeseriesCpp mpcaTimeseriesCpp pcaClosedFormCpp pcaEMupdateCpp
NULL
