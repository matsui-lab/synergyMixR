#' Extract subject-level features by single PCA
#'
#' This function applies PCA to the concatenated data across all subjects,
#' then for each subject \code{i}, computes the mean of the top \code{r_dim} principal
#' component scores as a single feature vector. This is useful for initializing
#' cluster assignments in mixture models.
#'
#' @param list_of_data A list of \code{(T_i x M)} matrices, where each element
#'   represents time-series data for one subject. \code{T_i} is the number of
#'   time points for subject \code{i}, and \code{M} is the number of variables
#'   (e.g., muscles in EMG data).
#' @param r_dim Number of principal components to extract. Defaults to 2.
#'
#' @return A numeric matrix of shape \code{(N x r_dim)}, where \code{N=length(list_of_data)}.
#'   Each row contains the mean PC scores for one subject.
#'
#' @details
#' The function first concatenates all subject data into a single matrix, performs
#' PCA on this combined dataset, and then computes the mean of the PC scores for
#' each subject's time points. This provides a low-dimensional representation of
#' each subject that can be used for clustering initialization.
#'
#' @examples
#' \dontrun{
#' # Create example data
#' N <- 10
#' M <- 8
#' list_of_data <- lapply(1:N, function(i) matrix(rnorm(100 * M), 100, M))
#'
#' # Extract 2D features for each subject
#' features <- extract_subject_features_by_singlePCA(list_of_data, r_dim = 2)
#' dim(features)  # Should be 10 x 2
#' }
#'
#' @export
extract_subject_features_by_singlePCA <- function(list_of_data, r_dim = 2)
{
  # Concatenate all subject data into one matrix
  bigX <- do.call(rbind, list_of_data)

  # Perform PCA on the combined data
  pca_res <- prcomp(bigX, center = TRUE, scale. = FALSE)
  scores_all <- pca_res$x[, seq_len(r_dim), drop = FALSE]

  # Compute mean PC scores for each subject
  N <- length(list_of_data)
  out_features <- matrix(NA, nrow = N, ncol = r_dim)

  start_idx <- 1
  for (i in seq_len(N)) {
    Ti <- nrow(list_of_data[[i]])
    end_idx <- start_idx + Ti - 1
    out_features[i, ] <- colMeans(scores_all[start_idx:end_idx, , drop = FALSE])
    start_idx <- end_idx + 1
  }

  out_features
}

#' K-means-based cluster assignment
#'
#' Performs k-means clustering on subject-level features to obtain initial
#' cluster assignments for mixture model fitting.
#'
#' @param features An \code{(N x d)} matrix of subject-level features, where
#'   \code{N} is the number of subjects and \code{d} is the feature dimension.
#' @param K Number of clusters to create.
#'
#' @return An integer vector of length \code{N} containing cluster labels (1..K).
#'   Each element indicates which cluster the corresponding subject is assigned to.
#'
#' @details
#' This function uses the standard k-means algorithm with 10 random starts
#' (\code{nstart=10}) to find a stable clustering solution. The resulting
#' cluster assignments can be used as initial values for EM algorithms in
#' mixture models.
#'
#' @examples
#' \dontrun{
#' # Create example features
#' features <- matrix(rnorm(100), nrow = 10, ncol = 10)
#'
#' # Assign to 3 clusters
#' z_init <- assign_by_kmeans(features, K = 3)
#' table(z_init)  # Check cluster sizes
#' }
#'
#' @export
assign_by_kmeans <- function(features, K) {
  km_res <- kmeans(features, centers = K, nstart = 10)
  km_res$cluster
}
