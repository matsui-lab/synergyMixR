#' Align Basis Vectors Across Subjects Using Hungarian Algorithm
#'
#' This function aligns the columns (basis vectors) of each subject's basis matrix
#' to a reference basis using the Hungarian algorithm for optimal assignment.
#' This ensures that corresponding components across subjects are in the same order,
#' which improves clustering accuracy in two-step methods.
#'
#' @param W_list A list of basis matrices, one per subject. Each matrix should be
#'   M x r (M channels, r components). These are typically the output from individual
#'   PCA or FA decompositions in Step 1 of the two-step baseline method.
#' @param reference Integer index of the reference subject (default: 1). The basis
#'   of this subject will be used as the reference for alignment. All other subjects'
#'   basis vectors will be reordered to match this reference.
#' @param method Character string specifying the similarity metric. Currently only
#'   "correlation" is supported (Pearson correlation between basis vectors).
#'   Future versions may support "cosine" similarity.
#' @param verbose Logical. If \code{TRUE}, prints informational messages during alignment.
#'   Default is \code{FALSE} to avoid spamming output during parameter sweeps.
#'
#' @return A list of aligned basis matrices with the same structure as \code{W_list},
#'   but with columns reordered to match the reference basis. The reference subject's
#'   basis remains unchanged. Each element is an M x r matrix where columns have been
#'   permuted to maximize correspondence with the reference.
#'
#' @details
#' ## When This Function Is Called
#'
#' This function is automatically called in \strong{Step 1.5} of the two-step baseline
#' method (\code{\link{fit_two_step_baseline}}), between individual decomposition
#' (Step 1) and clustering (Step 2). It is called regardless of the
#' \code{refine_assignments} parameter setting, as alignment improves clustering
#' accuracy in both modes.
#'
#' ## Why Alignment Is Necessary
#'
#' Individual PCA/FA produces basis vectors in arbitrary order. Without alignment,
#' the first component of subject A may represent a different pattern than the first
#' component of subject B. This arbitrary ordering causes k-means clustering in Step 2
#' to group unrelated components together, leading to poor cluster assignments and
#' low ARI (Adjusted Rand Index).
#'
#' ## How Alignment Works
#'
#' The alignment process works as follows:
#' \enumerate{
#'   \item Select a reference basis (default: first subject).
#'   \item For each subject, compute a similarity matrix between the reference
#'         basis columns and the subject's basis columns using Pearson correlation.
#'   \item Use the Hungarian algorithm (\code{clue::solve_LSAP}) to find the
#'         optimal column permutation that maximizes total similarity.
#'   \item Reorder the subject's basis columns according to this permutation.
#' }
#'
#' The Hungarian algorithm solves the linear sum assignment problem, finding the
#' permutation that maximizes the sum of similarities. The correlation matrix is
#' transformed from [-1, 1] to [0, 1] to satisfy \code{solve_LSAP}'s requirement
#' for nonnegative entries. This affine transformation preserves the optimal matching.
#' After permutation, column signs are corrected to ensure positive correlation with
#' the reference, preventing k-means from splitting components by arbitrary sign.
#'
#' ## Relationship with refine_assignments
#'
#' This alignment step is independent of the \code{refine_assignments} parameter
#' in \code{fit_two_step_baseline()}:
#' \itemize{
#'   \item \strong{refine_assignments = FALSE}: Alignment is applied, then k-means
#'         clustering is performed on aligned basis vectors. No reassignment after
#'         cluster-specific model fitting.
#'   \item \strong{refine_assignments = TRUE}: Alignment is applied, k-means clustering
#'         is performed, cluster-specific models are fitted, and then subjects are
#'         reassigned based on reconstruction error.
#' }
#'
#' In both cases, alignment ensures that corresponding components are in the same
#' order before clustering, which significantly improves clustering accuracy.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic usage with simulated data
#' library(synergyMixR)
#' 
#' # Simulate some basis matrices (8 channels, 2 components, 3 subjects)
#' set.seed(123)
#' W_list <- list(
#'   matrix(rnorm(8*2), 8, 2),
#'   matrix(rnorm(8*2), 8, 2),
#'   matrix(rnorm(8*2), 8, 2)
#' )
#' 
#' # Align to first subject
#' W_aligned <- align_basis_across_subjects(W_list, reference = 1)
#' 
#' # Check that dimensions are preserved
#' sapply(W_aligned, dim)  # Each should be 8 x 2
#' 
#' # Example 2: Usage within two-step baseline
#' # This is done automatically in fit_two_step_baseline()
#' sim_data <- generate_synergy_data(N = 50, K = 2, r = 2, M = 6, T_each = 100, seed = 123)
#' 
#' # Extract individual basis vectors (Step 1)
#' W_list_individual <- lapply(sim_data$list_of_data, function(X) {
#'   pca_res <- prcomp(X, center = TRUE, scale. = FALSE)
#'   pca_res$rotation[, 1:2]  # First 2 components
#' })
#' 
#' # Align basis vectors (Step 1.5) - this happens automatically in fit_two_step_baseline
#' W_list_aligned <- align_basis_across_subjects(W_list_individual, reference = 1)
#' 
#' # Now W_list_aligned can be used for clustering (Step 2)
#' # with improved accuracy due to aligned components
#' }
#'
#' @seealso \code{\link{fit_two_step_baseline}} for the two-step baseline method
#'   that uses this alignment function.
align_basis_across_subjects <- function(W_list, reference = 1, method = "correlation", verbose = FALSE) {
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("Package 'clue' is required for basis alignment. Please install it with: install.packages('clue')")
  }
  
  if (length(W_list) == 0) {
    stop("W_list cannot be empty")
  }
  
  if (reference < 1 || reference > length(W_list)) {
    stop(sprintf("reference must be between 1 and %d (length of W_list)", length(W_list)))
  }
  
  # Get reference basis
  ref <- W_list[[reference]]
  r <- ncol(ref)
  
  # Initialize aligned list
  aligned <- vector("list", length(W_list))
  
  for (i in seq_along(W_list)) {
    W <- W_list[[i]]
    
    # Check dimensions
    if (ncol(W) != r) {
      stop(sprintf("Subject %d has %d components, but reference has %d", i, ncol(W), r))
    }
    
    if (i == reference) {
      # Reference subject: no alignment needed
      aligned[[i]] <- W
    } else {
      # Compute similarity matrix
      if (method == "correlation") {
        # Correlation between reference columns and subject columns
        # sim[j, k] = correlation between ref[, j] and W[, k]
        sim <- cor(ref, W, use = "pairwise.complete.obs")
        
        # Handle NA/Inf values (replace with -1, the minimum correlation)
        sim[!is.finite(sim)] <- -1
        
        # Transform correlation from [-1, 1] to [0, 1] for solve_LSAP
        # This is required because solve_LSAP needs nonnegative entries
        # The transformation preserves optimal matching (affine transform)
        sim_nonneg <- (sim + 1) / 2
        
        if (verbose) {
          message(sprintf("Subject %d: Transformed correlation matrix to [0,1] range for Hungarian matching", i))
        }
      } else {
        stop(sprintf("Method '%s' is not supported. Currently only 'correlation' is available.", method))
      }
      
      # Hungarian algorithm: maximize similarity (use maximum = TRUE)
      # assignment[j] = k means ref column j matches W column k
      assignment <- clue::solve_LSAP(sim_nonneg, maximum = TRUE)
      
      # Reorder columns of W according to assignment
      aligned[[i]] <- W[, assignment, drop = FALSE]
      
      # Sign correction: flip columns so they have positive correlation with reference
      # This prevents k-means from splitting components by arbitrary sign
      for (j in seq_len(r)) {
        corr_sign <- sign(cor(ref[, j], aligned[[i]][, j], use = "pairwise.complete.obs"))
        if (is.na(corr_sign) || corr_sign == 0) {
          corr_sign <- 1  # Default to no flip if correlation is NA or zero
        }
        if (corr_sign < 0) {
          aligned[[i]][, j] <- -aligned[[i]][, j]
          if (verbose) {
            message(sprintf("Subject %d, component %d: Flipped sign to match reference", i, j))
          }
        }
      }
    }
  }
  
  aligned
}
