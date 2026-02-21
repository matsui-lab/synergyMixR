# Tests for subspace recovery metrics (Step 4)
# These functions evaluate how well estimated subspaces recover true subspaces

# =============================================================================
# Tests for compute_subspace_similarity
# =============================================================================

test_that("compute_subspace_similarity returns ~1 for identical subspaces", {
  # Create a simple orthonormal basis
  L_true <- matrix(c(1, 0, 0, 0,
                     0, 1, 0, 0), nrow = 4, ncol = 2)
  L_hat <- L_true
  
  sim <- compute_subspace_similarity(L_true, L_hat)
  
  expect_true(is.numeric(sim))
  expect_true(!is.na(sim))
  expect_equal(sim, 1, tolerance = 1e-10)
})

test_that("compute_subspace_similarity returns ~1 for rotated subspaces", {
  # Create a basis and rotate it
  L_true <- matrix(c(1, 0, 0, 0,
                     0, 1, 0, 0), nrow = 4, ncol = 2)
  
  # 45-degree rotation matrix
  theta <- pi / 4
  R <- matrix(c(cos(theta), -sin(theta),
                sin(theta), cos(theta)), nrow = 2, ncol = 2)
  
  L_hat <- L_true %*% R
  
  sim <- compute_subspace_similarity(L_true, L_hat)
  
  expect_true(is.numeric(sim))
  expect_true(!is.na(sim))
  expect_equal(sim, 1, tolerance = 1e-6)
})

test_that("compute_subspace_similarity returns ~1 for scaled subspaces", {
  # Create a basis and scale it
  L_true <- matrix(c(1, 0, 0, 0,
                     0, 1, 0, 0), nrow = 4, ncol = 2)
  
  L_hat <- L_true * 5  # Scale by factor of 5
  
  sim <- compute_subspace_similarity(L_true, L_hat)
  
  expect_true(is.numeric(sim))
  expect_true(!is.na(sim))
  expect_equal(sim, 1, tolerance = 1e-6)
})

test_that("compute_subspace_similarity handles rank differences", {
  # True basis has rank 2, estimated has rank 3
  L_true <- matrix(c(1, 0, 0, 0,
                     0, 1, 0, 0), nrow = 4, ncol = 2)
  
  L_hat <- matrix(c(1, 0, 0, 0,
                    0, 1, 0, 0,
                    0, 0, 1, 0), nrow = 4, ncol = 3)
  
  sim <- compute_subspace_similarity(L_true, L_hat)
  
  expect_true(is.numeric(sim))
  expect_true(!is.na(sim))
  # Should use min_dim = 2, and the first 2 dimensions should match
  expect_equal(sim, 1, tolerance = 1e-6)
})

test_that("compute_subspace_similarity returns low value for orthogonal subspaces", {
  # Two orthogonal 2D subspaces in 4D space
  L_true <- matrix(c(1, 0, 0, 0,
                     0, 1, 0, 0), nrow = 4, ncol = 2)
  
  L_hat <- matrix(c(0, 0, 1, 0,
                    0, 0, 0, 1), nrow = 4, ncol = 2)
  
  sim <- compute_subspace_similarity(L_true, L_hat)
  
  expect_true(is.numeric(sim))
  expect_true(!is.na(sim))
  expect_true(sim < 0.1)  # Should be close to 0
})

test_that("compute_subspace_similarity returns NA for NULL inputs", {
  L_true <- matrix(1:4, 2, 2)
  
  expect_true(is.na(compute_subspace_similarity(NULL, L_true)))
  expect_true(is.na(compute_subspace_similarity(L_true, NULL)))
  expect_true(is.na(compute_subspace_similarity(NULL, NULL)))
})

test_that("compute_subspace_similarity returns NA for empty matrices", {
  L_true <- matrix(1:4, 2, 2)
  L_empty <- matrix(nrow = 0, ncol = 0)
  
  expect_true(is.na(compute_subspace_similarity(L_empty, L_true)))
  expect_true(is.na(compute_subspace_similarity(L_true, L_empty)))
})

test_that("compute_subspace_similarity returns NA for matrices with NA/NaN/Inf", {
  L_true <- matrix(1:4, 2, 2)
  L_na <- matrix(c(1, NA, 3, 4), 2, 2)
  L_nan <- matrix(c(1, NaN, 3, 4), 2, 2)
  L_inf <- matrix(c(1, Inf, 3, 4), 2, 2)
  
  expect_true(is.na(compute_subspace_similarity(L_na, L_true)))
  expect_true(is.na(compute_subspace_similarity(L_nan, L_true)))
  expect_true(is.na(compute_subspace_similarity(L_inf, L_true)))
})

test_that("compute_subspace_similarity returns NA for dimension mismatch", {
  L_true <- matrix(1:6, 3, 2)  # 3 rows

  L_hat <- matrix(1:8, 4, 2)   # 4 rows
  
  expect_warning(
    result <- compute_subspace_similarity(L_true, L_hat),
    "Dimension mismatch"
  )
  expect_true(is.na(result))
})

# =============================================================================
# Tests for align_clusters_to_truth
# =============================================================================

test_that("align_clusters_to_truth handles perfect permutation", {
  z_true <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
  z_hat <- c(3, 3, 3, 1, 1, 1, 2, 2, 2)  # Labels are permuted
  
  result <- align_clusters_to_truth(z_true, z_hat)
  
  expect_true(is.list(result))
  expect_true(!is.na(result$accuracy))
  expect_equal(result$accuracy, 1.0)  # Perfect after alignment
  expect_equal(result$z_aligned, z_true)
})

test_that("align_clusters_to_truth handles identical labels", {
  z_true <- c(1, 1, 2, 2, 3, 3)
  z_hat <- z_true
  
  result <- align_clusters_to_truth(z_true, z_hat)
  
  expect_true(is.list(result))
  expect_equal(result$accuracy, 1.0)
  expect_equal(result$z_aligned, z_true)
})

test_that("align_clusters_to_truth handles K_hat < K_true", {
  z_true <- c(1, 1, 2, 2, 3, 3)  # 3 clusters
  z_hat <- c(1, 1, 2, 2, 2, 2)   # 2 clusters (cluster 3 merged with 2)
  
  result <- align_clusters_to_truth(z_true, z_hat)
  
  expect_true(is.list(result))
  expect_true(!is.na(result$accuracy))
  expect_true(result$accuracy < 1.0)  # Not perfect
})

test_that("align_clusters_to_truth handles K_hat > K_true", {
  z_true <- c(1, 1, 1, 2, 2, 2)  # 2 clusters
  z_hat <- c(1, 1, 3, 2, 2, 2)   # 3 clusters (cluster 1 split)
  
  result <- align_clusters_to_truth(z_true, z_hat)
  
  expect_true(is.list(result))
  expect_true(!is.na(result$accuracy))
})

test_that("align_clusters_to_truth returns NA for NULL inputs", {
  z_true <- c(1, 1, 2, 2)
  
  result_null1 <- align_clusters_to_truth(NULL, z_true)
  result_null2 <- align_clusters_to_truth(z_true, NULL)
  
  expect_true(is.na(result_null1$accuracy))
  expect_true(is.na(result_null2$accuracy))
})

test_that("align_clusters_to_truth returns NA for empty inputs", {
  z_true <- c(1, 1, 2, 2)
  z_empty <- integer(0)
  
  result <- align_clusters_to_truth(z_empty, z_true)
  expect_true(is.na(result$accuracy))
})

test_that("align_clusters_to_truth returns NA for inputs with NA values", {
  z_true <- c(1, 1, 2, 2)
  z_na <- c(1, NA, 2, 2)
  
  result <- align_clusters_to_truth(z_true, z_na)
  expect_true(is.na(result$accuracy))
})

test_that("align_clusters_to_truth returns NA for length mismatch", {
  z_true <- c(1, 1, 2, 2)
  z_hat <- c(1, 1, 2)
  
  expect_warning(
    result <- align_clusters_to_truth(z_true, z_hat),
    "Length mismatch"
  )
  expect_true(is.na(result$accuracy))
})

# =============================================================================
# Tests for compute_clusterwise_subspace_score
# =============================================================================

test_that("compute_clusterwise_subspace_score works with perfect recovery", {
  # Create true Lambda list
  true_Lambda_list <- list(
    matrix(c(1, 0, 0, 0, 0, 1, 0, 0), nrow = 4, ncol = 2),
    matrix(c(0, 0, 1, 0, 0, 0, 0, 1), nrow = 4, ncol = 2)
  )
  
  # Estimated is identical
  est_Lambda_list <- true_Lambda_list
  
  # Perfect mapping
  mapping <- c(1L, 2L)
  names(mapping) <- c("1", "2")
  
  result <- compute_clusterwise_subspace_score(true_Lambda_list, est_Lambda_list, mapping)
  
  expect_true(is.list(result))
  expect_true(!is.na(result$mean_score))
  expect_equal(result$mean_score, 1.0, tolerance = 1e-6)
  expect_equal(result$n_matched, 2L)
})

test_that("compute_clusterwise_subspace_score handles permuted clusters", {
  # Create true Lambda list
  true_Lambda_list <- list(
    matrix(c(1, 0, 0, 0, 0, 1, 0, 0), nrow = 4, ncol = 2),
    matrix(c(0, 0, 1, 0, 0, 0, 0, 1), nrow = 4, ncol = 2)
  )
  
  # Estimated has permuted order
  est_Lambda_list <- list(
    true_Lambda_list[[2]],  # Cluster 1 in est = Cluster 2 in true
    true_Lambda_list[[1]]   # Cluster 2 in est = Cluster 1 in true
  )
  
  # Mapping reflects the permutation: est cluster 1 -> true cluster 2, etc.
  mapping <- c(2L, 1L)
  names(mapping) <- c("1", "2")
  
  result <- compute_clusterwise_subspace_score(true_Lambda_list, est_Lambda_list, mapping)
  
  expect_true(is.list(result))
  expect_true(!is.na(result$mean_score))
  expect_equal(result$mean_score, 1.0, tolerance = 1e-6)
})

test_that("compute_clusterwise_subspace_score returns NA for NULL inputs", {
  true_Lambda_list <- list(matrix(1:4, 2, 2))
  est_Lambda_list <- list(matrix(1:4, 2, 2))
  mapping <- c(1L)
  names(mapping) <- "1"
  
  result_null1 <- compute_clusterwise_subspace_score(NULL, est_Lambda_list, mapping)
  result_null2 <- compute_clusterwise_subspace_score(true_Lambda_list, NULL, mapping)
  
  expect_true(is.na(result_null1$mean_score))
  expect_true(is.na(result_null2$mean_score))
})

test_that("compute_clusterwise_subspace_score returns NA for NA mapping", {
  true_Lambda_list <- list(matrix(1:4, 2, 2))
  est_Lambda_list <- list(matrix(1:4, 2, 2))
  
  result <- compute_clusterwise_subspace_score(true_Lambda_list, est_Lambda_list, NA)
  expect_true(is.na(result$mean_score))
})

test_that("compute_clusterwise_subspace_score handles missing cluster matches", {
  # 3 true clusters
  true_Lambda_list <- list(
    matrix(c(1, 0, 0, 0), nrow = 2, ncol = 2),
    matrix(c(0, 1, 0, 0), nrow = 2, ncol = 2),
    matrix(c(0, 0, 1, 0), nrow = 2, ncol = 2)
  )
  
  # Only 2 estimated clusters
  est_Lambda_list <- list(
    matrix(c(1, 0, 0, 0), nrow = 2, ncol = 2),
    matrix(c(0, 1, 0, 0), nrow = 2, ncol = 2)
  )
  
  # Mapping: est 1 -> true 1, est 2 -> true 2 (true 3 has no match)
  mapping <- c(1L, 2L)
  names(mapping) <- c("1", "2")
  
  result <- compute_clusterwise_subspace_score(true_Lambda_list, est_Lambda_list, mapping)
  
  expect_true(is.list(result))
  expect_true(!is.na(result$mean_score))
  expect_equal(result$n_matched, 2L)  # Only 2 of 3 clusters matched
  expect_true(is.na(result$scores["3"]))  # Cluster 3 has no match
})

# =============================================================================
# Integration test: Full workflow
# =============================================================================

test_that("subspace metrics work in full workflow", {
  skip_if_not_installed("clue")
  
  # Create synthetic data with known structure
  set.seed(123)
  
  # True cluster assignments
  z_true <- c(rep(1, 5), rep(2, 5))
  
  # Simulated estimated assignments (permuted)
  z_hat <- c(rep(2, 5), rep(1, 5))
  
  # True Lambda list
  true_Lambda_list <- list(
    matrix(rnorm(8), nrow = 4, ncol = 2),
    matrix(rnorm(8), nrow = 4, ncol = 2)
  )
  
  # Estimated Lambda list (same as true but permuted order)
  est_Lambda_list <- list(
    true_Lambda_list[[2]],
    true_Lambda_list[[1]]
  )
  
  # Step 1: Align clusters
  alignment <- align_clusters_to_truth(z_true, z_hat)
  expect_equal(alignment$accuracy, 1.0)
  
  # Step 2: Compute cluster-wise subspace scores
  scores <- compute_clusterwise_subspace_score(
    true_Lambda_list,
    est_Lambda_list,
    alignment$mapping
  )
  
  expect_equal(scores$mean_score, 1.0, tolerance = 1e-6)
  expect_equal(scores$n_matched, 2L)
})
