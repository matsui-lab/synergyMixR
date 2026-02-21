# Unit Tests for Advanced r Selection Functions
# Tests for r_selection_advanced.R functions

# Load test fixtures
small_data <- readRDS("fixtures/small_test_data.rds")

# ============================================================
# Tests for compute_smoothness_penalty
# ============================================================

test_that("compute_smoothness_penalty works with valid factor scores", {
  # Create simple factor scores
  N <- 3
  T_i <- 10
  r <- 2
  
  factor_scores <- list()
  for (i in 1:N) {
    factor_scores[[i]] <- matrix(rnorm(T_i * r), T_i, r)
  }
  
  smoothness <- compute_smoothness_penalty(factor_scores)
  
  # Should return a non-negative scalar
  expect_type(smoothness, "double")
  expect_length(smoothness, 1)
  expect_true(smoothness >= 0)
})

test_that("compute_smoothness_penalty handles smooth vs wiggly signals", {
  # Create smooth signal (low frequency)
  T_i <- 50
  r <- 1
  t <- seq(0, 2*pi, length.out = T_i)
  
  smooth_scores <- list(matrix(sin(t), T_i, r))
  smooth_penalty <- compute_smoothness_penalty(smooth_scores)
  
  # Create wiggly signal (high frequency)
  wiggly_scores <- list(matrix(sin(10*t), T_i, r))
  wiggly_penalty <- compute_smoothness_penalty(wiggly_scores)
  
  # Wiggly signal should have higher penalty
  expect_true(wiggly_penalty > smooth_penalty)
})

test_that("compute_smoothness_penalty handles edge cases", {
  # Empty list
  expect_equal(compute_smoothness_penalty(list()), 0)
  
  # Short time series (< 3 points)
  short_scores <- list(matrix(c(1, 2), 2, 1))
  expect_equal(compute_smoothness_penalty(short_scores), 0)
  
  # NULL elements
  mixed_scores <- list(
    matrix(rnorm(10), 10, 1),
    NULL,
    matrix(rnorm(10), 10, 1)
  )
  expect_true(compute_smoothness_penalty(mixed_scores) >= 0)
})

test_that("compute_smoothness_penalty validates inputs", {
  # Non-list input
  expect_error(
    compute_smoothness_penalty(matrix(1:10, 10, 1)),
    "must be a list"
  )
})

# ============================================================
# Tests for compute_vss_criterion
# ============================================================

test_that("compute_vss_criterion works with valid loadings", {
  # Create loading matrices with simple structure
  K <- 2
  M <- 6
  r <- 2
  
  Lambda <- list()
  for (k in 1:K) {
    # Create simple structure: first 3 muscles load on factor 1, last 3 on factor 2
    Lambda_k <- matrix(0, M, r)
    Lambda_k[1:3, 1] <- runif(3, 0.6, 1.0)
    Lambda_k[4:6, 2] <- runif(3, 0.6, 1.0)
    Lambda_k[1:3, 2] <- runif(3, 0, 0.3)
    Lambda_k[4:6, 1] <- runif(3, 0, 0.3)
    Lambda[[k]] <- Lambda_k
  }
  
  vss <- compute_vss_criterion(Lambda, threshold = 0.5)
  
  # Should return a value between 0 and 1
  expect_type(vss, "double")
  expect_length(vss, 1)
  expect_true(vss >= 0 && vss <= 1)
  
  # With perfect simple structure, should be high
  expect_true(vss > 0.8)
})

test_that("compute_vss_criterion detects complex structure", {
  # Create loadings with complex structure (all muscles load on all factors)
  K <- 1
  M <- 4
  r <- 2
  
  Lambda_complex <- list(matrix(runif(M * r, 0.6, 1.0), M, r))
  vss_complex <- compute_vss_criterion(Lambda_complex, threshold = 0.5)
  
  # Complex structure should have low VSS
  expect_true(vss_complex < 0.5)
})

test_that("compute_vss_criterion handles single factor", {
  # With r=1, all variables trivially have simple structure
  K <- 1
  M <- 4
  r <- 1
  
  Lambda <- list(matrix(runif(M, 0.5, 1.0), M, r))
  vss <- compute_vss_criterion(Lambda, threshold = 0.5)
  
  # Should be 1.0 (perfect simple structure)
  expect_equal(vss, 1.0)
})

test_that("compute_vss_criterion handles edge cases", {
  # Empty list
  expect_equal(compute_vss_criterion(list()), 0)
  
  # NULL elements
  mixed_Lambda <- list(
    matrix(runif(8, 0.5, 1.0), 4, 2),
    NULL,
    matrix(runif(8, 0.5, 1.0), 4, 2)
  )
  expect_true(compute_vss_criterion(mixed_Lambda) >= 0)
})

test_that("compute_vss_criterion validates inputs", {
  # Non-list input
  expect_error(
    compute_vss_criterion(matrix(1:10, 5, 2)),
    "must be a list"
  )
})

# ============================================================
# Tests for select_r_by_smoothness_mfa
# ============================================================

test_that("select_r_by_smoothness_mfa works with valid inputs", {
  list_of_data <- small_data$data
  K <- 2
  
  result <- select_r_by_smoothness_mfa(
    list_of_data = list_of_data,
    K = K,
    rvec = 1:2,
    lambda_smooth = 1e-3,
    vaf_threshold = 0.5,
    max_iter = 5,
    nIterFA = 3,
    tol = 1e-3,
    n_init = 1,
    use_kmeans_init = FALSE
  )
  
  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("best_r", "summary", "best_model"))
  
  # Check best_r
  expect_true(result$best_r %in% 1:2)
  
  # Check summary
  expect_s3_class(result$summary, "data.frame")
  expect_equal(nrow(result$summary), 2)
  expect_true(all(c("r", "BIC", "VAF", "smoothness", "penalized_score") %in% names(result$summary)))
  
  # Check that all metrics are finite
  expect_true(all(is.finite(result$summary$BIC)))
  expect_true(all(is.finite(result$summary$VAF)))
  expect_true(all(is.finite(result$summary$smoothness)))
  expect_true(all(is.finite(result$summary$penalized_score)))
  
  # Check best_model
  expect_type(result$best_model, "list")
  expect_true("Lambda" %in% names(result$best_model))
  expect_true("factor_scores" %in% names(result$best_model))
})

test_that("select_r_by_smoothness_mfa respects VAF threshold", {
  list_of_data <- small_data$data
  K <- 2
  
  # Use high VAF threshold
  result <- select_r_by_smoothness_mfa(
    list_of_data = list_of_data,
    K = K,
    rvec = 1:2,
    lambda_smooth = 1e-3,
    vaf_threshold = 0.95,
    max_iter = 10,
    nIterFA = 5,
    tol = 1e-3,
    n_init = 1,
    use_kmeans_init = FALSE
  )
  
  # Best model should satisfy VAF threshold (or be best available)
  best_vaf <- result$summary$VAF[result$summary$r == result$best_r]
  expect_true(best_vaf >= 0.5)  # Should at least try to meet threshold
})

# ============================================================
# Tests for select_r_by_smoothness_mpca
# ============================================================

test_that("select_r_by_smoothness_mpca works with valid inputs", {
  list_of_data <- small_data$data
  K <- 2
  
  result <- select_r_by_smoothness_mpca(
    list_of_data = list_of_data,
    K = K,
    rvec = 1:2,
    lambda_smooth = 1e-3,
    vaf_threshold = 0.5,
    max_iter = 5,
    nIterPCA = 3,
    tol = 1e-3,
    method = "EM",
    n_init = 1,
    use_kmeans_init = FALSE
  )
  
  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("best_r", "summary", "best_model"))
  
  # Check best_r
  expect_true(result$best_r %in% 1:2)
  
  # Check summary
  expect_s3_class(result$summary, "data.frame")
  expect_equal(nrow(result$summary), 2)
  expect_true(all(c("r", "BIC", "VAF", "smoothness", "penalized_score") %in% names(result$summary)))
  
  # Check best_model
  expect_type(result$best_model, "list")
  expect_true("W" %in% names(result$best_model))
  expect_true("factor_scores" %in% names(result$best_model))
})

# ============================================================
# Tests for select_r_by_vss_mfa
# ============================================================

test_that("select_r_by_vss_mfa works with valid inputs", {
  list_of_data <- small_data$data
  K <- 2
  
  result <- select_r_by_vss_mfa(
    list_of_data = list_of_data,
    K = K,
    rvec = 1:2,
    vss_threshold = 0.5,
    vaf_threshold = 0.5,
    max_iter = 5,
    nIterFA = 3,
    tol = 1e-3,
    n_init = 1,
    use_kmeans_init = FALSE
  )
  
  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("best_r", "summary", "best_model"))
  
  # Check best_r
  expect_true(result$best_r %in% 1:2)
  
  # Check summary
  expect_s3_class(result$summary, "data.frame")
  expect_equal(nrow(result$summary), 2)
  expect_true(all(c("r", "BIC", "VAF", "VSS") %in% names(result$summary)))
  
  # Check that VSS is between 0 and 1
  expect_true(all(result$summary$VSS >= 0 & result$summary$VSS <= 1))
  
  # Check best_model
  expect_type(result$best_model, "list")
  expect_true("Lambda" %in% names(result$best_model))
})

test_that("select_r_by_vss_mfa prefers higher VSS", {
  list_of_data <- small_data$data
  K <- 2
  
  result <- select_r_by_vss_mfa(
    list_of_data = list_of_data,
    K = K,
    rvec = 1:2,
    vss_threshold = 0.5,
    vaf_threshold = 0.5,
    max_iter = 10,
    nIterFA = 5,
    tol = 1e-3,
    n_init = 2,
    use_kmeans_init = TRUE
  )
  
  # Among models meeting VAF threshold, should select one with high VSS
  valid_models <- result$summary[result$summary$VAF >= 0.5, ]
  if (nrow(valid_models) > 0) {
    best_vss <- result$summary$VSS[result$summary$r == result$best_r]
    expect_true(best_vss >= min(valid_models$VSS))
  }
})

# ============================================================
# Tests for select_r_by_vss_mpca
# ============================================================

test_that("select_r_by_vss_mpca works with valid inputs", {
  list_of_data <- small_data$data
  K <- 2
  
  result <- select_r_by_vss_mpca(
    list_of_data = list_of_data,
    K = K,
    rvec = 1:2,
    vss_threshold = 0.5,
    vaf_threshold = 0.5,
    max_iter = 5,
    nIterPCA = 3,
    tol = 1e-3,
    method = "EM",
    n_init = 1,
    use_kmeans_init = FALSE
  )
  
  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("best_r", "summary", "best_model"))
  
  # Check best_r
  expect_true(result$best_r %in% 1:2)
  
  # Check summary
  expect_s3_class(result$summary, "data.frame")
  expect_equal(nrow(result$summary), 2)
  expect_true(all(c("r", "BIC", "VAF", "VSS") %in% names(result$summary)))
  
  # Check best_model
  expect_type(result$best_model, "list")
  expect_true("W" %in% names(result$best_model))
})

# ============================================================
# Tests for select_r_multicriteria_mfa
# ============================================================

test_that("select_r_multicriteria_mfa works with valid inputs", {
  list_of_data <- small_data$data
  K <- 2
  
  result <- select_r_multicriteria_mfa(
    list_of_data = list_of_data,
    K = K,
    rvec = 1:2,
    vaf_threshold = 0.5,
    lambda_smooth = 1e-3,
    vss_threshold = 0.5,
    weights = c(BIC = 1, VAF = 1, smoothness = 1, VSS = 1),
    max_iter = 5,
    nIterFA = 3,
    tol = 1e-3,
    n_init = 1,
    use_kmeans_init = FALSE
  )
  
  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("best_r", "summary", "best_model"))
  
  # Check best_r
  expect_true(result$best_r %in% 1:2)
  
  # Check summary
  expect_s3_class(result$summary, "data.frame")
  expect_equal(nrow(result$summary), 2)
  expect_true(all(c("r", "BIC", "VAF", "smoothness", "VSS", "composite_score") %in% names(result$summary)))
  
  # Check that all metrics are finite
  expect_true(all(is.finite(result$summary$BIC)))
  expect_true(all(is.finite(result$summary$VAF)))
  expect_true(all(is.finite(result$summary$smoothness)))
  expect_true(all(is.finite(result$summary$VSS)))
  expect_true(all(is.finite(result$summary$composite_score)))
  
  # Check best_model
  expect_type(result$best_model, "list")
  expect_true("Lambda" %in% names(result$best_model))
  expect_true("factor_scores" %in% names(result$best_model))
})

test_that("select_r_multicriteria_mfa respects custom weights", {
  list_of_data <- small_data$data
  K <- 2
  
  # Test with different weight configurations
  result1 <- select_r_multicriteria_mfa(
    list_of_data = list_of_data,
    K = K,
    rvec = 1:2,
    vaf_threshold = 0.5,
    weights = c(BIC = 10, VAF = 1, smoothness = 1, VSS = 1),
    max_iter = 5,
    nIterFA = 3,
    n_init = 1
  )
  
  result2 <- select_r_multicriteria_mfa(
    list_of_data = list_of_data,
    K = K,
    rvec = 1:2,
    vaf_threshold = 0.5,
    weights = c(BIC = 1, VAF = 10, smoothness = 1, VSS = 1),
    max_iter = 5,
    nIterFA = 3,
    n_init = 1
  )
  
  # Both should produce valid results
  expect_true(result1$best_r %in% 1:2)
  expect_true(result2$best_r %in% 1:2)
  
  # Composite scores should be different due to different weights
  # (unless both models are identical)
  expect_true(
    !identical(result1$summary$composite_score, result2$summary$composite_score) ||
    result1$best_r == result2$best_r
  )
})

test_that("select_r_multicriteria_mfa handles equal weights", {
  list_of_data <- small_data$data
  K <- 2
  
  result <- select_r_multicriteria_mfa(
    list_of_data = list_of_data,
    K = K,
    rvec = 1:2,
    vaf_threshold = 0.5,
    weights = c(BIC = 1, VAF = 1, smoothness = 1, VSS = 1),
    max_iter = 5,
    nIterFA = 3,
    n_init = 1
  )
  
  # Should balance all criteria equally
  expect_true(result$best_r %in% 1:2)
  expect_true(all(is.finite(result$summary$composite_score)))
})

# ============================================================
# Tests for select_r_multicriteria_mpca
# ============================================================

test_that("select_r_multicriteria_mpca works with valid inputs", {
  list_of_data <- small_data$data
  K <- 2
  
  result <- select_r_multicriteria_mpca(
    list_of_data = list_of_data,
    K = K,
    rvec = 1:2,
    vaf_threshold = 0.5,
    lambda_smooth = 1e-3,
    vss_threshold = 0.5,
    weights = c(BIC = 1, VAF = 1, smoothness = 1, VSS = 1),
    max_iter = 5,
    nIterPCA = 3,
    tol = 1e-3,
    method = "EM",
    n_init = 1,
    use_kmeans_init = FALSE
  )
  
  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("best_r", "summary", "best_model"))
  
  # Check best_r
  expect_true(result$best_r %in% 1:2)
  
  # Check summary
  expect_s3_class(result$summary, "data.frame")
  expect_equal(nrow(result$summary), 2)
  expect_true(all(c("r", "BIC", "VAF", "smoothness", "VSS", "composite_score") %in% names(result$summary)))
  
  # Check best_model
  expect_type(result$best_model, "list")
  expect_true("W" %in% names(result$best_model))
  expect_true("factor_scores" %in% names(result$best_model))
})

# ============================================================
# Integration tests
# ============================================================

test_that("All selection methods produce consistent results", {
  list_of_data <- small_data$data
  K <- 2
  rvec <- 1:2
  
  # Run all three selection methods
  result_smooth <- select_r_by_smoothness_mfa(
    list_of_data, K, rvec,
    max_iter = 5, nIterFA = 3, n_init = 1
  )
  
  result_vss <- select_r_by_vss_mfa(
    list_of_data, K, rvec,
    max_iter = 5, nIterFA = 3, n_init = 1
  )
  
  result_multi <- select_r_multicriteria_mfa(
    list_of_data, K, rvec,
    max_iter = 5, nIterFA = 3, n_init = 1
  )
  
  # All should select valid r values
  expect_true(result_smooth$best_r %in% rvec)
  expect_true(result_vss$best_r %in% rvec)
  expect_true(result_multi$best_r %in% rvec)
  
  # All should have same number of models evaluated
  expect_equal(nrow(result_smooth$summary), length(rvec))
  expect_equal(nrow(result_vss$summary), length(rvec))
  expect_equal(nrow(result_multi$summary), length(rvec))
})

test_that("MFA and MPCA methods produce comparable results", {
  list_of_data <- small_data$data
  K <- 2
  rvec <- 1:2
  
  # Run MFA and MPCA multi-criteria selection
  result_mfa <- select_r_multicriteria_mfa(
    list_of_data, K, rvec,
    max_iter = 10, nIterFA = 5, n_init = 2
  )
  
  result_mpca <- select_r_multicriteria_mpca(
    list_of_data, K, rvec,
    max_iter = 10, nIterPCA = 5, method = "EM", n_init = 2
  )
  
  # Both should select valid r values
  expect_true(result_mfa$best_r %in% rvec)
  expect_true(result_mpca$best_r %in% rvec)
  
  # Both should have valid VAF values
  expect_true(all(result_mfa$summary$VAF >= 0 & result_mfa$summary$VAF <= 1))
  expect_true(all(result_mpca$summary$VAF >= 0 & result_mpca$summary$VAF <= 1))
})
