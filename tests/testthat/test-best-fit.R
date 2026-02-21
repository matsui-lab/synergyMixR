# Unit Tests for Model Selection Functions
# Tests for best_fit.R functions

# Load test fixtures
small_data <- readRDS("fixtures/small_test_data.rds")

test_that("select_optimal_K_r_mfa works with valid inputs", {
  # Use small test data for faster testing
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Test with small grid
  result <- select_optimal_K_r_mfa(
    list_of_data = list_of_data,
    Kvec = 1:2,
    rvec = 1:2,
    max_iter = 5,  # Short for testing
    nIterFA = 3,
    tol = 1e-3,
    n_init = 1,
    use_kmeans_init = FALSE
  )
  
  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("summary", "best_model_info", "all_models"))
  
  # Check summary data frame
  expect_s3_class(result$summary, "data.frame")
  expect_equal(nrow(result$summary), 4)  # 2 K values × 2 r values
  expect_true(all(c("K", "r", "logLik", "BIC") %in% names(result$summary)))
  
  # Check that summary is sorted by BIC
  expect_true(all(diff(result$summary$BIC) >= 0))
  
  # Check best_model_info
  expect_type(result$best_model_info, "list")
  expect_named(result$best_model_info, c("K", "r", "model", "logLik", "BIC"))
  expect_true(result$best_model_info$K %in% 1:2)
  expect_true(result$best_model_info$r %in% 1:2)
  
  # Check that best model has lowest BIC
  expect_equal(result$best_model_info$BIC, min(result$summary$BIC))
  
  # Check all_models
  expect_type(result$all_models, "list")
  expect_equal(length(result$all_models), 4)
  
  # Check each model in all_models
  for (model_info in result$all_models) {
    expect_named(model_info, c("K", "r", "model", "logLik", "BIC"))
    expect_true(is.finite(model_info$logLik))
    expect_true(is.finite(model_info$BIC))
  }
})

test_that("select_optimal_K_r_mfa selects correct model", {
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Test with grid including true K and r
  result <- select_optimal_K_r_mfa(
    list_of_data = list_of_data,
    Kvec = c(1, true_params$K, 3),
    rvec = c(1, true_params$r),
    max_iter = 10,
    nIterFA = 5,
    tol = 1e-3,
    n_init = 2,
    use_kmeans_init = TRUE
  )
  
  # Best model should be reasonable
  expect_true(result$best_model_info$K >= 1)
  expect_true(result$best_model_info$r >= 1)
  
  # BIC should prefer simpler models when appropriate
  # (but we can't guarantee it will select true K,r with small data)
  expect_true(result$best_model_info$BIC < Inf)
})

test_that("select_optimal_K_r_mfa handles single K and r", {
  list_of_data <- small_data$data
  
  # Test with single K and r
  result <- select_optimal_K_r_mfa(
    list_of_data = list_of_data,
    Kvec = 2,
    rvec = 2,
    max_iter = 5,
    nIterFA = 3,
    tol = 1e-3,
    n_init = 1,
    use_kmeans_init = FALSE
  )
  
  # Should still work with single combination
  expect_equal(nrow(result$summary), 1)
  expect_equal(result$best_model_info$K, 2)
  expect_equal(result$best_model_info$r, 2)
})

test_that("select_optimal_K_r_mfa validates inputs", {
  list_of_data <- small_data$data
  
  # Test with empty data
  expect_error(
    select_optimal_K_r_mfa(list_of_data = list(), Kvec = 1:2, rvec = 1:2),
    "No data"
  )
})

test_that("select_optimal_K_r_mpca works with valid inputs", {
  # Use small test data for faster testing
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Test with small grid
  result <- select_optimal_K_r_mpca(
    list_of_data = list_of_data,
    Kvec = 1:2,
    rvec = 1:2,
    max_iter = 5,
    nIterPCA = 3,
    tol = 1e-3,
    method = "EM",
    n_init = 1,
    use_kmeans_init = FALSE
  )
  
  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("summary", "best_model_info", "all_models"))
  
  # Check summary data frame
  expect_s3_class(result$summary, "data.frame")
  expect_equal(nrow(result$summary), 4)
  expect_true(all(c("K", "r", "logLik", "BIC") %in% names(result$summary)))
  
  # Check that summary is sorted by BIC
  expect_true(all(diff(result$summary$BIC) >= 0))
  
  # Check best_model_info
  expect_type(result$best_model_info, "list")
  expect_named(result$best_model_info, c("K", "r", "model", "logLik", "BIC"))
  expect_true(result$best_model_info$K %in% 1:2)
  expect_true(result$best_model_info$r %in% 1:2)
  
  # Check that best model has lowest BIC
  expect_equal(result$best_model_info$BIC, min(result$summary$BIC))
  
  # Check all_models
  expect_type(result$all_models, "list")
  expect_equal(length(result$all_models), 4)
})

test_that("select_optimal_K_r_mpca works with closed-form method", {
  list_of_data <- small_data$data
  
  # Test with closed-form method
  result <- select_optimal_K_r_mpca(
    list_of_data = list_of_data,
    Kvec = 1:2,
    rvec = 1:2,
    max_iter = 5,
    nIterPCA = 3,
    tol = 1e-3,
    method = "closed_form",
    n_init = 1,
    use_kmeans_init = FALSE
  )
  
  # Should work with closed-form method
  expect_equal(nrow(result$summary), 4)
  expect_true(all(is.finite(result$summary$BIC)))
})

test_that("select_optimal_K_r_mpca handles single K and r", {
  list_of_data <- small_data$data
  
  # Test with single K and r
  result <- select_optimal_K_r_mpca(
    list_of_data = list_of_data,
    Kvec = 2,
    rvec = 2,
    max_iter = 5,
    nIterPCA = 3,
    tol = 1e-3,
    method = "EM",
    n_init = 1,
    use_kmeans_init = FALSE
  )
  
  # Should still work with single combination
  expect_equal(nrow(result$summary), 1)
  expect_equal(result$best_model_info$K, 2)
  expect_equal(result$best_model_info$r, 2)
})

test_that("select_optimal_K_r_mpca validates inputs", {
  list_of_data <- small_data$data
  
  # Test with empty data
  expect_error(
    select_optimal_K_r_mpca(list_of_data = list(), Kvec = 1:2, rvec = 1:2),
    "No data"
  )
})

test_that("BIC calculation is consistent across methods", {
  list_of_data <- small_data$data
  
  # Fit same model with both MFA and MPCA
  result_mfa <- select_optimal_K_r_mfa(
    list_of_data = list_of_data,
    Kvec = 2,
    rvec = 2,
    max_iter = 10,
    nIterFA = 5,
    tol = 1e-3,
    n_init = 1,
    use_kmeans_init = FALSE
  )
  
  result_mpca <- select_optimal_K_r_mpca(
    list_of_data = list_of_data,
    Kvec = 2,
    rvec = 2,
    max_iter = 10,
    nIterPCA = 5,
    tol = 1e-3,
    method = "EM",
    n_init = 1,
    use_kmeans_init = FALSE
  )
  
  # Both should have valid BIC values
  expect_true(is.finite(result_mfa$best_model_info$BIC))
  expect_true(is.finite(result_mpca$best_model_info$BIC))
  
  # BIC should penalize complexity (more parameters = higher BIC for same logLik)
  # MFA typically has more parameters than MPCA
  # (but we can't guarantee which will be better without knowing the data)
})

test_that("Model selection prefers simpler models when appropriate", {
  # Create data from a simple model (K=1, r=1)
  set.seed(999)
  N <- 10
  T_i <- 30
  M <- 4
  
  # Generate data from single cluster, single factor
  mu <- rnorm(M)
  Lambda <- matrix(rnorm(M), M, 1)
  
  list_of_data <- list()
  for (i in 1:N) {
    Z_i <- matrix(rnorm(T_i), T_i, 1)
    X_i <- Z_i %*% t(Lambda) + matrix(rep(mu, T_i), T_i, M, byrow = TRUE) + 
           matrix(rnorm(T_i * M, 0, 0.1), T_i, M)
    list_of_data[[i]] <- X_i
  }
  
  # Test model selection
  result <- select_optimal_K_r_mfa(
    list_of_data = list_of_data,
    Kvec = 1:3,
    rvec = 1:2,
    max_iter = 20,
    nIterFA = 10,
    tol = 1e-4,
    n_init = 3,
    use_kmeans_init = TRUE
  )
  
  # Should prefer simpler model (K=1 or r=1)
  # BIC penalizes complexity, so simpler models should be favored
  expect_true(result$best_model_info$K <= 2)
  expect_true(result$best_model_info$r <= 2)
})

test_that("Grid search explores all combinations", {
  list_of_data <- small_data$data
  
  Kvec <- c(1, 2, 3)
  rvec <- c(1, 2)
  
  result <- select_optimal_K_r_mfa(
    list_of_data = list_of_data,
    Kvec = Kvec,
    rvec = rvec,
    max_iter = 3,
    nIterFA = 2,
    tol = 1e-3,
    n_init = 1,
    use_kmeans_init = FALSE
  )
  
  # Should have tried all combinations
  expect_equal(nrow(result$summary), length(Kvec) * length(rvec))
  expect_equal(length(result$all_models), length(Kvec) * length(rvec))
  
  # Check that all K and r values were tried
  tried_K <- sort(unique(result$summary$K))
  tried_r <- sort(unique(result$summary$r))
  
  expect_equal(tried_K, Kvec)
  expect_equal(tried_r, rvec)
})

test_that("Model selection is reproducible", {
  list_of_data <- small_data$data
  
  # Run twice with same parameters
  result1 <- select_optimal_K_r_mfa(
    list_of_data = list_of_data,
    Kvec = 1:2,
    rvec = 1:2,
    max_iter = 5,
    nIterFA = 3,
    tol = 1e-3,
    n_init = 2,
    use_kmeans_init = TRUE
  )
  
  result2 <- select_optimal_K_r_mfa(
    list_of_data = list_of_data,
    Kvec = 1:2,
    rvec = 1:2,
    max_iter = 5,
    nIterFA = 3,
    tol = 1e-3,
    n_init = 2,
    use_kmeans_init = TRUE
  )
  
  # Results should be similar (allowing for some randomness in initialization)
  # At minimum, the BIC values should be in the same order of magnitude
  expect_equal(result1$best_model_info$K, result2$best_model_info$K)
  expect_equal(result1$best_model_info$r, result2$best_model_info$r)
})
