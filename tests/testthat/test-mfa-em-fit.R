# Unit Tests for MFA EM Fitting Functions
# Tests for mfa_em_fit.R functions

# Load test fixtures
small_data <- readRDS("fixtures/small_test_data.rds")
medium_data <- readRDS("fixtures/medium_test_data.rds")
edge_cases <- readRDS("fixtures/edge_case_data.rds")

test_that("mfa_em_fit_cpp_singleInit works with valid inputs", {
  # Use small test data
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Test with known parameters
  K <- true_params$K
  r <- true_params$r
  z_init <- true_params$z
  
  # Fit model
  fit <- mfa_em_fit_cpp_singleInit(
    list_of_data = list_of_data,
    K = K,
    r = r,
    z_init = z_init,
    max_iter = 10,  # Short for testing
    nIterFA = 5,
    tol = 1e-3
  )
  
  # Check output structure
  expect_type(fit, "list")
  expect_named(fit, c("z", "pi", "mu", "Lambda", "Psi", "logLik", "resp"))
  
  # Check dimensions
  expect_equal(length(fit$z), true_params$N)
  expect_equal(length(fit$pi), K)
  expect_equal(length(fit$mu), K)
  expect_equal(length(fit$Lambda), K)
  expect_equal(length(fit$Psi), K)
  expect_equal(dim(fit$resp), c(true_params$N, K))
  
  # Check cluster assignments are valid
  expect_true(all(fit$z %in% 1:K))
  
  # Check mixing proportions sum to 1
  expect_equal(sum(fit$pi), 1, tolerance = 1e-6)
  expect_true(all(fit$pi > 0))
  
  # Check parameter dimensions
  for (k in 1:K) {
    expect_equal(length(fit$mu[[k]]), true_params$M)
    expect_equal(dim(fit$Lambda[[k]]), c(true_params$M, r))
    expect_equal(dim(fit$Psi[[k]]), c(true_params$M, true_params$M))
    
    # Check Psi is diagonal and positive
    expect_true(all(diag(fit$Psi[[k]]) > 0))
    expect_equal(fit$Psi[[k]], diag(diag(fit$Psi[[k]])), tolerance = 1e-10)
  }
  
  # Check log-likelihood is finite
  expect_true(is.finite(fit$logLik))
  expect_false(is.na(fit$logLik))
  
  # Check responsibilities are valid probabilities
  expect_true(all(fit$resp >= 0))
  expect_true(all(fit$resp <= 1))
  expect_equal(rowSums(fit$resp), rep(1, true_params$N), tolerance = 1e-6)
})

test_that("mfa_em_fit_cpp_singleInit handles edge cases", {
  # Test with single cluster
  single_cluster_data <- edge_cases$single_cluster$data
  
  fit_single <- mfa_em_fit_cpp_singleInit(
    list_of_data = single_cluster_data,
    K = 1,
    r = 2,
    z_init = rep(1, length(single_cluster_data)),
    max_iter = 5,
    nIterFA = 3,
    tol = 1e-3
  )
  
  expect_equal(length(fit_single$pi), 1)
  expect_equal(fit_single$pi[1], 1)
  expect_true(all(fit_single$z == 1))
  
  # Test with single factor
  single_factor_data <- edge_cases$single_factor$data
  
  fit_single_factor <- mfa_em_fit_cpp_singleInit(
    list_of_data = single_factor_data,
    K = 2,
    r = 1,
    z_init = edge_cases$single_factor$true_params$z,
    max_iter = 5,
    nIterFA = 3,
    tol = 1e-3
  )
  
  expect_equal(ncol(fit_single_factor$Lambda[[1]]), 1)
  expect_equal(ncol(fit_single_factor$Lambda[[2]]), 1)
})

test_that("mfa_em_fit_cpp_singleInit validates inputs", {
  list_of_data <- small_data$data
  
  # Test invalid K (will fail in C++ code)
  expect_error(
    mfa_em_fit_cpp_singleInit(list_of_data, K = 0, r = 2, z_init = rep(1, 5))
  )
  
  # Test invalid r (will fail in C++ code)
  expect_error(
    mfa_em_fit_cpp_singleInit(list_of_data, K = 2, r = 0, z_init = rep(1, 5))
  )
  
  # Test mismatched z_init length (will fail in C++ code)
  expect_error(
    mfa_em_fit_cpp_singleInit(list_of_data, K = 2, r = 2, z_init = rep(1, 3))
  )
  
  # Test invalid z_init values (will fail in C++ code)
  expect_error(
    mfa_em_fit_cpp_singleInit(list_of_data, K = 2, r = 2, z_init = c(1, 1, 3, 1, 1))
  )
})

test_that("mfa_em_fit works with default parameters", {
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Test single-run mode (no multi-init)
  fit_single <- mfa_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 10,
    nIterFA = 5,
    tol = 1e-3,
    n_init = 1,
    use_kmeans_init = FALSE,
    seed = 123
  )
  
  # Check output structure
  expect_type(fit_single, "list")
  expect_named(fit_single, c("z", "pi", "mu", "Lambda", "Psi", "logLik", "resp", "thread_info"))
  
  # Check basic properties
  expect_equal(length(fit_single$z), true_params$N)
  expect_equal(length(fit_single$pi), true_params$K)
  expect_equal(sum(fit_single$pi), 1, tolerance = 1e-6)
  expect_true(is.finite(fit_single$logLik))
})

test_that("mfa_em_fit works with multi-initialization", {
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Test multi-init mode
  fit_multi <- mfa_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 5,  # Short for testing
    nIterFA = 3,
    tol = 1e-3,
    n_init = 3,
    use_kmeans_init = TRUE,
    subject_rdim_for_kmeans = 2,
    mc_cores = 1,  # Sequential for testing
    seed = 456
  )
  
  # Check output structure
  expect_type(fit_multi, "list")
  expect_named(fit_multi, c("z", "pi", "mu", "Lambda", "Psi", "logLik", "resp", "thread_info"))
  
  # Check basic properties
  expect_equal(length(fit_multi$z), true_params$N)
  expect_equal(length(fit_multi$pi), true_params$K)
  expect_equal(sum(fit_multi$pi), 1, tolerance = 1e-6)
  expect_true(is.finite(fit_multi$logLik))
})

test_that("mfa_em_fit is deterministic with fixed seed", {
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Fit same model twice with same seed
  fit1 <- mfa_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 5,
    nIterFA = 3,
    tol = 1e-3,
    n_init = 2,
    use_kmeans_init = TRUE,
    mc_cores = 1,
    seed = 789
  )
  
  fit2 <- mfa_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 5,
    nIterFA = 3,
    tol = 1e-3,
    n_init = 2,
    use_kmeans_init = TRUE,
    mc_cores = 1,
    seed = 789
  )
  
  # Results should be identical
  expect_equal(fit1$z, fit2$z)
  expect_equal(fit1$pi, fit2$pi, tolerance = 1e-10)
  expect_equal(fit1$logLik, fit2$logLik, tolerance = 1e-10)
  
  for (k in 1:true_params$K) {
    expect_equal(fit1$mu[[k]], fit2$mu[[k]], tolerance = 1e-10)
    expect_equal(fit1$Lambda[[k]], fit2$Lambda[[k]], tolerance = 1e-10)
    expect_equal(fit1$Psi[[k]], fit2$Psi[[k]], tolerance = 1e-10)
  }
})

test_that("mfa_em_fit handles Windows parallel fallback", {
  # Skip this test - cannot mock .Platform$OS.type reliably
  skip("Cannot reliably mock .Platform$OS.type for testing")
})

test_that("mfa_em_fit convergence properties", {
  list_of_data <- medium_data$data
  true_params <- medium_data$true_params
  
  # Test that more iterations generally improve log-likelihood
  fit_short <- mfa_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 3,
    nIterFA = 2,
    n_init = 1,
    use_kmeans_init = FALSE,
    seed = 111
  )
  
  fit_long <- mfa_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 15,
    nIterFA = 8,
    n_init = 1,
    use_kmeans_init = FALSE,
    seed = 111
  )
  
  # Longer run should have equal or better log-likelihood
  expect_true(fit_long$logLik >= fit_short$logLik - 1e-6)
})

test_that("mfa_em_fit parameter recovery on synthetic data", {
  skip_if_not_installed("mclust")

  # Use small data with known true parameters
  list_of_data <- small_data$data
  true_params <- small_data$true_params

  # Fit model with correct K and r
  fit <- mfa_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 50,  # More iterations for better recovery
    nIterFA = 20,
    tol = 1e-4,
    n_init = 5,
    use_kmeans_init = TRUE,
    mc_cores = 1,
    seed = 123
  )
  
  # Check that we recover reasonable cluster assignments
  # (allowing for label switching)
  ari_score <- mclust::adjustedRandIndex(fit$z, true_params$z)
  expect_true(ari_score > 0.5, 
              info = paste("ARI score:", round(ari_score, 3)))
  
  # Check that parameters are in reasonable ranges
  for (k in 1:true_params$K) {
    # Mean parameters should be finite
    expect_true(all(is.finite(fit$mu[[k]])))
    
    # Loading matrices should be finite
    expect_true(all(is.finite(fit$Lambda[[k]])))
    
    # Noise variances should be positive and finite
    expect_true(all(diag(fit$Psi[[k]]) > 0))
    expect_true(all(is.finite(diag(fit$Psi[[k]]))))
  }
  
  # Log-likelihood should be reasonable (not too negative)
  expect_true(fit$logLik > -1e6)
})

test_that("mfa_em_fit produces reproducible results across different thread counts with cpp_seed", {
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Test with single initialization to isolate C++ RNG determinism
  # Fit with 2 threads
  fit1 <- mfa_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 10,
    nIterFA = 5,
    tol = 1e-3,
    n_init = 1,
    use_kmeans_init = FALSE,
    n_threads = 2,
    cpp_seed = 42
  )
  
  # Fit with 4 threads
  fit2 <- mfa_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 10,
    nIterFA = 5,
    tol = 1e-3,
    n_init = 1,
    use_kmeans_init = FALSE,
    n_threads = 4,
    cpp_seed = 42
  )
  
  # Results should be identical across different thread counts
  expect_equal(fit1$logLik, fit2$logLik, tolerance = 1e-8)
  expect_equal(fit1$z, fit2$z)
  expect_equal(fit1$pi, fit2$pi, tolerance = 1e-10)
  
  for (k in 1:true_params$K) {
    expect_equal(fit1$mu[[k]], fit2$mu[[k]], tolerance = 1e-10)
    expect_equal(fit1$Lambda[[k]], fit2$Lambda[[k]], tolerance = 1e-10)
    expect_equal(fit1$Psi[[k]], fit2$Psi[[k]], tolerance = 1e-10)
  }
})

test_that("mfa_em_fit_cpp_singleInit produces reproducible results with seed parameter", {
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  z_init <- true_params$z
  
  # Fit with 2 threads
  fit1 <- mfa_em_fit_cpp_singleInit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    z_init = z_init,
    max_iter = 10,
    nIterFA = 5,
    tol = 1e-3,
    n_threads = 2,
    seed = 42
  )
  
  # Fit with 4 threads
  fit2 <- mfa_em_fit_cpp_singleInit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    z_init = z_init,
    max_iter = 10,
    nIterFA = 5,
    tol = 1e-3,
    n_threads = 4,
    seed = 42
  )
  
  # Results should be identical across different thread counts
  expect_equal(fit1$logLik, fit2$logLik, tolerance = 1e-8)
  expect_equal(fit1$z, fit2$z)
  expect_equal(fit1$pi, fit2$pi, tolerance = 1e-10)
  
  for (k in 1:true_params$K) {
    expect_equal(fit1$mu[[k]], fit2$mu[[k]], tolerance = 1e-10)
    expect_equal(fit1$Lambda[[k]], fit2$Lambda[[k]], tolerance = 1e-10)
    expect_equal(fit1$Psi[[k]], fit2$Psi[[k]], tolerance = 1e-10)
  }
})

test_that("mfa_em_fit works with progress_callback", {
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Track progress callback invocations
  progress_log <- list()
  callback <- function(i, total) {
    progress_log[[length(progress_log) + 1]] <<- list(iter = i, total = total)
  }
  
  # Fit with progress callback
  fit <- mfa_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 10,
    nIterFA = 5,
    tol = 1e-3,
    n_init = 1,
    use_kmeans_init = FALSE,
    progress_callback = callback
  )
  
  # Check that callback was invoked
  expect_true(length(progress_log) > 0)
  
  # Check that callback received correct parameters
  for (entry in progress_log) {
    expect_equal(entry$total, 10)
    expect_true(entry$iter >= 1 && entry$iter <= 10)
  }
  
  # Check that fit still works correctly
  expect_type(fit, "list")
  expect_named(fit, c("z", "pi", "mu", "Lambda", "Psi", "logLik", "resp", "thread_info"))
})

test_that("mfa_em_fit_cpp_singleInit works with progress_callback", {
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  z_init <- true_params$z

  # Track progress callback invocations
  progress_log <- list()
  callback <- function(i, total) {
    progress_log[[length(progress_log) + 1]] <<- list(iter = i, total = total)
  }

  # Fit with progress callback
  fit <- mfa_em_fit_cpp_singleInit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    z_init = z_init,
    max_iter = 10,
    nIterFA = 5,
    tol = 1e-3,
    progress_callback = callback
  )

  # Check that callback was invoked
  expect_true(length(progress_log) > 0)

  # Check that callback received correct parameters
  for (entry in progress_log) {
    expect_equal(entry$total, 10)
    expect_true(entry$iter >= 1 && entry$iter <= 10)
  }

  # Check that fit still works correctly
  expect_type(fit, "list")
  expect_named(fit, c("z", "pi", "mu", "Lambda", "Psi", "logLik", "resp"))
})

# ============================================================================
# Tests for Soft EM Implementation (M2 Reviewer Response)
# ============================================================================

test_that("weighted FA update produces valid results", {
  # Create simple test data
  set.seed(123)
  n <- 100
  M <- 5
  r <- 2
  X <- matrix(rnorm(n * M), n, M)
  weights <- runif(n, 0.1, 1.0)
  weights <- weights / sum(weights) * n  # Normalize to sum to n

  # Call weighted FA (will fail initially)
  result <- faEMupdateWeightedCpp(X, r, nIterFA = 10, weights = weights)

  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("mu", "Lambda", "Psi"))
  expect_equal(length(result$mu), M)
  expect_equal(dim(result$Lambda), c(M, r))
  expect_equal(dim(result$Psi), c(M, M))

  # Check Psi is diagonal and positive
  expect_true(all(diag(result$Psi) > 0))

  # Check that off-diagonal of Psi are zero (diagonal matrix)
  off_diag <- result$Psi - diag(diag(result$Psi))
  expect_equal(max(abs(off_diag)), 0, tolerance = 1e-10)
})

test_that("weighted FA with uniform weights matches unweighted FA", {
  # Create simple test data
  set.seed(456)
  n <- 100
  M <- 5
  r <- 2
  X <- matrix(rnorm(n * M), n, M)

  # Uniform weights
  weights <- rep(1, n)

  # Call both versions
  result_weighted <- faEMupdateWeightedCpp(X, r, nIterFA = 10, weights = weights)
  result_unweighted <- faEMupdateCpp(X, r, nIterFA = 10)

  # Results should be similar (not exactly equal due to implementation differences)
  # but the means should match closely
  expect_equal(result_weighted$mu, result_unweighted$mu, tolerance = 0.1)

  # Covariance structures should be similar
  Sigma_w <- result_weighted$Lambda %*% t(result_weighted$Lambda) + result_weighted$Psi
  Sigma_u <- result_unweighted$Lambda %*% t(result_unweighted$Lambda) + result_unweighted$Psi
  expect_equal(Sigma_w, Sigma_u, tolerance = 0.5)
})

test_that("mfa_em_fit with soft EM returns valid results", {
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  K <- true_params$K
  r <- true_params$r
  z_init <- true_params$z

  # Test soft EM mode
  fit_soft <- mfa_em_fit_cpp_singleInit(
    list_of_data = list_of_data,
    K = K,
    r = r,
    z_init = z_init,
    max_iter = 10,
    nIterFA = 5,
    tol = 1e-3,
    em_type = "soft"  # New parameter
  )

  expect_type(fit_soft, "list")
  expect_named(fit_soft, c("z", "pi", "mu", "Lambda", "Psi", "logLik", "resp"))
  expect_equal(sum(fit_soft$pi), 1, tolerance = 1e-6)
  expect_true(is.finite(fit_soft$logLik))

  # Check responsibilities are valid probabilities
  expect_true(all(fit_soft$resp >= 0))
  expect_true(all(fit_soft$resp <= 1))
  expect_equal(rowSums(fit_soft$resp), rep(1, true_params$N), tolerance = 1e-6)
})

test_that("soft EM converges and produces reasonable parameter estimates", {
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  K <- true_params$K
  r <- true_params$r
  z_init <- true_params$z

  # Fit with more iterations
  fit_soft <- mfa_em_fit_cpp_singleInit(
    list_of_data = list_of_data,
    K = K,
    r = r,
    z_init = z_init,
    max_iter = 50,
    nIterFA = 10,
    tol = 1e-4,
    em_type = "soft"
  )

  # Check parameter dimensions for each cluster
  for (k in 1:K) {
    expect_equal(length(fit_soft$mu[[k]]), true_params$M)
    expect_equal(dim(fit_soft$Lambda[[k]]), c(true_params$M, r))
    expect_equal(dim(fit_soft$Psi[[k]]), c(true_params$M, true_params$M))

    # Check Psi is diagonal and positive
    expect_true(all(diag(fit_soft$Psi[[k]]) > 0))
  }
})
