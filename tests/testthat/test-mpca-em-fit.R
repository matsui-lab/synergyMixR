# Unit Tests for MPCA EM Fitting Functions
# Tests for mpca_em_fit.R functions

# Load test fixtures
small_data <- readRDS("fixtures/small_test_data.rds")
medium_data <- readRDS("fixtures/medium_test_data.rds")
edge_cases <- readRDS("fixtures/edge_case_data.rds")

test_that("mixture_pca_em_fit_cpp_singleInit works with valid inputs", {
  # Use small test data
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Test with known parameters
  K <- true_params$K
  r <- true_params$r
  z_init <- true_params$z
  
  # Fit model
  fit <- mixture_pca_em_fit_cpp_singleInit(
    list_of_data = list_of_data,
    K = K,
    r = r,
    max_iter = 10,
    nIterPCA = 5,
    tol = 1e-3,
    method = "EM",
    z_init = z_init
  )
  
  # Check output structure
  expect_type(fit, "list")
  expect_named(fit, c("z", "pi", "mu", "W", "P", "D", "Psi", "logLik", "sigma2", "resp"))
  
  # Check dimensions
  expect_equal(length(fit$z), true_params$N)
  expect_equal(length(fit$pi), K)
  expect_equal(length(fit$mu), K)
  expect_equal(length(fit$W), K)
  expect_equal(length(fit$P), K)
  expect_equal(length(fit$D), K)
  expect_equal(length(fit$Psi), K)
  expect_equal(length(fit$sigma2), K)
  expect_equal(dim(fit$resp), c(true_params$N, K))
  
  # Check cluster assignments are valid
  expect_true(all(fit$z %in% 1:K))
  
  # Check mixing proportions sum to 1
  expect_equal(sum(fit$pi), 1, tolerance = 1e-6)
  expect_true(all(fit$pi > 0))
  
  # Check parameter dimensions
  for (k in 1:K) {
    expect_equal(length(fit$mu[[k]]), true_params$M)
    expect_equal(dim(fit$W[[k]]), c(true_params$M, r))
    expect_equal(dim(fit$P[[k]]), c(true_params$M, r))
    expect_equal(dim(fit$D[[k]]), c(r, r))
    expect_equal(dim(fit$Psi[[k]]), c(true_params$M, true_params$M))
    
    # Check sigma2 is positive
    expect_true(fit$sigma2[k] > 0)
    
    # Check P has orthonormal columns (approximately)
    P_k <- fit$P[[k]]
    PtP <- t(P_k) %*% P_k
    expect_equal(PtP, diag(r), tolerance = 1e-4)
    
    # Check D is diagonal
    D_k <- fit$D[[k]]
    expect_equal(D_k, diag(diag(D_k)), tolerance = 1e-10)
    expect_true(all(diag(D_k) > 0))
    
    # Check Psi is diagonal with sigma2 on diagonal
    expect_equal(fit$Psi[[k]], diag(fit$sigma2[k], true_params$M), tolerance = 1e-10)
  }
  
  # Check log-likelihood is finite
  expect_true(is.finite(fit$logLik))
  expect_false(is.na(fit$logLik))
  
  # Check responsibilities are valid probabilities
  expect_true(all(fit$resp >= 0))
  expect_true(all(fit$resp <= 1))
  expect_equal(rowSums(fit$resp), rep(1, true_params$N), tolerance = 1e-6)
})

test_that("mixture_pca_em_fit_cpp_singleInit works with closed-form method", {
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Fit with closed-form method
  fit_closed <- mixture_pca_em_fit_cpp_singleInit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 10,
    nIterPCA = 5,
    tol = 1e-3,
    method = "closed_form",
    z_init = true_params$z
  )
  
  # Check output structure
  expect_type(fit_closed, "list")
  expect_true(is.finite(fit_closed$logLik))
  
  # Compare with EM method
  fit_em <- mixture_pca_em_fit_cpp_singleInit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 10,
    nIterPCA = 5,
    tol = 1e-3,
    method = "EM",
    z_init = true_params$z
  )
  
  # Both should produce valid results
  expect_true(is.finite(fit_em$logLik))
  expect_equal(length(fit_closed$z), length(fit_em$z))
})

test_that("mixture_pca_em_fit_cpp_singleInit handles edge cases", {
  # Test with single cluster
  single_cluster_data <- edge_cases$single_cluster$data
  
  fit_single <- mixture_pca_em_fit_cpp_singleInit(
    list_of_data = single_cluster_data,
    K = 1,
    r = 2,
    max_iter = 5,
    nIterPCA = 3,
    tol = 1e-3,
    method = "EM",
    z_init = rep(1, length(single_cluster_data))
  )
  
  expect_equal(length(fit_single$pi), 1)
  expect_equal(fit_single$pi[1], 1)
  expect_true(all(fit_single$z == 1))
  
  # Test with single factor
  single_factor_data <- edge_cases$single_factor$data
  
  fit_single_factor <- mixture_pca_em_fit_cpp_singleInit(
    list_of_data = single_factor_data,
    K = 2,
    r = 1,
    max_iter = 5,
    nIterPCA = 3,
    tol = 1e-3,
    method = "EM",
    z_init = edge_cases$single_factor$true_params$z
  )
  
  expect_equal(ncol(fit_single_factor$W[[1]]), 1)
  expect_equal(ncol(fit_single_factor$W[[2]]), 1)
  expect_equal(ncol(fit_single_factor$P[[1]]), 1)
  expect_equal(ncol(fit_single_factor$P[[2]]), 1)
})

test_that("mixture_pca_em_fit works with default parameters", {
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Test single-run mode (no multi-init)
  fit_single <- mixture_pca_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 10,
    nIterPCA = 5,
    tol = 1e-3,
    method = "EM",
    n_init = 1,
    use_kmeans_init = FALSE,
    seed = 123
  )
  
  # Check output structure
  expect_type(fit_single, "list")
  expect_named(fit_single, c("z", "pi", "mu", "P", "D", "Psi", "logLik", "sigma2", "resp"))
  
  # Check basic properties
  expect_equal(length(fit_single$z), true_params$N)
  expect_equal(length(fit_single$pi), true_params$K)
  expect_equal(sum(fit_single$pi), 1, tolerance = 1e-6)
  expect_true(is.finite(fit_single$logLik))
})

test_that("mixture_pca_em_fit works with multi-initialization", {
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Test multi-init mode
  fit_multi <- mixture_pca_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 5,
    nIterPCA = 3,
    tol = 1e-3,
    method = "EM",
    n_init = 3,
    use_kmeans_init = TRUE,
    subject_rdim_for_kmeans = 2,
    mc_cores = 1,
    seed = 456
  )
  
  # Check output structure (W may or may not be present depending on implementation)
  expect_type(fit_multi, "list")
  expect_true(all(c("z", "pi", "mu", "P", "D", "Psi", "logLik", "sigma2", "resp") %in% names(fit_multi)))
  
  # Check basic properties
  expect_equal(length(fit_multi$z), true_params$N)
  expect_equal(length(fit_multi$pi), true_params$K)
  expect_equal(sum(fit_multi$pi), 1, tolerance = 1e-6)
  expect_true(is.finite(fit_multi$logLik))
})

test_that("mixture_pca_em_fit is deterministic with fixed seed", {
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Fit same model twice with same seed
  fit1 <- mixture_pca_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 5,
    nIterPCA = 3,
    tol = 1e-3,
    method = "EM",
    n_init = 2,
    use_kmeans_init = TRUE,
    mc_cores = 1,
    seed = 789
  )
  
  fit2 <- mixture_pca_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 5,
    nIterPCA = 3,
    tol = 1e-3,
    method = "EM",
    n_init = 2,
    use_kmeans_init = TRUE,
    mc_cores = 1,
    seed = 789
  )
  
  # Results should be identical
  expect_equal(fit1$z, fit2$z)
  expect_equal(fit1$pi, fit2$pi, tolerance = 1e-10)
  expect_equal(fit1$logLik, fit2$logLik, tolerance = 1e-10)
  expect_equal(fit1$sigma2, fit2$sigma2, tolerance = 1e-10)
  
  for (k in 1:true_params$K) {
    expect_equal(fit1$mu[[k]], fit2$mu[[k]], tolerance = 1e-10)
    expect_equal(fit1$P[[k]], fit2$P[[k]], tolerance = 1e-10)
    expect_equal(fit1$D[[k]], fit2$D[[k]], tolerance = 1e-10)
  }
})

test_that("mixture_pca_em_fit handles Windows parallel fallback", {
  # Skip this test - cannot mock .Platform$OS.type reliably
  skip("Cannot reliably mock .Platform$OS.type for testing")
})

test_that("mixture_pca_em_fit convergence properties", {
  list_of_data <- medium_data$data
  true_params <- medium_data$true_params
  
  # Test that more iterations generally improve log-likelihood
  fit_short <- mixture_pca_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 3,
    nIterPCA = 2,
    method = "EM",
    n_init = 1,
    use_kmeans_init = FALSE,
    seed = 111
  )
  
  fit_long <- mixture_pca_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 15,
    nIterPCA = 8,
    method = "EM",
    n_init = 1,
    use_kmeans_init = FALSE,
    seed = 111
  )
  
  # Longer run should have equal or better log-likelihood
  expect_true(fit_long$logLik >= fit_short$logLik - 1e-6)
})

test_that("mixture_pca_em_fit parameter recovery on synthetic data", {
  skip_if_not_installed("mclust")

  # Use small data with known true parameters
  list_of_data <- small_data$data
  true_params <- small_data$true_params

  # Fit model with correct K and r
  fit <- mixture_pca_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 50,
    nIterPCA = 20,
    tol = 1e-4,
    method = "EM",
    n_init = 5,
    use_kmeans_init = TRUE,
    mc_cores = 1,
    seed = 123
  )
  
  # Check that we recover reasonable cluster assignments
  ari_score <- mclust::adjustedRandIndex(fit$z, true_params$z)
  expect_true(ari_score > 0.5,
              info = paste("ARI score:", round(ari_score, 3)))
  
  # Check that parameters are in reasonable ranges
  for (k in 1:true_params$K) {
    # Mean parameters should be finite
    expect_true(all(is.finite(fit$mu[[k]])))
    
    # Loading matrices should be finite
    expect_true(all(is.finite(fit$P[[k]])))
    expect_true(all(is.finite(fit$D[[k]])))
    
    # Noise variance should be positive and finite
    expect_true(fit$sigma2[k] > 0)
    expect_true(is.finite(fit$sigma2[k]))
    
    # P should have orthonormal columns
    P_k <- fit$P[[k]]
    PtP <- t(P_k) %*% P_k
    expect_equal(PtP, diag(true_params$r), tolerance = 1e-3)
  }
  
  # Log-likelihood should be reasonable
  expect_true(fit$logLik > -1e6)
})

test_that("mixture_pca_em_fit W matrix consistency", {
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Fit model
  fit <- mixture_pca_em_fit_cpp_singleInit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 10,
    nIterPCA = 5,
    tol = 1e-3,
    method = "EM",
    z_init = true_params$z
  )
  
  # Check that W = P * D (approximately)
  for (k in 1:true_params$K) {
    W_reconstructed <- fit$P[[k]] %*% fit$D[[k]]
    expect_equal(fit$W[[k]], W_reconstructed, tolerance = 1e-8)
  }
})

test_that("mixture_pca_em_fit produces reproducible results across different thread counts with cpp_seed", {
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Test with single initialization to isolate C++ RNG determinism
  # Fit with 2 threads
  fit1 <- mixture_pca_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 10,
    nIterPCA = 5,
    tol = 1e-3,
    method = "EM",
    n_init = 1,
    use_kmeans_init = FALSE,
    n_threads = 2,
    cpp_seed = 42
  )
  
  # Fit with 4 threads
  fit2 <- mixture_pca_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 10,
    nIterPCA = 5,
    tol = 1e-3,
    method = "EM",
    n_init = 1,
    use_kmeans_init = FALSE,
    n_threads = 4,
    cpp_seed = 42
  )
  
  # Results should be identical across different thread counts
  expect_equal(fit1$logLik, fit2$logLik, tolerance = 1e-8)
  expect_equal(fit1$z, fit2$z)
  expect_equal(fit1$pi, fit2$pi, tolerance = 1e-10)
  expect_equal(fit1$sigma2, fit2$sigma2, tolerance = 1e-10)
  
  for (k in 1:true_params$K) {
    expect_equal(fit1$mu[[k]], fit2$mu[[k]], tolerance = 1e-10)
    expect_equal(fit1$P[[k]], fit2$P[[k]], tolerance = 1e-10)
    expect_equal(fit1$D[[k]], fit2$D[[k]], tolerance = 1e-10)
  }
})

test_that("mixture_pca_em_fit_cpp_singleInit produces reproducible results with seed parameter", {
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  z_init <- true_params$z
  
  # Fit with 2 threads
  fit1 <- mixture_pca_em_fit_cpp_singleInit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 10,
    nIterPCA = 5,
    tol = 1e-3,
    method = "EM",
    z_init = z_init,
    n_threads = 2,
    seed = 42
  )
  
  # Fit with 4 threads
  fit2 <- mixture_pca_em_fit_cpp_singleInit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 10,
    nIterPCA = 5,
    tol = 1e-3,
    method = "EM",
    z_init = z_init,
    n_threads = 4,
    seed = 42
  )
  
  # Results should be identical across different thread counts
  expect_equal(fit1$logLik, fit2$logLik, tolerance = 1e-8)
  expect_equal(fit1$z, fit2$z)
  expect_equal(fit1$pi, fit2$pi, tolerance = 1e-10)
  expect_equal(fit1$sigma2, fit2$sigma2, tolerance = 1e-10)
  
  for (k in 1:true_params$K) {
    expect_equal(fit1$mu[[k]], fit2$mu[[k]], tolerance = 1e-10)
    expect_equal(fit1$P[[k]], fit2$P[[k]], tolerance = 1e-10)
    expect_equal(fit1$D[[k]], fit2$D[[k]], tolerance = 1e-10)
    expect_equal(fit1$W[[k]], fit2$W[[k]], tolerance = 1e-10)
  }
})

test_that("mixture_pca_em_fit works with progress_callback", {
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Track progress callback invocations
  progress_log <- list()
  callback <- function(i, total) {
    progress_log[[length(progress_log) + 1]] <<- list(iter = i, total = total)
  }
  
  # Fit with progress callback
  fit <- mixture_pca_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 10,
    nIterPCA = 5,
    tol = 1e-3,
    method = "EM",
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
  expect_true(all(c("z", "pi", "mu", "P", "D", "sigma2", "logLik", "resp") %in% names(fit)))
})

test_that("mixture_pca_em_fit_cpp_singleInit works with progress_callback", {
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  z_init <- true_params$z
  
  # Track progress callback invocations
  progress_log <- list()
  callback <- function(i, total) {
    progress_log[[length(progress_log) + 1]] <<- list(iter = i, total = total)
  }
  
  # Fit with progress callback
  fit <- mixture_pca_em_fit_cpp_singleInit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 10,
    nIterPCA = 5,
    tol = 1e-3,
    method = "EM",
    z_init = z_init,
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
  expect_true(all(c("z", "pi", "mu", "P", "D", "W", "sigma2", "logLik", "resp") %in% names(fit)))
})

# ============================================================================
# Tests for Soft EM Implementation (M2 Reviewer Response)
# ============================================================================

test_that("mpca_em_fit with soft EM returns valid results", {
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  K <- true_params$K
  r <- true_params$r
  z_init <- true_params$z

  fit_soft <- mixture_pca_em_fit_cpp_singleInit(
    list_of_data = list_of_data,
    K = K,
    r = r,
    z_init = z_init,
    max_iter = 10,
    nIterPCA = 5,
    tol = 1e-3,
    em_type = "soft"
  )

  expect_type(fit_soft, "list")
  expect_true(all(c("z", "pi", "mu", "W", "sigma2", "logLik", "resp") %in% names(fit_soft)))
  expect_equal(sum(fit_soft$pi), 1, tolerance = 1e-6)
  expect_true(is.finite(fit_soft$logLik))

  # Check responsibilities are valid probabilities
  expect_true(all(fit_soft$resp >= 0))
  expect_true(all(fit_soft$resp <= 1))
  expect_equal(rowSums(fit_soft$resp), rep(1, true_params$N), tolerance = 1e-6)
})

test_that("mpca soft EM converges and produces reasonable parameter estimates", {
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  K <- true_params$K
  r <- true_params$r
  z_init <- true_params$z

  # Fit with more iterations
  fit_soft <- mixture_pca_em_fit_cpp_singleInit(
    list_of_data = list_of_data,
    K = K,
    r = r,
    z_init = z_init,
    max_iter = 50,
    nIterPCA = 10,
    tol = 1e-4,
    em_type = "soft"
  )

  # Check parameter dimensions for each cluster
  for (k in 1:K) {
    expect_equal(length(fit_soft$mu[[k]]), true_params$M)
    expect_equal(dim(fit_soft$W[[k]]), c(true_params$M, r))

    # Check sigma2 is positive
    expect_true(fit_soft$sigma2[k] > 0)
  }
})
