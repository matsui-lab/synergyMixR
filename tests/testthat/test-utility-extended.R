# Extended Unit Tests for Utility Functions
# Tests for utility.R functions with focus on edge cases and numerical precision

# Load test fixtures
small_data <- readRDS("fixtures/small_test_data.rds")
medium_data <- readRDS("fixtures/medium_test_data.rds")
edge_cases <- readRDS("fixtures/edge_case_data.rds")

test_that("compute_logLik_mfa handles well-separated clusters", {
  # Create data with well-separated clusters
  set.seed(123)
  N <- 10
  T_i <- 20
  M <- 4
  r <- 2
  K <- 2
  
  # Create very different cluster means
  mu_list <- list(
    rep(-2, M),
    rep(2, M)
  )
  
  Lambda_list <- list(
    matrix(rnorm(M * r, sd = 0.5), M, r),
    matrix(rnorm(M * r, sd = 0.5), M, r)
  )
  
  Psi_list <- list(
    diag(0.1, M),
    diag(0.1, M)
  )
  
  z <- rep(1:K, each = N/K)
  
  # Generate data
  list_of_data <- list()
  for (i in 1:N) {
    k <- z[i]
    Z_i <- matrix(rnorm(T_i * r), T_i, r)
    X_i <- Z_i %*% t(Lambda_list[[k]]) + 
           matrix(rep(mu_list[[k]], T_i), T_i, M, byrow = TRUE) +
           matrix(rnorm(T_i * M, 0, sqrt(diag(Psi_list[[k]]))), T_i, M)
    list_of_data[[i]] <- X_i
  }
  
  mfa_fit <- list(
    Lambda = Lambda_list,
    mu = mu_list,
    Psi = Psi_list,
    pi = c(0.5, 0.5)
  )
  
  # Compute log-likelihood
  loglik <- compute_logLik_mfa(list_of_data, mfa_fit)
  
  # Should be finite and negative
  expect_true(is.finite(loglik))
  expect_true(loglik < 0)
})

test_that("compute_logLik_mfa handles single cluster", {
  # Use edge case data with single cluster
  list_of_data <- edge_cases$single_cluster$data
  true_params <- edge_cases$single_cluster$true_params
  
  mfa_fit <- list(
    Lambda = true_params$Lambda,
    mu = true_params$mu,
    Psi = true_params$Psi,
    pi = c(1.0)
  )
  
  # Compute log-likelihood
  loglik <- compute_logLik_mfa(list_of_data, mfa_fit)
  
  # Should be finite and negative
  expect_true(is.finite(loglik))
  expect_true(loglik < 0)
})

test_that("compute_logLik_mfa handles very small noise variance", {
  # Create data with very small noise
  set.seed(456)
  N <- 5
  T_i <- 15
  M <- 3
  r <- 2
  K <- 2
  
  mu_list <- list(rnorm(M), rnorm(M))
  Lambda_list <- list(
    matrix(rnorm(M * r), M, r),
    matrix(rnorm(M * r), M, r)
  )
  
  # Very small noise
  Psi_list <- list(
    diag(0.001, M),
    diag(0.001, M)
  )
  
  z <- rep(1:K, length.out = N)
  
  # Generate data
  list_of_data <- list()
  for (i in 1:N) {
    k <- z[i]
    Z_i <- matrix(rnorm(T_i * r), T_i, r)
    X_i <- Z_i %*% t(Lambda_list[[k]]) + 
           matrix(rep(mu_list[[k]], T_i), T_i, M, byrow = TRUE) +
           matrix(rnorm(T_i * M, 0, sqrt(diag(Psi_list[[k]]))), T_i, M)
    list_of_data[[i]] <- X_i
  }
  
  mfa_fit <- list(
    Lambda = Lambda_list,
    mu = mu_list,
    Psi = Psi_list,
    pi = c(0.5, 0.5)
  )
  
  # Compute log-likelihood
  loglik <- compute_logLik_mfa(list_of_data, mfa_fit)
  
  # Should be finite
  expect_true(is.finite(loglik))
})

test_that("compute_logLik_mfa is consistent with BIC computation", {
  # Use small test data
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  mfa_fit <- list(
    Lambda = true_params$Lambda,
    mu = true_params$mu,
    Psi = true_params$Psi,
    pi = rep(1/true_params$K, true_params$K)
  )
  
  # Compute log-likelihood
  loglik <- compute_logLik_mfa(list_of_data, mfa_fit)
  
  # Compute BIC
  bic <- compute_BIC_mfa(loglik, true_params$K, true_params$r, 
                         true_params$M, true_params$N)
  
  # BIC should be finite and positive (since it's -2*loglik + penalty)
  expect_true(is.finite(bic))
  expect_true(bic > 0)
  
  # Check BIC formula (matches compute_BIC_mfa implementation)
  # For each cluster: M*r (loadings) + r*(r-1)/2 (factor covariances) + M (unique variances) + M (means)
  # Plus (K-1) for mixing proportions
  K <- true_params$K
  M <- true_params$M
  r <- true_params$r
  param_count <- K * (M*r + r*(r-1)/2 + M + M) + (K - 1)
  expected_bic <- -2 * loglik + param_count * log(true_params$N)
  expect_equal(bic, expected_bic, tolerance = 1e-10)
})

test_that("compute_BIC_mfa increases with model complexity", {
  # For fixed data, BIC should increase with K and r (due to penalty term)
  N <- 10
  M <- 4
  
  # Compute BIC for different model complexities (same loglik for comparison)
  loglik <- -100
  
  bic_K2_r2 <- compute_BIC_mfa(loglik, K = 2, r = 2, M = M, N = N)
  bic_K3_r2 <- compute_BIC_mfa(loglik, K = 3, r = 2, M = M, N = N)
  bic_K2_r3 <- compute_BIC_mfa(loglik, K = 2, r = 3, M = M, N = N)
  
  # Higher K or r should give higher BIC (more penalty)
  expect_true(bic_K3_r2 > bic_K2_r2)
  expect_true(bic_K2_r3 > bic_K2_r2)
})

test_that("compute_logLik_mpca handles well-separated clusters", {
  # Create MPCA data with well-separated clusters
  set.seed(789)
  N <- 10
  T_i <- 20
  M <- 4
  r <- 2
  K <- 2
  
  mu_list <- list(
    rep(-2, M),
    rep(2, M)
  )
  
  P_list <- list(
    qr.Q(qr(matrix(rnorm(M * r), M, r))),
    qr.Q(qr(matrix(rnorm(M * r), M, r)))
  )
  
  sigma2vec <- c(0.5, 0.5)
  z <- rep(1:K, each = N/K)
  
  # Generate data
  list_of_data <- list()
  for (i in 1:N) {
    k <- z[i]
    Z_i <- matrix(rnorm(T_i * r), T_i, r)
    X_i <- Z_i %*% t(P_list[[k]]) + 
           matrix(rep(mu_list[[k]], T_i), T_i, M, byrow = TRUE) +
           matrix(rnorm(T_i * M, 0, sqrt(sigma2vec[k])), T_i, M)
    list_of_data[[i]] <- X_i
  }
  
  mpca_fit <- list(
    P = P_list,
    mu = mu_list,
    sigma2 = sigma2vec,
    pi = c(0.5, 0.5)
  )
  
  # Compute log-likelihood
  loglik <- compute_logLik_mpca(list_of_data, mpca_fit)
  
  # Should be finite and negative
  expect_true(is.finite(loglik))
  expect_true(loglik < 0)
})

test_that("compute_logLik_mpca handles single cluster", {
  # Create single cluster MPCA data
  set.seed(111)
  N <- 5
  T_i <- 15
  M <- 4
  r <- 2
  
  P <- qr.Q(qr(matrix(rnorm(M * r), M, r)))
  mu <- rnorm(M)
  sigma2 <- 0.5
  
  # Generate data
  list_of_data <- list()
  for (i in 1:N) {
    Z_i <- matrix(rnorm(T_i * r), T_i, r)
    X_i <- Z_i %*% t(P) + 
           matrix(rep(mu, T_i), T_i, M, byrow = TRUE) +
           matrix(rnorm(T_i * M, 0, sqrt(sigma2)), T_i, M)
    list_of_data[[i]] <- X_i
  }
  
  mpca_fit <- list(
    P = list(P),
    mu = list(mu),
    sigma2 = c(sigma2),
    pi = c(1.0)
  )
  
  # Compute log-likelihood
  loglik <- compute_logLik_mpca(list_of_data, mpca_fit)
  
  # Should be finite and negative
  expect_true(is.finite(loglik))
  expect_true(loglik < 0)
})

test_that("compute_BIC_mpca increases with model complexity", {
  # For fixed data, BIC should increase with K and r
  N <- 10
  M <- 4
  
  loglik <- -100
  
  bic_K2_r2 <- compute_BIC_mpca(loglik, K = 2, r = 2, M = M, N = N)
  bic_K3_r2 <- compute_BIC_mpca(loglik, K = 3, r = 2, M = M, N = N)
  bic_K2_r3 <- compute_BIC_mpca(loglik, K = 2, r = 3, M = M, N = N)
  
  # Higher K or r should give higher BIC
  expect_true(bic_K3_r2 > bic_K2_r2)
  expect_true(bic_K2_r3 > bic_K2_r2)
})

test_that("compute_global_vaf_mfa handles zero variance case", {
  # Create data with zero variance (constant values)
  set.seed(222)
  N <- 5
  T_i <- 10
  M <- 3
  r <- 2
  K <- 2
  
  # Create constant data (zero variance)
  list_of_data <- list()
  for (i in 1:N) {
    X_i <- matrix(1.0, T_i, M)
    list_of_data[[i]] <- X_i
  }
  
  z <- rep(1:K, length.out = N)
  mu_list <- list(rep(1.0, M), rep(1.0, M))
  Lambda_list <- list(
    matrix(0, M, r),
    matrix(0, M, r)
  )
  Psi_list <- list(
    diag(0.01, M),
    diag(0.01, M)
  )
  
  mfa_fit <- list(
    z = z,
    mu = mu_list,
    Lambda = Lambda_list,
    Psi = Psi_list
  )
  
  # Compute VAF (should handle zero variance gracefully)
  vaf <- compute_global_vaf_mfa(list_of_data, mfa_fit)
  
  # VAF should be defined (could be NaN or 1 depending on implementation)
  expect_true(is.numeric(vaf))
})

test_that("compute_global_vaf_mfa handles near-perfect fit", {
  # Create data from model with very small noise
  set.seed(333)
  N <- 5
  T_i <- 15
  M <- 4
  r <- 2
  K <- 2
  
  z <- rep(1:K, length.out = N)
  mu_list <- list(rnorm(M), rnorm(M))
  Lambda_list <- list(
    matrix(rnorm(M * r), M, r),
    matrix(rnorm(M * r), M, r)
  )
  Psi_list <- list(
    diag(0.001, M),
    diag(0.001, M)
  )
  
  # Generate data with very small noise
  list_of_data <- list()
  for (i in 1:N) {
    k <- z[i]
    Z_i <- matrix(rnorm(T_i * r), T_i, r)
    X_i <- Z_i %*% t(Lambda_list[[k]]) + 
           matrix(rep(mu_list[[k]], T_i), T_i, M, byrow = TRUE) +
           matrix(rnorm(T_i * M, 0, 0.01), T_i, M)
    list_of_data[[i]] <- X_i
  }
  
  mfa_fit <- list(
    z = z,
    mu = mu_list,
    Lambda = Lambda_list,
    Psi = Psi_list
  )
  
  # Compute VAF
  vaf <- compute_global_vaf_mfa(list_of_data, mfa_fit)
  
  # Should be very high (close to 1)
  expect_true(vaf > 0.95)
  expect_true(vaf <= 1.0)
})

test_that("compute_global_vaf_mfa handles poor fit", {
  # Create data and fit with mismatched parameters
  set.seed(444)
  N <- 5
  T_i <- 15
  M <- 4
  r <- 2
  K <- 2
  
  z_true <- rep(1:K, length.out = N)
  mu_true <- list(rnorm(M), rnorm(M))
  Lambda_true <- list(
    matrix(rnorm(M * r), M, r),
    matrix(rnorm(M * r), M, r)
  )
  
  # Generate data
  list_of_data <- list()
  for (i in 1:N) {
    k <- z_true[i]
    Z_i <- matrix(rnorm(T_i * r), T_i, r)
    X_i <- Z_i %*% t(Lambda_true[[k]]) + 
           matrix(rep(mu_true[[k]], T_i), T_i, M, byrow = TRUE) +
           matrix(rnorm(T_i * M, 0, 0.5), T_i, M)
    list_of_data[[i]] <- X_i
  }
  
  # Use wrong parameters (random)
  mfa_fit <- list(
    z = z_true,
    mu = list(rnorm(M), rnorm(M)),
    Lambda = list(
      matrix(rnorm(M * r), M, r),
      matrix(rnorm(M * r), M, r)
    ),
    Psi = list(diag(1, M), diag(1, M))
  )
  
  # Compute VAF
  vaf <- compute_global_vaf_mfa(list_of_data, mfa_fit)
  
  # Should be lower (but still between 0 and 1)
  expect_true(vaf >= 0)
  expect_true(vaf <= 1)
})

test_that("compute_global_vaf_mpca handles edge cases", {
  # Use edge case data
  list_of_data <- edge_cases$single_cluster$data
  true_params <- edge_cases$single_cluster$true_params
  
  # Create MPCA structure
  W_list <- list()
  for (k in 1:true_params$K) {
    W_list[[k]] <- true_params$Lambda[[k]]
  }
  
  mpca_fit <- list(
    z = true_params$z,
    mu = true_params$mu,
    W = W_list,
    sigma2 = rep(0.5, true_params$K)
  )
  
  # Compute VAF
  vaf <- compute_global_vaf_mpca(list_of_data, mpca_fit)
  
  # Should be between 0 and 1
  expect_true(vaf >= 0)
  expect_true(vaf <= 1)
})

test_that("compute_global_vaf_mpca handles short time series", {
  # Use short time series edge case
  list_of_data <- edge_cases$short_timeseries$data
  true_params <- edge_cases$short_timeseries$true_params
  
  # Create MPCA structure
  W_list <- list()
  for (k in 1:true_params$K) {
    W_list[[k]] <- true_params$Lambda[[k]]
  }
  
  mpca_fit <- list(
    z = true_params$z,
    mu = true_params$mu,
    W = W_list,
    sigma2 = rep(0.5, true_params$K)
  )
  
  # Compute VAF
  vaf <- compute_global_vaf_mpca(list_of_data, mpca_fit)
  
  # Should be between 0 and 1
  expect_true(vaf >= 0)
  expect_true(vaf <= 1)
})

test_that("compute_factor_scores_optimized handles extreme values", {
  # Test with very large values
  set.seed(555)
  T_i <- 10
  M <- 4
  r <- 2
  
  X <- matrix(rnorm(T_i * M, mean = 100, sd = 10), T_i, M)
  mu <- rnorm(M, mean = 100, sd = 5)
  Lambda <- matrix(rnorm(M * r), M, r)
  Psi <- diag(runif(M, 0.1, 0.5))
  
  # Should handle large values
  Z <- compute_factor_scores_optimized(X, mu, Lambda, Psi)
  
  expect_equal(dim(Z), c(T_i, r))
  expect_true(all(is.finite(Z)))
})

test_that("compute_factor_scores_optimized handles very small variance", {
  # Test with very small Psi values
  set.seed(666)
  T_i <- 10
  M <- 4
  r <- 2
  
  X <- matrix(rnorm(T_i * M), T_i, M)
  mu <- rnorm(M)
  Lambda <- matrix(rnorm(M * r), M, r)
  Psi <- diag(rep(0.001, M))
  
  # Should handle small variance
  Z <- compute_factor_scores_optimized(X, mu, Lambda, Psi)
  
  expect_equal(dim(Z), c(T_i, r))
  expect_true(all(is.finite(Z)))
})

test_that("compute_factor_scores_mpca_single handles edge cases", {
  # Test with single time point
  set.seed(777)
  M <- 4
  r <- 2
  
  X <- matrix(rnorm(M), 1, M)
  mu <- rnorm(M)
  W <- matrix(rnorm(M * r), M, r)
  sigma2 <- 0.5
  
  Z <- compute_factor_scores_mpca_single(X, mu, W, sigma2)
  
  expect_equal(dim(Z), c(1, r))
  expect_true(all(is.finite(Z)))
})

test_that("compute_factor_scores_mpca_single handles very small sigma2", {
  # Test with very small noise variance
  set.seed(888)
  T_i <- 10
  M <- 4
  r <- 2
  
  X <- matrix(rnorm(T_i * M), T_i, M)
  mu <- rnorm(M)
  W <- matrix(rnorm(M * r), M, r)
  sigma2 <- 0.001
  
  Z <- compute_factor_scores_mpca_single(X, mu, W, sigma2)
  
  expect_equal(dim(Z), c(T_i, r))
  expect_true(all(is.finite(Z)))
})

test_that("get_model_by_K_r_mpca works correctly", {
  # Create mock selection object for MPCA
  model1 <- list(params = "mpca_model1")
  model2 <- list(params = "mpca_model2")
  
  selection_obj <- list(
    all_models = list(
      list(K = 2, r = 2, model = model1),
      list(K = 3, r = 2, model = model2)
    )
  )
  
  # Test retrieval
  retrieved1 <- get_model_by_K_r_mpca(selection_obj, K_target = 2, r_target = 2)
  expect_equal(retrieved1, model1)
  
  retrieved2 <- get_model_by_K_r_mpca(selection_obj, K_target = 3, r_target = 2)
  expect_equal(retrieved2, model2)
  
  # Test non-existent combination
  expect_message(
    retrieved_none <- get_model_by_K_r_mpca(selection_obj, K_target = 4, r_target = 4),
    "No PCA model found"
  )
  expect_null(retrieved_none)
})

test_that("posthoc evaluations handle unbalanced clusters", {
  # Use unbalanced cluster edge case
  list_of_data <- edge_cases$unbalanced$data
  true_params <- edge_cases$unbalanced$true_params
  
  model1 <- list(
    z = true_params$z,
    mu = true_params$mu,
    Lambda = true_params$Lambda,
    Psi = true_params$Psi
  )
  
  selection_obj <- list(
    all_models = list(
      list(K = 2, r = 2, model = model1, logLik = -100, BIC = 250)
    )
  )
  
  # Compute post-hoc evaluation
  df_posthoc <- posthoc_mfa_evaluation(list_of_data, selection_obj)
  
  # Check that it handles unbalanced clusters
  expect_s3_class(df_posthoc, "data.frame")
  expect_equal(nrow(df_posthoc), 1)
  expect_true(df_posthoc$minClusterSize < df_posthoc$maxClusterSize)
})

test_that("VAF computation is numerically stable", {
  # Test with data that might cause numerical issues
  set.seed(999)
  N <- 5
  T_i <- 20
  M <- 4
  r <- 2
  K <- 2
  
  z <- rep(1:K, length.out = N)
  mu_list <- list(rnorm(M, mean = 0, sd = 0.1), rnorm(M, mean = 0, sd = 0.1))
  Lambda_list <- list(
    matrix(rnorm(M * r, sd = 0.1), M, r),
    matrix(rnorm(M * r, sd = 0.1), M, r)
  )
  Psi_list <- list(
    diag(runif(M, 0.01, 0.1)),
    diag(runif(M, 0.01, 0.1))
  )
  
  # Generate data
  list_of_data <- list()
  for (i in 1:N) {
    k <- z[i]
    Z_i <- matrix(rnorm(T_i * r, sd = 0.1), T_i, r)
    X_i <- Z_i %*% t(Lambda_list[[k]]) + 
           matrix(rep(mu_list[[k]], T_i), T_i, M, byrow = TRUE) +
           matrix(rnorm(T_i * M, 0, sqrt(diag(Psi_list[[k]]))), T_i, M)
    list_of_data[[i]] <- X_i
  }
  
  mfa_fit <- list(
    z = z,
    mu = mu_list,
    Lambda = Lambda_list,
    Psi = Psi_list
  )
  
  # Compute VAF multiple times
  vaf1 <- compute_global_vaf_mfa(list_of_data, mfa_fit)
  vaf2 <- compute_global_vaf_mfa(list_of_data, mfa_fit)
  
  # Should be identical (numerically stable)
  expect_equal(vaf1, vaf2)
  expect_true(is.finite(vaf1))
})
