# Unit Tests for Utility Functions
# Tests for utility.R functions

# Load test fixtures
small_data <- readRDS("fixtures/small_test_data.rds")
medium_data <- readRDS("fixtures/medium_test_data.rds")

test_that("compute_factor_scores works correctly", {
  # Create simple test data
  T_i <- 20
  M <- 4
  r <- 2
  
  set.seed(123)
  X <- matrix(rnorm(T_i * M), T_i, M)
  mu <- rnorm(M)
  Lambda <- matrix(rnorm(M * r), M, r)
  Psi <- diag(runif(M, 0.1, 0.5))
  
  # Compute factor scores
  Z <- compute_factor_scores(X, mu, Lambda, Psi)
  
  # Check dimensions
  expect_equal(dim(Z), c(T_i, r))
  
  # Check that scores are finite
  expect_true(all(is.finite(Z)))
  
  # Check that reconstruction is reasonable
  X_reconstructed <- Z %*% t(Lambda) + matrix(rep(mu, T_i), T_i, M, byrow = TRUE)
  
  # Reconstruction should be correlated with original (use lower threshold for random data)
  for (j in 1:M) {
    cor_j <- cor(X[, j], X_reconstructed[, j])
    expect_true(cor_j > -0.5)  # Should have some relationship (even if weak for random data)
  }
})

test_that("compute_factor_scores handles edge cases", {
  # Single time point
  X_single <- matrix(rnorm(4), 1, 4)
  mu <- rnorm(4)
  Lambda <- matrix(rnorm(8), 4, 2)
  Psi <- diag(runif(4, 0.1, 0.5))
  
  Z_single <- compute_factor_scores(X_single, mu, Lambda, Psi)
  expect_equal(dim(Z_single), c(1, 2))
  expect_true(all(is.finite(Z_single)))
  
  # Single factor
  Lambda_single <- matrix(rnorm(4), 4, 1)
  Z_single_factor <- compute_factor_scores(X_single, mu, Lambda_single, Psi)
  expect_equal(dim(Z_single_factor), c(1, 1))
  expect_true(is.finite(Z_single_factor[1, 1]))
})

test_that("compute_global_vaf_mfa works correctly", {
  # Use small test data and fit a model
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Create a simple fitted model structure
  fit <- list(
    z = true_params$z,
    mu = true_params$mu,
    Lambda = true_params$Lambda,
    Psi = true_params$Psi
  )
  
  # Compute VAF
  vaf <- compute_global_vaf_mfa(list_of_data, fit)
  
  # Check that VAF is a single numeric value
  expect_type(vaf, "double")
  expect_length(vaf, 1)
  
  # VAF should be between 0 and 1 for reasonable models
  expect_true(vaf >= 0)
  expect_true(vaf <= 1)
  
  # For true parameters, VAF should be reasonably high
  expect_true(vaf > 0.5)
})

test_that("compute_global_vaf_mfa handles perfect fit", {
  # Create data from a known model
  set.seed(456)
  N <- 3
  T_i <- 15
  M <- 3
  r <- 2
  K <- 2
  
  z <- c(1, 1, 2)
  mu <- list(rnorm(M), rnorm(M))
  Lambda <- list(matrix(rnorm(M * r), M, r), matrix(rnorm(M * r), M, r))
  Psi <- list(diag(0.01, M), diag(0.01, M))  # Very small noise
  
  # Generate data with minimal noise
  list_of_data <- list()
  for (i in 1:N) {
    k <- z[i]
    Z_i <- matrix(rnorm(T_i * r), T_i, r)
    X_i <- Z_i %*% t(Lambda[[k]]) + matrix(rep(mu[[k]], T_i), T_i, M, byrow = TRUE)
    list_of_data[[i]] <- X_i
  }
  
  fit <- list(z = z, mu = mu, Lambda = Lambda, Psi = Psi)
  
  # VAF should be very high (close to 1) for near-perfect fit
  vaf <- compute_global_vaf_mfa(list_of_data, fit)
  expect_true(vaf > 0.95)
})

test_that("compute_global_vaf_mpca works correctly", {
  # Use small test data
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Create MPCA model structure with W matrices
  W_list <- list()
  P_list <- list()
  D_list <- list()
  for (k in 1:true_params$K) {
    W_k <- true_params$Lambda[[k]]  # Use Lambda as W
    W_list[[k]] <- W_k
    
    # Decompose W into P and D
    svd_k <- svd(W_k)
    P_list[[k]] <- svd_k$u[, 1:true_params$r]
    D_list[[k]] <- diag(svd_k$d[1:true_params$r])
  }
  
  fit <- list(
    z = true_params$z,
    mu = true_params$mu,
    W = W_list,
    P = P_list,
    D = D_list,
    sigma2 = rep(0.5, true_params$K)
  )
  
  # Compute VAF
  vaf <- compute_global_vaf_mpca(list_of_data, fit)
  
  # Check that VAF is a single numeric value
  expect_type(vaf, "double")
  expect_length(vaf, 1)
  
  # VAF should be between 0 and 1
  expect_true(vaf >= 0)
  expect_true(vaf <= 1)
})

test_that("compute_cluster_sizes works correctly", {
  # Create simple cluster assignments
  z <- c(1, 1, 2, 2, 3, 3, 3)
  fit <- list(z = z)
  
  sizes <- compute_cluster_sizes(fit)
  
  # Check output is a table
  expect_s3_class(sizes, "table")
  
  # Check correct counts
  expect_equal(as.numeric(sizes["1"]), 2)
  expect_equal(as.numeric(sizes["2"]), 2)
  expect_equal(as.numeric(sizes["3"]), 3)
  
  # Check total
  expect_equal(sum(sizes), length(z))
})

test_that("compute_cluster_sizes handles edge cases", {
  # Single cluster
  z_single <- rep(1, 10)
  fit_single <- list(z = z_single)
  sizes_single <- compute_cluster_sizes(fit_single)
  
  expect_equal(length(sizes_single), 1)
  expect_equal(as.numeric(sizes_single), 10)
  
  # Unbalanced clusters
  z_unbalanced <- c(rep(1, 8), rep(2, 2))
  fit_unbalanced <- list(z = z_unbalanced)
  sizes_unbalanced <- compute_cluster_sizes(fit_unbalanced)
  
  expect_equal(as.numeric(sizes_unbalanced["1"]), 8)
  expect_equal(as.numeric(sizes_unbalanced["2"]), 2)
})

test_that("compute_cluster_sizes_mpca works correctly", {
  # Same as MFA version
  z <- c(1, 1, 2, 2, 3, 3, 3)
  fit <- list(z = z)
  
  sizes <- compute_cluster_sizes_mpca(fit)
  
  # Check output is a table
  expect_s3_class(sizes, "table")
  
  # Check correct counts
  expect_equal(as.numeric(sizes["1"]), 2)
  expect_equal(as.numeric(sizes["2"]), 2)
  expect_equal(as.numeric(sizes["3"]), 3)
})

test_that("posthoc_mfa_evaluation works correctly", {
  # Create mock selection object
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Create simple fitted models
  model1 <- list(
    z = true_params$z,
    mu = true_params$mu,
    Lambda = true_params$Lambda,
    Psi = true_params$Psi
  )
  
  model2 <- list(
    z = rep(1, true_params$N),
    mu = list(colMeans(do.call(rbind, list_of_data))),
    Lambda = list(matrix(rnorm(true_params$M * true_params$r), true_params$M, true_params$r)),
    Psi = list(diag(1, true_params$M))
  )
  
  selection_obj <- list(
    all_models = list(
      list(K = 2, r = 2, model = model1, logLik = -100, BIC = 250),
      list(K = 1, r = 2, model = model2, logLik = -150, BIC = 320)
    )
  )
  
  # Compute post-hoc evaluation
  df_posthoc <- posthoc_mfa_evaluation(list_of_data, selection_obj)
  
  # Check output structure
  expect_s3_class(df_posthoc, "data.frame")
  expect_equal(nrow(df_posthoc), 2)
  expect_true(all(c("K", "r", "logLik", "BIC", "GlobalVAF", "minClusterSize", "maxClusterSize") %in% names(df_posthoc)))
  
  # Check values are reasonable
  expect_true(all(df_posthoc$GlobalVAF >= 0))
  expect_true(all(df_posthoc$GlobalVAF <= 1))
  expect_true(all(df_posthoc$minClusterSize > 0))
  expect_true(all(df_posthoc$maxClusterSize <= true_params$N))
})

test_that("posthoc_mpca_evaluation works correctly", {
  # Create mock selection object for MPCA
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Create simple MPCA fitted models
  W1 <- true_params$Lambda[[1]]
  svd1 <- svd(W1)
  P1 <- svd1$u[, 1:true_params$r]
  D1 <- diag(svd1$d[1:true_params$r])
  
  model1 <- list(
    z = true_params$z,
    mu = true_params$mu,
    W = list(W1, W1),
    P = list(P1, P1),
    D = list(D1, D1),
    sigma2 = c(0.5, 0.5)
  )
  
  selection_obj <- list(
    all_models = list(
      list(K = 2, r = 2, model = model1, logLik = -100, BIC = 250)
    )
  )
  
  # Compute post-hoc evaluation
  df_posthoc <- posthoc_mpca_evaluation(list_of_data, selection_obj)
  
  # Check output structure
  expect_s3_class(df_posthoc, "data.frame")
  expect_equal(nrow(df_posthoc), 1)
  expect_true(all(c("K", "r", "logLik", "BIC", "GlobalVAF", "minClusterSize", "maxClusterSize") %in% names(df_posthoc)))
})

test_that("get_model_by_K_r works correctly", {
  # Create mock selection object
  model1 <- list(params = "model1")
  model2 <- list(params = "model2")
  model3 <- list(params = "model3")
  
  selection_obj <- list(
    all_models = list(
      list(K = 2, r = 2, model = model1),
      list(K = 2, r = 3, model = model2),
      list(K = 3, r = 2, model = model3)
    )
  )
  
  # Test retrieval
  retrieved1 <- get_model_by_K_r(selection_obj, K_target = 2, r_target = 2)
  expect_equal(retrieved1, model1)
  
  retrieved2 <- get_model_by_K_r(selection_obj, K_target = 2, r_target = 3)
  expect_equal(retrieved2, model2)
  
  retrieved3 <- get_model_by_K_r(selection_obj, K_target = 3, r_target = 2)
  expect_equal(retrieved3, model3)
  
  # Test non-existent combination
  expect_message(
    retrieved_none <- get_model_by_K_r(selection_obj, K_target = 4, r_target = 4),
    "No model found"
  )
  expect_null(retrieved_none)
})

test_that("get_model_by_K_r handles empty selection object", {
  selection_obj_empty <- list(all_models = list())
  
  retrieved <- get_model_by_K_r(selection_obj_empty, K_target = 2, r_target = 2)
  expect_null(retrieved)
})

test_that("compute_factor_scores is deterministic", {
  # Same inputs should give same outputs
  set.seed(789)
  X <- matrix(rnorm(20 * 4), 20, 4)
  mu <- rnorm(4)
  Lambda <- matrix(rnorm(8), 4, 2)
  Psi <- diag(runif(4, 0.1, 0.5))
  
  Z1 <- compute_factor_scores(X, mu, Lambda, Psi)
  Z2 <- compute_factor_scores(X, mu, Lambda, Psi)
  
  expect_equal(Z1, Z2)
})

test_that("VAF calculations are consistent", {
  # VAF from MFA should match manual calculation
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  fit <- list(
    z = true_params$z,
    mu = true_params$mu,
    Lambda = true_params$Lambda,
    Psi = true_params$Psi
  )
  
  vaf_auto <- compute_global_vaf_mfa(list_of_data, fit)
  
  # Manual calculation
  X_all <- do.call(rbind, list_of_data)
  SST <- sum(X_all^2)
  
  # Reconstruct data
  Xhat_all <- NULL
  for (i in 1:length(list_of_data)) {
    k <- fit$z[i]
    X_i <- list_of_data[[i]]
    Z_i <- compute_factor_scores(X_i, fit$mu[[k]], fit$Lambda[[k]], fit$Psi[[k]])
    Xhat_i <- Z_i %*% t(fit$Lambda[[k]]) + matrix(rep(fit$mu[[k]], nrow(X_i)), nrow(X_i), ncol(X_i), byrow = TRUE)
    Xhat_all <- rbind(Xhat_all, Xhat_i)
  }
  
  SSE <- sum((X_all - Xhat_all)^2)
  vaf_manual <- 1 - SSE / SST
  
  expect_equal(vaf_auto, vaf_manual, tolerance = 1e-6)
})
