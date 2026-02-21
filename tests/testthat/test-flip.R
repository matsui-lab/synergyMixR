# Unit Tests for Sign Alignment Functions
# Tests for flip.R functions

# Load test fixtures
small_data <- readRDS("fixtures/small_test_data.rds")
medium_data <- readRDS("fixtures/medium_test_data.rds")

test_that("align_factor_signs_mfa flips signs correctly", {
  # Create a mock MFA fit with negative max loadings
  set.seed(123)
  M <- 6
  r <- 3
  K <- 2
  
  Lambda_list <- list()
  for (k in 1:K) {
    Lambda_k <- matrix(rnorm(M * r), M, r)
    # Force the max absolute loading in column 1 to be negative
    max_idx <- which.max(abs(Lambda_k[, 1]))
    Lambda_k[max_idx, 1] <- -abs(Lambda_k[max_idx, 1])
    Lambda_list[[k]] <- Lambda_k
  }
  
  mfa_fit <- list(
    z = c(1, 1, 2, 2),
    mu = list(rnorm(M), rnorm(M)),
    Lambda = Lambda_list,
    Psi = list(diag(0.5, M), diag(0.5, M))
  )
  
  # Align signs
  aligned <- align_factor_signs_mfa(mfa_fit, method = "max_abs")
  
  # Check that the max absolute loading is now positive in column 1
  for (k in 1:K) {
    max_idx <- which.max(abs(aligned$Lambda[[k]][, 1]))
    expect_true(aligned$Lambda[[k]][max_idx, 1] > 0)
  }
})

test_that("align_factor_signs_mfa preserves magnitude", {
  # Create a mock MFA fit
  set.seed(456)
  M <- 6
  r <- 3
  K <- 2
  
  Lambda_list <- list()
  for (k in 1:K) {
    Lambda_list[[k]] <- matrix(rnorm(M * r), M, r)
  }
  
  mfa_fit <- list(
    z = c(1, 1, 2, 2),
    mu = list(rnorm(M), rnorm(M)),
    Lambda = Lambda_list,
    Psi = list(diag(0.5, M), diag(0.5, M))
  )
  
  # Align signs
  aligned <- align_factor_signs_mfa(mfa_fit, method = "max_abs")
  
  # Check that magnitudes are preserved (only signs change)
  for (k in 1:K) {
    for (j in 1:r) {
      orig_abs <- abs(mfa_fit$Lambda[[k]][, j])
      aligned_abs <- abs(aligned$Lambda[[k]][, j])
      expect_equal(orig_abs, aligned_abs, tolerance = 1e-10)
    }
  }
})

test_that("align_factor_signs_mfa with factor_scores maintains consistency", {
  # Create a mock MFA fit with factor scores
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  mfa_fit <- list(
    z = true_params$z,
    mu = true_params$mu,
    Lambda = true_params$Lambda,
    Psi = true_params$Psi
  )
  
  # Compute factor scores
  mfa_fit <- compute_factor_scores_mfa(mfa_fit, list_of_data)
  
  # Align signs
  aligned <- align_factor_signs_mfa(mfa_fit, method = "max_abs")
  
  # Check that reconstruction is preserved
  for (i in 1:true_params$N) {
    k <- mfa_fit$z[i]
    
    # Original reconstruction
    X_orig <- mfa_fit$factor_scores[[i]] %*% t(mfa_fit$Lambda[[k]])
    X_orig <- sweep(X_orig, 2, mfa_fit$mu[[k]], FUN = "+")
    
    # Aligned reconstruction
    X_aligned <- aligned$factor_scores[[i]] %*% t(aligned$Lambda[[k]])
    X_aligned <- sweep(X_aligned, 2, aligned$mu[[k]], FUN = "+")
    
    # Should be identical (up to numerical precision)
    expect_equal(X_orig, X_aligned, tolerance = 1e-10)
  }
})

test_that("align_factor_signs_mfa flips factor_scores consistently with Lambda", {
  # Create a mock MFA fit with factor scores
  set.seed(789)
  M <- 4
  r <- 2
  K <- 2
  N <- 4
  T_i <- 10
  
  Lambda_list <- list()
  for (k in 1:K) {
    Lambda_k <- matrix(rnorm(M * r), M, r)
    # Force negative max loading in first column
    max_idx <- which.max(abs(Lambda_k[, 1]))
    Lambda_k[max_idx, 1] <- -abs(Lambda_k[max_idx, 1])
    Lambda_list[[k]] <- Lambda_k
  }
  
  z <- c(1, 1, 2, 2)
  factor_scores_list <- list()
  for (i in 1:N) {
    factor_scores_list[[i]] <- matrix(rnorm(T_i * r), T_i, r)
  }
  
  mfa_fit <- list(
    z = z,
    mu = list(rnorm(M), rnorm(M)),
    Lambda = Lambda_list,
    Psi = list(diag(0.5, M), diag(0.5, M)),
    factor_scores = factor_scores_list
  )
  
  # Store original values
  orig_Lambda_1_col1 <- mfa_fit$Lambda[[1]][, 1]
  orig_scores_1_col1 <- mfa_fit$factor_scores[[1]][, 1]
  
  # Align signs
  aligned <- align_factor_signs_mfa(mfa_fit, method = "max_abs")
  
  # If Lambda column was flipped, factor scores should also be flipped
  max_idx <- which.max(abs(orig_Lambda_1_col1))
  if (orig_Lambda_1_col1[max_idx] < 0) {
    # Lambda should be flipped
    expect_equal(aligned$Lambda[[1]][, 1], -orig_Lambda_1_col1, tolerance = 1e-10)
    # Factor scores for subjects in cluster 1 should also be flipped
    expect_equal(aligned$factor_scores[[1]][, 1], -orig_scores_1_col1, tolerance = 1e-10)
  }
})

test_that("align_factor_signs_mfa without factor_scores works", {
  # Create a mock MFA fit without factor scores
  set.seed(111)
  M <- 4
  r <- 2
  K <- 2
  
  Lambda_list <- list()
  for (k in 1:K) {
    Lambda_list[[k]] <- matrix(rnorm(M * r), M, r)
  }
  
  mfa_fit <- list(
    z = c(1, 1, 2, 2),
    mu = list(rnorm(M), rnorm(M)),
    Lambda = Lambda_list,
    Psi = list(diag(0.5, M), diag(0.5, M))
  )
  
  # Should work without error
  expect_no_error(
    aligned <- align_factor_signs_mfa(mfa_fit, method = "max_abs")
  )
  
  # Check that Lambda is modified
  expect_equal(length(aligned$Lambda), K)
})

test_that("align_factor_signs_mfa with cluster_ids aligns only specified clusters", {
  # Create a mock MFA fit
  set.seed(222)
  M <- 4
  r <- 2
  K <- 3
  
  Lambda_list <- list()
  for (k in 1:K) {
    Lambda_k <- matrix(rnorm(M * r), M, r)
    # Force negative max loading in first column
    max_idx <- which.max(abs(Lambda_k[, 1]))
    Lambda_k[max_idx, 1] <- -abs(Lambda_k[max_idx, 1])
    Lambda_list[[k]] <- Lambda_k
  }
  
  mfa_fit <- list(
    z = c(1, 1, 2, 2, 3, 3),
    mu = list(rnorm(M), rnorm(M), rnorm(M)),
    Lambda = Lambda_list,
    Psi = list(diag(0.5, M), diag(0.5, M), diag(0.5, M))
  )
  
  # Store original cluster 1 and 3
  orig_Lambda_1 <- mfa_fit$Lambda[[1]]
  orig_Lambda_3 <- mfa_fit$Lambda[[3]]
  
  # Align only cluster 2
  aligned <- align_factor_signs_mfa(mfa_fit, method = "max_abs", cluster_ids = 2)
  
  # Clusters 1 and 3 should be unchanged
  expect_equal(aligned$Lambda[[1]], orig_Lambda_1)
  expect_equal(aligned$Lambda[[3]], orig_Lambda_3)
})

test_that("align_factor_signs_mpca flips signs correctly", {
  # Create a mock MPCA fit with negative max loadings
  set.seed(333)
  M <- 6
  r <- 3
  K <- 2
  
  W_list <- list()
  for (k in 1:K) {
    W_k <- matrix(rnorm(M * r), M, r)
    # Force the max absolute loading in column 1 to be negative
    max_idx <- which.max(abs(W_k[, 1]))
    W_k[max_idx, 1] <- -abs(W_k[max_idx, 1])
    W_list[[k]] <- W_k
  }
  
  mpca_fit <- list(
    z = c(1, 1, 2, 2),
    mu = list(rnorm(M), rnorm(M)),
    W = W_list,
    sigma2 = c(0.5, 0.5)
  )
  
  # Align signs
  aligned <- align_factor_signs_mpca(mpca_fit, method = "max_abs")
  
  # Check that the max absolute loading is now positive in column 1
  for (k in 1:K) {
    max_idx <- which.max(abs(aligned$W[[k]][, 1]))
    expect_true(aligned$W[[k]][max_idx, 1] > 0)
  }
})

test_that("align_factor_signs_mpca preserves magnitude", {
  # Create a mock MPCA fit
  set.seed(444)
  M <- 6
  r <- 3
  K <- 2
  
  W_list <- list()
  for (k in 1:K) {
    W_list[[k]] <- matrix(rnorm(M * r), M, r)
  }
  
  mpca_fit <- list(
    z = c(1, 1, 2, 2),
    mu = list(rnorm(M), rnorm(M)),
    W = W_list,
    sigma2 = c(0.5, 0.5)
  )
  
  # Align signs
  aligned <- align_factor_signs_mpca(mpca_fit, method = "max_abs")
  
  # Check that magnitudes are preserved (only signs change)
  for (k in 1:K) {
    for (j in 1:r) {
      orig_abs <- abs(mpca_fit$W[[k]][, j])
      aligned_abs <- abs(aligned$W[[k]][, j])
      expect_equal(orig_abs, aligned_abs, tolerance = 1e-10)
    }
  }
})

test_that("align_factor_signs_mpca with factor_scores maintains consistency", {
  # Create a mock MPCA fit with factor scores
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Create W matrices
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
  
  # Compute factor scores
  factor_scores_list <- list()
  for (i in 1:true_params$N) {
    k <- mpca_fit$z[i]
    X_i <- list_of_data[[i]]
    
    T_i <- nrow(X_i)
    r <- ncol(W_list[[k]])
    Xc <- sweep(X_i, 2, mpca_fit$mu[[k]], FUN = "-")
    
    A <- t(W_list[[k]]) %*% W_list[[k]] + mpca_fit$sigma2[k] * diag(r)
    A_inv <- solve(A)
    W_t <- t(W_list[[k]])
    
    Z_i <- matrix(0, nrow = T_i, ncol = r)
    for (t in 1:T_i) {
      Z_i[t, ] <- A_inv %*% W_t %*% Xc[t, ]
    }
    factor_scores_list[[i]] <- Z_i
  }
  mpca_fit$factor_scores <- factor_scores_list
  
  # Align signs
  aligned <- align_factor_signs_mpca(mpca_fit, method = "max_abs")
  
  # Check that reconstruction is preserved
  for (i in 1:true_params$N) {
    k <- mpca_fit$z[i]
    
    # Original reconstruction
    X_orig <- mpca_fit$factor_scores[[i]] %*% t(mpca_fit$W[[k]])
    X_orig <- sweep(X_orig, 2, mpca_fit$mu[[k]], FUN = "+")
    
    # Aligned reconstruction
    X_aligned <- aligned$factor_scores[[i]] %*% t(aligned$W[[k]])
    X_aligned <- sweep(X_aligned, 2, aligned$mu[[k]], FUN = "+")
    
    # Should be identical (up to numerical precision)
    expect_equal(X_orig, X_aligned, tolerance = 1e-10)
  }
})

test_that("align_factor_signs_mpca flips factor_scores consistently with W", {
  # Create a mock MPCA fit with factor scores
  set.seed(555)
  M <- 4
  r <- 2
  K <- 2
  N <- 4
  T_i <- 10
  
  W_list <- list()
  for (k in 1:K) {
    W_k <- matrix(rnorm(M * r), M, r)
    # Force negative max loading in first column
    max_idx <- which.max(abs(W_k[, 1]))
    W_k[max_idx, 1] <- -abs(W_k[max_idx, 1])
    W_list[[k]] <- W_k
  }
  
  z <- c(1, 1, 2, 2)
  factor_scores_list <- list()
  for (i in 1:N) {
    factor_scores_list[[i]] <- matrix(rnorm(T_i * r), T_i, r)
  }
  
  mpca_fit <- list(
    z = z,
    mu = list(rnorm(M), rnorm(M)),
    W = W_list,
    sigma2 = c(0.5, 0.5),
    factor_scores = factor_scores_list
  )
  
  # Store original values
  orig_W_1_col1 <- mpca_fit$W[[1]][, 1]
  orig_scores_1_col1 <- mpca_fit$factor_scores[[1]][, 1]
  
  # Align signs
  aligned <- align_factor_signs_mpca(mpca_fit, method = "max_abs")
  
  # If W column was flipped, factor scores should also be flipped
  max_idx <- which.max(abs(orig_W_1_col1))
  if (orig_W_1_col1[max_idx] < 0) {
    # W should be flipped
    expect_equal(aligned$W[[1]][, 1], -orig_W_1_col1, tolerance = 1e-10)
    # Factor scores for subjects in cluster 1 should also be flipped
    expect_equal(aligned$factor_scores[[1]][, 1], -orig_scores_1_col1, tolerance = 1e-10)
  }
})

test_that("align_factor_signs_mpca without factor_scores works", {
  # Create a mock MPCA fit without factor scores
  set.seed(666)
  M <- 4
  r <- 2
  K <- 2
  
  W_list <- list()
  for (k in 1:K) {
    W_list[[k]] <- matrix(rnorm(M * r), M, r)
  }
  
  mpca_fit <- list(
    z = c(1, 1, 2, 2),
    mu = list(rnorm(M), rnorm(M)),
    W = W_list,
    sigma2 = c(0.5, 0.5)
  )
  
  # Should work without error
  expect_no_error(
    aligned <- align_factor_signs_mpca(mpca_fit, method = "max_abs")
  )
  
  # Check that W is modified
  expect_equal(length(aligned$W), K)
})

test_that("align_factor_signs_mpca with cluster_ids aligns only specified clusters", {
  # Create a mock MPCA fit
  set.seed(777)
  M <- 4
  r <- 2
  K <- 3
  
  W_list <- list()
  for (k in 1:K) {
    W_k <- matrix(rnorm(M * r), M, r)
    # Force negative max loading in first column
    max_idx <- which.max(abs(W_k[, 1]))
    W_k[max_idx, 1] <- -abs(W_k[max_idx, 1])
    W_list[[k]] <- W_k
  }
  
  mpca_fit <- list(
    z = c(1, 1, 2, 2, 3, 3),
    mu = list(rnorm(M), rnorm(M), rnorm(M)),
    W = W_list,
    sigma2 = c(0.5, 0.5, 0.5)
  )
  
  # Store original cluster 1 and 3
  orig_W_1 <- mpca_fit$W[[1]]
  orig_W_3 <- mpca_fit$W[[3]]
  
  # Align only cluster 2
  aligned <- align_factor_signs_mpca(mpca_fit, method = "max_abs", cluster_ids = 2)
  
  # Clusters 1 and 3 should be unchanged
  expect_equal(aligned$W[[1]], orig_W_1)
  expect_equal(aligned$W[[3]], orig_W_3)
})

test_that("sign alignment is idempotent", {
  # Create a mock MFA fit
  set.seed(888)
  M <- 4
  r <- 2
  K <- 2
  
  Lambda_list <- list()
  for (k in 1:K) {
    Lambda_list[[k]] <- matrix(rnorm(M * r), M, r)
  }
  
  mfa_fit <- list(
    z = c(1, 1, 2, 2),
    mu = list(rnorm(M), rnorm(M)),
    Lambda = Lambda_list,
    Psi = list(diag(0.5, M), diag(0.5, M))
  )
  
  # Align once
  aligned1 <- align_factor_signs_mfa(mfa_fit, method = "max_abs")
  
  # Align again
  aligned2 <- align_factor_signs_mfa(aligned1, method = "max_abs")
  
  # Should be identical (idempotent)
  expect_equal(aligned1$Lambda, aligned2$Lambda)
})

test_that("sign alignment handles all positive loadings", {
  # Create a mock MFA fit with all positive max loadings
  set.seed(999)
  M <- 4
  r <- 2
  K <- 2
  
  Lambda_list <- list()
  for (k in 1:K) {
    Lambda_k <- matrix(abs(rnorm(M * r)), M, r)
    Lambda_list[[k]] <- Lambda_k
  }
  
  mfa_fit <- list(
    z = c(1, 1, 2, 2),
    mu = list(rnorm(M), rnorm(M)),
    Lambda = Lambda_list,
    Psi = list(diag(0.5, M), diag(0.5, M))
  )
  
  # Align signs (should not change anything)
  aligned <- align_factor_signs_mfa(mfa_fit, method = "max_abs")
  
  # Should be unchanged
  expect_equal(aligned$Lambda, mfa_fit$Lambda)
})
