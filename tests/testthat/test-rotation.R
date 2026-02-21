# Unit Tests for Rotation Functions
# Tests for rotation.R functions

# Load test fixtures
small_data <- readRDS("fixtures/small_test_data.rds")
medium_data <- readRDS("fixtures/medium_test_data.rds")

test_that("rotate_mfa_model preserves factor subspace dimensions", {
  # Create a mock MFA fit
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  mfa_fit <- list(
    z = true_params$z,
    mu = true_params$mu,
    Lambda = true_params$Lambda,
    Psi = true_params$Psi
  )
  
  # Rotate the model
  rotated <- rotate_mfa_model(mfa_fit, rotation = "varimax", rotate_scores = FALSE)
  
  # Check that dimensions are preserved
  expect_equal(length(rotated$Lambda), length(mfa_fit$Lambda))
  for (k in 1:true_params$K) {
    expect_equal(dim(rotated$Lambda[[k]]), dim(mfa_fit$Lambda[[k]]))
  }
  
  # Check that Lambda is still a valid matrix
  for (k in 1:true_params$K) {
    expect_true(all(is.finite(rotated$Lambda[[k]])))
  }
})

test_that("rotate_mfa_model preserves factor subspace (orthogonality)", {
  # Create a mock MFA fit with well-conditioned loadings
  set.seed(123)
  M <- 6
  r <- 3
  K <- 2
  
  Lambda_list <- list()
  for (k in 1:K) {
    # Create orthogonal loadings
    Q <- qr.Q(qr(matrix(rnorm(M * r), M, r)))
    Lambda_list[[k]] <- Q
  }
  
  mfa_fit <- list(
    z = c(1, 1, 2, 2),
    mu = list(rnorm(M), rnorm(M)),
    Lambda = Lambda_list,
    Psi = list(diag(0.5, M), diag(0.5, M))
  )
  
  rotated <- rotate_mfa_model(mfa_fit, rotation = "varimax", rotate_scores = FALSE)
  
  # Check that the column space is preserved (up to rotation)
  # The Gram matrix Lambda^T Lambda should have similar eigenvalues
  for (k in 1:K) {
    gram_orig <- t(mfa_fit$Lambda[[k]]) %*% mfa_fit$Lambda[[k]]
    gram_rot <- t(rotated$Lambda[[k]]) %*% rotated$Lambda[[k]]
    
    eig_orig <- sort(eigen(gram_orig)$values, decreasing = TRUE)
    eig_rot <- sort(eigen(gram_rot)$values, decreasing = TRUE)
    
    expect_equal(eig_orig, eig_rot, tolerance = 1e-6)
  }
})

test_that("rotate_mfa_model with rotate_scores maintains consistency", {
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
  
  # Rotate with scores
  rotated <- rotate_mfa_model(mfa_fit, rotation = "varimax", rotate_scores = TRUE)
  
  # Check that factor scores are rotated
  expect_true(!is.null(rotated$factor_scores))
  expect_equal(length(rotated$factor_scores), length(mfa_fit$factor_scores))
  
  # Check that reconstruction is preserved
  # For each subject, X_reconstructed should be the same before and after rotation
  for (i in 1:true_params$N) {
    k <- mfa_fit$z[i]
    
    # Original reconstruction
    X_orig <- mfa_fit$factor_scores[[i]] %*% t(mfa_fit$Lambda[[k]])
    X_orig <- sweep(X_orig, 2, mfa_fit$mu[[k]], FUN = "+")
    
    # Rotated reconstruction
    X_rot <- rotated$factor_scores[[i]] %*% t(rotated$Lambda[[k]])
    X_rot <- sweep(X_rot, 2, rotated$mu[[k]], FUN = "+")
    
    # Should be identical (up to numerical precision)
    expect_equal(X_orig, X_rot, tolerance = 1e-10)
  }
})

test_that("rotate_mfa_model handles single factor gracefully", {
  # Create a mock MFA fit with r=1
  set.seed(456)
  M <- 4
  r <- 1
  K <- 2
  
  mfa_fit <- list(
    z = c(1, 1, 2),
    mu = list(rnorm(M), rnorm(M)),
    Lambda = list(matrix(rnorm(M * r), M, r), matrix(rnorm(M * r), M, r)),
    Psi = list(diag(0.5, M), diag(0.5, M))
  )
  
  # Rotation should skip single-factor clusters
  rotated <- rotate_mfa_model(mfa_fit, rotation = "varimax", rotate_scores = FALSE)
  
  # Lambda should be unchanged for single-factor case
  expect_equal(rotated$Lambda[[1]], mfa_fit$Lambda[[1]])
  expect_equal(rotated$Lambda[[2]], mfa_fit$Lambda[[2]])
})

test_that("rotate_mfa_model errors when factor_scores missing but rotate_scores=TRUE", {
  # Create a mock MFA fit without factor scores
  mfa_fit <- list(
    z = c(1, 1, 2),
    mu = list(rnorm(4), rnorm(4)),
    Lambda = list(matrix(rnorm(8), 4, 2), matrix(rnorm(8), 4, 2)),
    Psi = list(diag(0.5, 4), diag(0.5, 4))
  )
  
  # Should error when trying to rotate scores that don't exist
  expect_error(
    rotate_mfa_model(mfa_fit, rotation = "varimax", rotate_scores = TRUE),
    "factor_scores not found"
  )
})

test_that("rotate_mfa_model with cluster_ids rotates only specified clusters", {
  # Create a mock MFA fit
  set.seed(789)
  M <- 4
  r <- 2
  K <- 3
  
  Lambda_list <- list()
  for (k in 1:K) {
    Lambda_list[[k]] <- matrix(rnorm(M * r), M, r)
  }
  
  mfa_fit <- list(
    z = c(1, 1, 2, 2, 3, 3),
    mu = list(rnorm(M), rnorm(M), rnorm(M)),
    Lambda = Lambda_list,
    Psi = list(diag(0.5, M), diag(0.5, M), diag(0.5, M))
  )
  
  # Rotate only cluster 2
  rotated <- rotate_mfa_model(mfa_fit, rotation = "varimax", 
                              rotate_scores = FALSE, cluster_ids = 2)
  
  # Clusters 1 and 3 should be unchanged
  expect_equal(rotated$Lambda[[1]], mfa_fit$Lambda[[1]])
  expect_equal(rotated$Lambda[[3]], mfa_fit$Lambda[[3]])
  
  # Cluster 2 should be different (rotated)
  expect_false(isTRUE(all.equal(rotated$Lambda[[2]], mfa_fit$Lambda[[2]])))
})

test_that("rotate_mpca_model preserves factor subspace dimensions", {
  # Create a mock MPCA fit
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Create W matrices from Lambda
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
  
  # Rotate the model
  rotated <- rotate_mpca_model(mpca_fit, rotation = "varimax", rotate_scores = FALSE)
  
  # Check that dimensions are preserved
  expect_equal(length(rotated$W), length(mpca_fit$W))
  for (k in 1:true_params$K) {
    expect_equal(dim(rotated$W[[k]]), dim(mpca_fit$W[[k]]))
  }
  
  # Check that W is still a valid matrix
  for (k in 1:true_params$K) {
    expect_true(all(is.finite(rotated$W[[k]])))
  }
})

test_that("rotate_mpca_model with rotate_scores maintains consistency", {
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
    
    # Compute PPCA scores
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
  
  # Rotate with scores
  rotated <- rotate_mpca_model(mpca_fit, rotation = "varimax", rotate_scores = TRUE)
  
  # Check that factor scores are rotated
  expect_true(!is.null(rotated$factor_scores))
  expect_equal(length(rotated$factor_scores), length(mpca_fit$factor_scores))
  
  # Check that reconstruction is preserved
  for (i in 1:true_params$N) {
    k <- mpca_fit$z[i]
    
    # Original reconstruction
    X_orig <- mpca_fit$factor_scores[[i]] %*% t(mpca_fit$W[[k]])
    X_orig <- sweep(X_orig, 2, mpca_fit$mu[[k]], FUN = "+")
    
    # Rotated reconstruction
    X_rot <- rotated$factor_scores[[i]] %*% t(rotated$W[[k]])
    X_rot <- sweep(X_rot, 2, rotated$mu[[k]], FUN = "+")
    
    # Should be identical (up to numerical precision)
    expect_equal(X_orig, X_rot, tolerance = 1e-10)
  }
})

test_that("rotate_mpca_model handles single factor with warning", {
  # Create a mock MPCA fit with r=1
  set.seed(456)
  M <- 4
  r <- 1
  K <- 2
  
  mpca_fit <- list(
    z = c(1, 1, 2),
    mu = list(rnorm(M), rnorm(M)),
    W = list(matrix(rnorm(M * r), M, r), matrix(rnorm(M * r), M, r)),
    sigma2 = c(0.5, 0.5)
  )
  
  # Should warn about single factor
  expect_warning(
    rotated <- rotate_mpca_model(mpca_fit, rotation = "varimax", rotate_scores = FALSE),
    "only 1 factor"
  )
  
  # W should be unchanged for single-factor case
  expect_equal(rotated$W[[1]], mpca_fit$W[[1]])
  expect_equal(rotated$W[[2]], mpca_fit$W[[2]])
})

test_that("rotate_mpca_model errors when factor_scores missing but rotate_scores=TRUE", {
  # Create a mock MPCA fit without factor scores
  mpca_fit <- list(
    z = c(1, 1, 2),
    mu = list(rnorm(4), rnorm(4)),
    W = list(matrix(rnorm(8), 4, 2), matrix(rnorm(8), 4, 2)),
    sigma2 = c(0.5, 0.5)
  )
  
  # Should error when trying to rotate scores that don't exist
  expect_error(
    rotate_mpca_model(mpca_fit, rotation = "varimax", rotate_scores = TRUE),
    "factor_scores not found"
  )
})

test_that("rotate_mpca_model with cluster_ids rotates only specified clusters", {
  # Create a mock MPCA fit
  set.seed(789)
  M <- 4
  r <- 2
  K <- 3
  
  W_list <- list()
  for (k in 1:K) {
    W_list[[k]] <- matrix(rnorm(M * r), M, r)
  }
  
  mpca_fit <- list(
    z = c(1, 1, 2, 2, 3, 3),
    mu = list(rnorm(M), rnorm(M), rnorm(M)),
    W = W_list,
    sigma2 = rep(0.5, K)
  )
  
  # Rotate only cluster 2
  rotated <- rotate_mpca_model(mpca_fit, rotation = "varimax", 
                               rotate_scores = FALSE, cluster_ids = 2)
  
  # Clusters 1 and 3 should be unchanged
  expect_equal(rotated$W[[1]], mpca_fit$W[[1]])
  expect_equal(rotated$W[[3]], mpca_fit$W[[3]])
  
  # Cluster 2 should be different (rotated)
  expect_false(isTRUE(all.equal(rotated$W[[2]], mpca_fit$W[[2]])))
})

test_that("varimax rotation is reproducible", {
  # Create a mock MFA fit
  set.seed(999)
  M <- 6
  r <- 3
  
  Lambda <- matrix(rnorm(M * r), M, r)
  
  mfa_fit <- list(
    z = c(1, 1),
    mu = list(rnorm(M)),
    Lambda = list(Lambda),
    Psi = list(diag(0.5, M))
  )
  
  # Rotate twice
  rotated1 <- rotate_mfa_model(mfa_fit, rotation = "varimax", rotate_scores = FALSE)
  rotated2 <- rotate_mfa_model(mfa_fit, rotation = "varimax", rotate_scores = FALSE)
  
  # Should be identical
  expect_equal(rotated1$Lambda[[1]], rotated2$Lambda[[1]])
})

test_that("rotation preserves total variance explained", {
  # Create a mock MFA fit
  set.seed(111)
  M <- 6
  r <- 3
  
  Lambda <- matrix(rnorm(M * r), M, r)
  
  mfa_fit <- list(
    z = c(1, 1),
    mu = list(rnorm(M)),
    Lambda = list(Lambda),
    Psi = list(diag(0.5, M))
  )
  
  # Rotate
  rotated <- rotate_mfa_model(mfa_fit, rotation = "varimax", rotate_scores = FALSE)
  
  # Total variance (sum of squared loadings) should be preserved
  var_orig <- sum(mfa_fit$Lambda[[1]]^2)
  var_rot <- sum(rotated$Lambda[[1]]^2)
  
  expect_equal(var_orig, var_rot, tolerance = 1e-10)
})
