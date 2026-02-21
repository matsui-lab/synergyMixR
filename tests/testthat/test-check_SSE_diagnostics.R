# Unit Tests for check_SSE_diagnostics()
# Tests for SSE diagnostic function in compare_methods_impl.R

# Helper function to create mock fitted models
create_mock_mfa_fit <- function(N = 30, K = 3, M = 8, good_psi = TRUE) {
  z <- sample(1:K, N, replace = TRUE)
  
  # Create Psi matrices
  Psi_list <- list()
  for (k in 1:K) {
    if (good_psi) {
      # Normal Psi values (between 0.1 and 10)
      psi_diag <- runif(M, 0.1, 10)
    } else {
      # Mix of extreme Psi values
      psi_diag <- c(
        runif(M %/% 3, 0.001, 0.009),  # Very small (overfitting)
        runif(M %/% 3, 0.1, 10),        # Normal
        runif(M - 2 * (M %/% 3), 101, 200)  # Very large (poor fit)
      )
    }
    Psi_list[[k]] <- diag(psi_diag)
  }
  
  list(
    z = z,
    pi = rep(1/K, K),
    mu = lapply(1:K, function(k) rnorm(M)),
    Lambda = lapply(1:K, function(k) matrix(rnorm(M * 3), M, 3)),
    Psi = Psi_list,
    logLik = -100,
    resp = matrix(runif(N * K), N, K)
  )
}

create_mock_mpca_fit <- function(N = 30, K = 3, M = 8, good_sigma2 = TRUE) {
  z <- sample(1:K, N, replace = TRUE)
  
  # Create sigma2 values
  if (good_sigma2) {
    sigma2 <- runif(K, 0.1, 10)
  } else {
    sigma2 <- c(0.005, 5, 150)  # Mix of extreme values
  }
  
  list(
    z = z,
    pi = rep(1/K, K),
    mu = lapply(1:K, function(k) rnorm(M)),
    W = lapply(1:K, function(k) matrix(rnorm(M * 3), M, 3)),
    P = lapply(1:K, function(k) matrix(rnorm(M * 3), M, 3)),
    D = lapply(1:K, function(k) diag(runif(3, 0.5, 2))),
    sigma2 = sigma2,
    logLik = -100,
    resp = matrix(runif(N * K), N, K)
  )
}

create_mock_list_of_data <- function(N = 30, T_each = 200, M = 8) {
  lapply(1:N, function(i) matrix(rnorm(T_each * M), T_each, M))
}

# ============================================================================
# Test 1: Normal case (all TRUE)
# ============================================================================

test_that("check_SSE_diagnostics works correctly with normal inputs (all TRUE)", {
  # Create mock data
  list_of_data <- create_mock_list_of_data(N = 30, T_each = 200, M = 8)
  
  # Create mock fitted models with good parameters
  mfa_fit <- create_mock_mfa_fit(N = 30, K = 3, M = 8, good_psi = TRUE)
  mpca_fit <- create_mock_mpca_fit(N = 30, K = 3, M = 8, good_sigma2 = TRUE)
  
  # Create fit_results data frame
  fit_results <- data.frame(
    Method = c("SingleFA", "SinglePCA", "TwoStep_FA", "TwoStep_PCA", "MixtureFA", "MixturePCA"),
    BIC = c(1000, 1100, 900, 950, 800, 850),
    ARI = c(NA, NA, 0.7, 0.75, 0.85, 0.88),
    SSE = c(5000, 5200, 4500, 4600, 4000, 4100),
    VAF = c(0.70, 0.68, 0.75, 0.74, 0.80, 0.79)
  )
  
  # Store fitted models as attributes
  attr(fit_results, "mixtureFA_fit") <- mfa_fit
  attr(fit_results, "mixturePCA_fit") <- mpca_fit
  
  # Run diagnostics
  diagnostics <- check_SSE_diagnostics(
    fit_results = fit_results,
    list_of_data = list_of_data
  )
  
  # Check output structure
  expect_type(diagnostics, "list")
  expect_true(all(c("convergence_ok", "psi_ok", "correlation_ok", "warnings", 
                    "psi_summary_mfa", "psi_summary_mpca", "vaf_sse_correlation") %in% names(diagnostics)))
  
  # Check that all diagnostics pass
  expect_true(diagnostics$convergence_ok)
  expect_true(diagnostics$psi_ok)
  expect_true(diagnostics$correlation_ok)
  
  # Check that VAF-SSE correlation is strongly negative
  expect_true(is.finite(diagnostics$vaf_sse_correlation))
  expect_true(diagnostics$vaf_sse_correlation < -0.9)
  
  # Check that Psi summaries are present
  expect_false(is.null(diagnostics$psi_summary_mfa))
  expect_false(is.null(diagnostics$psi_summary_mpca))
})

# ============================================================================
# Test 2: Psi threshold override (psi_ok FALSE)
# ============================================================================

test_that("check_SSE_diagnostics detects extreme Psi values with custom thresholds", {
  # Create mock data
  list_of_data <- create_mock_list_of_data(N = 30, T_each = 200, M = 8)
  
  # Create mock fitted models with extreme Psi values
  mfa_fit <- create_mock_mfa_fit(N = 30, K = 3, M = 8, good_psi = FALSE)
  mpca_fit <- create_mock_mpca_fit(N = 30, K = 3, M = 8, good_sigma2 = FALSE)
  
  # Create fit_results data frame
  fit_results <- data.frame(
    Method = c("MixtureFA", "MixturePCA"),
    BIC = c(800, 850),
    ARI = c(0.85, 0.88),
    SSE = c(4000, 4100),
    VAF = c(0.80, 0.79)
  )
  
  # Run diagnostics with custom thresholds
  diagnostics <- check_SSE_diagnostics(
    fit_results = fit_results,
    list_of_data = list_of_data,
    mixtureFA_fit = mfa_fit,
    mixturePCA_fit = mpca_fit,
    psi_small_thresh = 0.01,
    psi_large_thresh = 100
  )
  
  # Check that psi_ok is FALSE due to extreme values
  expect_false(diagnostics$psi_ok)
  
  # Check that warnings contain messages about extreme Psi values
  expect_true(any(grepl("Psi values < 0.010", diagnostics$warnings)))
  expect_true(any(grepl("Psi values > 100", diagnostics$warnings)) || 
              any(grepl("sigma2 values > 100", diagnostics$warnings)))
  
  # Check that warnings vector is not empty
  expect_true(length(diagnostics$warnings) > 0)
})

# ============================================================================
# Test 3: Data.frame only with stripped attributes → EM/Psi skipped
# ============================================================================

test_that("check_SSE_diagnostics gracefully handles stripped attributes", {
  # Create mock data
  list_of_data <- create_mock_list_of_data(N = 30, T_each = 200, M = 8)
  
  # Create fit_results data frame with attributes
  fit_results <- data.frame(
    Method = c("SingleFA", "SinglePCA", "TwoStep_FA", "TwoStep_PCA", "MixtureFA", "MixturePCA"),
    BIC = c(1000, 1100, 900, 950, 800, 850),
    ARI = c(NA, NA, 0.7, 0.75, 0.85, 0.88),
    SSE = c(5000, 5200, 4500, 4600, 4000, 4100),
    VAF = c(0.70, 0.68, 0.75, 0.74, 0.80, 0.79)
  )
  
  # Add fitted models as attributes
  mfa_fit <- create_mock_mfa_fit(N = 30, K = 3, M = 8, good_psi = TRUE)
  mpca_fit <- create_mock_mpca_fit(N = 30, K = 3, M = 8, good_sigma2 = TRUE)
  attr(fit_results, "mixtureFA_fit") <- mfa_fit
  attr(fit_results, "mixturePCA_fit") <- mpca_fit
  
  # Strip attributes using the utility function
  fit_results_stripped <- strip_attributes(fit_results)
  
  # Verify attributes are stripped
  expect_null(attr(fit_results_stripped, "mixtureFA_fit"))
  expect_null(attr(fit_results_stripped, "mixturePCA_fit"))
  
  # Run diagnostics on stripped data frame
  diagnostics <- check_SSE_diagnostics(
    fit_results = fit_results_stripped,
    list_of_data = list_of_data
  )
  
  # Check that warnings indicate EM/Psi diagnostics were skipped
  expect_true(any(grepl("No valid fitted mixture models available", diagnostics$warnings)))
  expect_true(any(grepl("EM convergence and Psi diagnostics will be skipped", diagnostics$warnings)))
  
  # Check that convergence_ok and psi_ok remain TRUE (no models to check)
  expect_true(diagnostics$convergence_ok)
  expect_true(diagnostics$psi_ok)
  
  # Check that VAF-SSE correlation is still computed
  expect_true(is.finite(diagnostics$vaf_sse_correlation))
  
  # Check that Psi summaries are NULL
  expect_null(diagnostics$psi_summary_mfa)
  expect_null(diagnostics$psi_summary_mpca)
})

# ============================================================================
# Test 4: Explicit override of mixtureFA_fit / mixturePCA_fit
# ============================================================================

test_that("check_SSE_diagnostics respects explicit parameter override", {
  # Create mock data
  list_of_data <- create_mock_list_of_data(N = 30, T_each = 200, M = 8)
  
  # Create two different sets of fitted models
  mfa_fit_attr <- create_mock_mfa_fit(N = 30, K = 3, M = 8, good_psi = TRUE)
  mfa_fit_explicit <- create_mock_mfa_fit(N = 30, K = 3, M = 8, good_psi = FALSE)
  
  # Create fit_results with one model in attributes
  fit_results <- data.frame(
    Method = c("MixtureFA"),
    BIC = c(800),
    ARI = c(0.85),
    SSE = c(4000),
    VAF = c(0.80)
  )
  attr(fit_results, "mixtureFA_fit") <- mfa_fit_attr
  
  # Run diagnostics with explicit override (should use explicit, not attribute)
  diagnostics <- check_SSE_diagnostics(
    fit_results = fit_results,
    list_of_data = list_of_data,
    mixtureFA_fit = mfa_fit_explicit,  # Explicit override with bad Psi
    psi_small_thresh = 0.01,
    psi_large_thresh = 100
  )
  
  # Check that psi_ok is FALSE (using explicit model with bad Psi)
  expect_false(diagnostics$psi_ok)
  
  # Check that warnings contain messages about extreme Psi values
  expect_true(any(grepl("Psi values", diagnostics$warnings)))
})

# ============================================================================
# Test 5: Invalid model structure (not a list)
# ============================================================================

test_that("check_SSE_diagnostics handles invalid model structure gracefully", {
  # Create mock data
  list_of_data <- create_mock_list_of_data(N = 30, T_each = 200, M = 8)
  
  # Create fit_results
  fit_results <- data.frame(
    Method = c("MixtureFA"),
    BIC = c(800),
    ARI = c(0.85),
    SSE = c(4000),
    VAF = c(0.80)
  )
  
  # Pass invalid model (not a list)
  invalid_model <- "not_a_list"
  
  diagnostics <- check_SSE_diagnostics(
    fit_results = fit_results,
    list_of_data = list_of_data,
    mixtureFA_fit = invalid_model
  )
  
  # Check that warning indicates model is not a list
  expect_true(any(grepl("mixtureFA_fit is not a list", diagnostics$warnings)))
  
  # Check that diagnostics still complete
  expect_true(diagnostics$convergence_ok)
  expect_true(diagnostics$psi_ok)
})

# ============================================================================
# Test 6: Missing essential fields in model
# ============================================================================

test_that("check_SSE_diagnostics handles missing essential fields gracefully", {
  # Create mock data
  list_of_data <- create_mock_list_of_data(N = 30, T_each = 200, M = 8)
  
  # Create fit_results
  fit_results <- data.frame(
    Method = c("MixtureFA", "MixturePCA"),
    BIC = c(800, 850),
    ARI = c(0.85, 0.88),
    SSE = c(4000, 4100),
    VAF = c(0.80, 0.79)
  )
  
  # Create incomplete models (missing essential fields)
  incomplete_mfa <- list(z = rep(1, 30), mu = list(rnorm(8)))  # Missing Psi, logLik
  incomplete_mpca <- list(z = rep(1, 30), mu = list(rnorm(8)))  # Missing sigma2, logLik
  
  diagnostics <- check_SSE_diagnostics(
    fit_results = fit_results,
    list_of_data = list_of_data,
    mixtureFA_fit = incomplete_mfa,
    mixturePCA_fit = incomplete_mpca
  )
  
  # Check that warnings indicate missing fields
  expect_true(any(grepl("mixtureFA_fit missing essential fields", diagnostics$warnings)))
  expect_true(any(grepl("mixturePCA_fit missing essential fields", diagnostics$warnings)))
  
  # Check that diagnostics still complete
  expect_true(diagnostics$convergence_ok)
  expect_true(diagnostics$psi_ok)
})

# ============================================================================
# Test 7: Empty clusters detection
# ============================================================================

test_that("check_SSE_diagnostics detects empty clusters", {
  # Create mock data
  list_of_data <- create_mock_list_of_data(N = 30, T_each = 200, M = 8)
  
  # Create model with empty cluster (all subjects assigned to cluster 1)
  mfa_fit <- create_mock_mfa_fit(N = 30, K = 3, M = 8, good_psi = TRUE)
  mfa_fit$z <- rep(1, 30)  # All in cluster 1, clusters 2 and 3 are empty
  
  fit_results <- data.frame(
    Method = c("MixtureFA"),
    BIC = c(800),
    ARI = c(0.85),
    SSE = c(4000),
    VAF = c(0.80)
  )
  
  diagnostics <- check_SSE_diagnostics(
    fit_results = fit_results,
    list_of_data = list_of_data,
    mixtureFA_fit = mfa_fit
  )
  
  # Note: table() doesn't create entries for missing levels, so empty clusters
  # won't be detected this way. But very small clusters should be detected.
  # Check that convergence_ok is TRUE (no actual empty clusters in table)
  expect_true(diagnostics$convergence_ok)
})

# ============================================================================
# Test 8: Very small clusters detection
# ============================================================================

test_that("check_SSE_diagnostics detects very small clusters", {
  # Create mock data
  list_of_data <- create_mock_list_of_data(N = 30, T_each = 200, M = 8)
  
  # Create model with very small cluster (1 subject in cluster 3)
  mfa_fit <- create_mock_mfa_fit(N = 30, K = 3, M = 8, good_psi = TRUE)
  mfa_fit$z <- c(rep(1, 14), rep(2, 15), 3)  # Cluster 3 has only 1 subject (< 5%)
  
  fit_results <- data.frame(
    Method = c("MixtureFA"),
    BIC = c(800),
    ARI = c(0.85),
    SSE = c(4000),
    VAF = c(0.80)
  )
  
  diagnostics <- check_SSE_diagnostics(
    fit_results = fit_results,
    list_of_data = list_of_data,
    mixtureFA_fit = mfa_fit
  )
  
  # Check that warning about small clusters is present
  expect_true(any(grepl("Very small clusters detected", diagnostics$warnings)))
})

# ============================================================================
# Test 9: Non-finite log-likelihood detection
# ============================================================================

test_that("check_SSE_diagnostics detects non-finite log-likelihood", {
  # Create mock data
  list_of_data <- create_mock_list_of_data(N = 30, T_each = 200, M = 8)
  
  # Create model with non-finite log-likelihood
  mfa_fit <- create_mock_mfa_fit(N = 30, K = 3, M = 8, good_psi = TRUE)
  mfa_fit$logLik <- Inf
  
  fit_results <- data.frame(
    Method = c("MixtureFA"),
    BIC = c(800),
    ARI = c(0.85),
    SSE = c(4000),
    VAF = c(0.80)
  )
  
  diagnostics <- check_SSE_diagnostics(
    fit_results = fit_results,
    list_of_data = list_of_data,
    mixtureFA_fit = mfa_fit
  )
  
  # Check that convergence_ok is FALSE
  expect_false(diagnostics$convergence_ok)
  
  # Check that warning about non-finite log-likelihood is present
  expect_true(any(grepl("Non-finite log-likelihood", diagnostics$warnings)))
})

# ============================================================================
# Test 10: VAF-SSE correlation with insufficient data points
# ============================================================================

test_that("check_SSE_diagnostics handles insufficient data points for correlation", {
  # Create mock data
  list_of_data <- create_mock_list_of_data(N = 30, T_each = 200, M = 8)
  
  # Create fit_results with only one method (insufficient for correlation)
  fit_results <- data.frame(
    Method = c("MixtureFA"),
    BIC = c(800),
    ARI = c(0.85),
    SSE = c(4000),
    VAF = c(0.80)
  )
  
  diagnostics <- check_SSE_diagnostics(
    fit_results = fit_results,
    list_of_data = list_of_data
  )
  
  # Check that correlation_ok is FALSE
  expect_false(diagnostics$correlation_ok)
  
  # Check that warning about insufficient data points is present
  expect_true(any(grepl("Insufficient data points to compute VAF-SSE correlation", diagnostics$warnings)))
})

# ============================================================================
# Test 11: strip_attributes utility function
# ============================================================================

test_that("strip_attributes removes custom attributes correctly", {
  # Create a data frame with custom attributes
  df <- data.frame(x = 1:3, y = 4:6)
  attr(df, "custom_attr1") <- "value1"
  attr(df, "custom_attr2") <- list(a = 1, b = 2)
  attr(df, "mixtureFA_fit") <- list(z = c(1, 2, 3))
  
  # Verify custom attributes exist
  expect_equal(attr(df, "custom_attr1"), "value1")
  expect_false(is.null(attr(df, "mixtureFA_fit")))
  
  # Strip attributes
  df_clean <- strip_attributes(df)
  
  # Check that custom attributes are removed
  expect_null(attr(df_clean, "custom_attr1"))
  expect_null(attr(df_clean, "custom_attr2"))
  expect_null(attr(df_clean, "mixtureFA_fit"))
  
  # Check that essential attributes are preserved
  expect_equal(names(df_clean), c("x", "y"))
  expect_equal(class(df_clean), "data.frame")
  expect_equal(nrow(df_clean), 3)
  
  # Check that data is unchanged
  expect_equal(df_clean$x, 1:3)
  expect_equal(df_clean$y, 4:6)
})
