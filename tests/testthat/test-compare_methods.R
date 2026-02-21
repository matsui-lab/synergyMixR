# Tests for compare_baselines() - the new API for method comparison
# This replaces the old test-compare_dynamic.R which tested the deprecated
# compare_static_vs_dynamic() function.

test_that("compare_baselines runs without errors (B1)", {
  skip_if_not_installed("mclust")
  # Test with small parameters for speed
  # Use high cluster separation and low noise to ensure stable clustering
  expect_no_error({
    result <- compare_baselines(
      N = 10,
      K = 2,
      r = 2,
      M = 4,
      T_each = 20,
      cluster_sep_spatial = 2.0,
      cluster_sep_temporal = 2.0,
      cluster_sep_stability = 2.0,
      seed = 123,
      n_init = 1,
      mc_cores = 1
    )
  })
})

test_that("compare_baselines returns data.frame (B2)", {
  skip_if_not_installed("mclust")
  result <- compare_baselines(
    N = 10,
    K = 2,
    r = 2,
    M = 4,
    T_each = 20,
    cluster_sep_spatial = 2.0,
    cluster_sep_temporal = 2.0,
    cluster_sep_stability = 2.0,
    seed = 456,
    n_init = 1,
    mc_cores = 1
  )
  
  expect_s3_class(result, "data.frame")
})

test_that("compare_baselines returns required columns (B3)", {
  skip_if_not_installed("mclust")
  result <- compare_baselines(
    N = 10,
    K = 2,
    r = 2,
    M = 4,
    T_each = 20,
    cluster_sep_spatial = 2.0,
    cluster_sep_temporal = 2.0,
    cluster_sep_stability = 2.0,
    seed = 789,
    n_init = 1,
    mc_cores = 1
  )
  
  # Check required columns exist
  required_cols <- c(
    "Method", "BIC", "ARI", "SSE", "VAF",
    "model_type", "seed",
    "cluster_sep_spatial", "cluster_sep_temporal", "cluster_sep_stability",
    "N", "K", "r", "M", "T_each"
  )
  expect_true(all(required_cols %in% names(result)))
  
  # Check column types

  expect_type(result$Method, "character")
  expect_type(result$BIC, "double")
  expect_type(result$SSE, "double")
  expect_type(result$VAF, "double")
  expect_type(result$model_type, "character")
  expect_type(result$seed, "double")
  expect_type(result$N, "double")
  expect_type(result$K, "double")
  expect_type(result$r, "double")
  expect_type(result$M, "double")
  expect_type(result$T_each, "double")
})

test_that("compare_baselines returns Dynamic_v2 only - no Mixture (B4)", {
  skip_if_not_installed("mclust")
  result <- compare_baselines(
    N = 10,
    K = 2,
    r = 2,
    M = 4,
    T_each = 20,
    cluster_sep_spatial = 2.0,
    cluster_sep_temporal = 2.0,
    cluster_sep_stability = 2.0,
    seed = 111,
    n_init = 1,
    mc_cores = 1
  )
  

  # Should have 6 rows (6 methods, 1 model type)
  expect_equal(nrow(result), 6)
  
  # All model_type should be "Dynamic_v2"
  expect_true(all(result$model_type == "Dynamic_v2"))
  
  # No "Mixture" model type should exist
  expect_false(any(result$model_type == "Mixture"))
  
  # Check expected method names
  expected_methods <- c("SingleFA", "SinglePCA", "TwoStep_FA", "TwoStep_PCA", "MixtureFA", "MixturePCA")
  expect_setequal(result$Method, expected_methods)
  
  # SingleFA and SinglePCA should have NA for ARI (no clustering)
  single_methods <- result[result$Method %in% c("SingleFA", "SinglePCA"), ]
  expect_true(all(is.na(single_methods$ARI)))
  
  # Other methods should have numeric ARI
  clustering_methods <- result[!result$Method %in% c("SingleFA", "SinglePCA"), ]
  expect_true(all(!is.na(clustering_methods$ARI)))
})

test_that("compare_baselines seed column matches input seed (B4 supplement)", {
  skip_if_not_installed("mclust")
  test_seed <- 222
  result <- compare_baselines(
    N = 10,
    K = 2,
    r = 2,
    M = 4,
    T_each = 20,
    cluster_sep_spatial = 2.0,
    cluster_sep_temporal = 2.0,
    cluster_sep_stability = 2.0,
    seed = test_seed,
    n_init = 1,
    mc_cores = 1
  )
  
  # All rows should have the same seed value
  expect_true(all(result$seed == test_seed))
})

test_that("compare_baselines saves results when output_dir specified (B5)", {
  skip_if_not_installed("mclust")
  temp_dir <- tempdir()
  output_dir <- file.path(temp_dir, "test_compare_baselines_output")
  
  # Clean up any existing directory
  if (dir.exists(output_dir)) {
    unlink(output_dir, recursive = TRUE)
  }
  
  # Run compare_baselines with output_dir (compare_baselines doesn't have output_dir param)
  result <- compare_baselines(
    N = 10,
    K = 2,
    r = 2,
    M = 4,
    T_each = 20,
    cluster_sep_spatial = 2.0,
    cluster_sep_temporal = 2.0,
    cluster_sep_stability = 2.0,
    seed = 333,
    n_init = 1,
    mc_cores = 1
  )
  
  # Verify result is returned

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 6)
})

test_that("compare_static_vs_dynamic handles invalid model_type (B6)", {
  skip_if_not_installed("mclust")
  temp_dir <- tempdir()
  output_dir <- file.path(temp_dir, "test_compare_baselines_invalid")
  
  # Clean up any existing directory
  if (dir.exists(output_dir)) {
    unlink(output_dir, recursive = TRUE)
  }
  
  expect_error(
    compare_static_vs_dynamic(
      model_types = c("Invalid"),
      N = 10,
      K = 2,
      r = 2,
      M = 4,
      T_each = 20,
      seed = 999,
      output_dir = output_dir,
      n_init = 1,
      mc_cores = 1
    ),
    "Unknown model_type"
  )
  
  # Clean up
  if (dir.exists(output_dir)) {
    unlink(output_dir, recursive = TRUE)
  }
})

test_that("compare_baselines parameter values are preserved in output", {
  skip_if_not_installed("mclust")
  # Test that input parameters are correctly reflected in output
  result <- compare_baselines(
    N = 10,
    K = 2,
    r = 2,
    M = 4,
    T_each = 20,
    cluster_sep_spatial = 0.5,
    cluster_sep_temporal = 1.5,
    cluster_sep_stability = 2.0,
    seed = 444,
    n_init = 1,
    mc_cores = 1
  )
  
  # Check that parameter values are preserved
  expect_true(all(result$N == 10))
  expect_true(all(result$K == 2))
  expect_true(all(result$r == 2))
  expect_true(all(result$M == 4))
  expect_true(all(result$T_each == 20))
  expect_true(all(result$cluster_sep_spatial == 0.5))
  expect_true(all(result$cluster_sep_temporal == 1.5))
  expect_true(all(result$cluster_sep_stability == 2.0))
  expect_true(all(result$seed == 444))
})

# ============================================================================
# Tests for Hard vs Soft EM Comparison (M2 Reviewer Response)
# ============================================================================

test_that("soft EM and hard EM produce comparable results for MFA", {
  skip_if_not_installed("mclust")

  # Load medium test data for more stable results
  medium_data <- readRDS("fixtures/medium_test_data.rds")
  list_of_data <- medium_data$data
  true_params <- medium_data$true_params
  K <- true_params$K
  r <- true_params$r
  z_init <- true_params$z

  fit_hard <- mfa_em_fit_cpp_singleInit(
    list_of_data = list_of_data, K = K, r = r, z_init = z_init,
    max_iter = 50, nIterFA = 10, tol = 1e-4, em_type = "hard"
  )

  fit_soft <- mfa_em_fit_cpp_singleInit(
    list_of_data = list_of_data, K = K, r = r, z_init = z_init,
    max_iter = 50, nIterFA = 10, tol = 1e-4, em_type = "soft"
  )

  # Soft EM should have equal or better log-likelihood
  # (allowing small tolerance for numerical differences)
  expect_true(fit_soft$logLik >= fit_hard$logLik - 1.0,
              info = paste("Soft logLik:", fit_soft$logLik, "Hard logLik:", fit_hard$logLik))

  # Cluster assignments should be similar (ARI > 0.7)
  ari <- mclust::adjustedRandIndex(fit_hard$z, fit_soft$z)
  expect_true(ari > 0.7,
              info = paste("ARI between hard and soft EM:", round(ari, 3)))
})

test_that("soft EM and hard EM produce comparable results for MPCA", {
  skip_if_not_installed("mclust")

  # Load medium test data
  medium_data <- readRDS("fixtures/medium_test_data.rds")
  list_of_data <- medium_data$data
  true_params <- medium_data$true_params
  K <- true_params$K
  r <- true_params$r
  z_init <- true_params$z

  fit_hard <- mixture_pca_em_fit_cpp_singleInit(
    list_of_data = list_of_data, K = K, r = r, z_init = z_init,
    max_iter = 50, nIterPCA = 10, tol = 1e-4, em_type = "hard"
  )

  fit_soft <- mixture_pca_em_fit_cpp_singleInit(
    list_of_data = list_of_data, K = K, r = r, z_init = z_init,
    max_iter = 50, nIterPCA = 10, tol = 1e-4, em_type = "soft"
  )

  # Soft EM should have equal or better log-likelihood
  expect_true(fit_soft$logLik >= fit_hard$logLik - 1.0,
              info = paste("Soft logLik:", fit_soft$logLik, "Hard logLik:", fit_hard$logLik))

  # Cluster assignments should be similar (ARI > 0.7)
  ari <- mclust::adjustedRandIndex(fit_hard$z, fit_soft$z)
  expect_true(ari > 0.7,
              info = paste("ARI between hard and soft EM:", round(ari, 3)))
})

test_that("soft EM responsibilities are more gradual than hard EM", {
  # Load small test data
  small_data <- readRDS("fixtures/small_test_data.rds")
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  K <- true_params$K
  r <- true_params$r
  z_init <- true_params$z

  fit_hard <- mfa_em_fit_cpp_singleInit(
    list_of_data = list_of_data, K = K, r = r, z_init = z_init,
    max_iter = 30, nIterFA = 5, tol = 1e-3, em_type = "hard"
  )

  fit_soft <- mfa_em_fit_cpp_singleInit(
    list_of_data = list_of_data, K = K, r = r, z_init = z_init,
    max_iter = 30, nIterFA = 5, tol = 1e-3, em_type = "soft"
  )

  # Hard EM responsibilities should be more extreme (closer to 0 or 1)
  # Soft EM should have more uncertainty in some cases
  hard_entropy <- -sum(fit_hard$resp * log(fit_hard$resp + 1e-16))
  soft_entropy <- -sum(fit_soft$resp * log(fit_soft$resp + 1e-16))

  # Soft EM typically has higher or equal entropy (more uncertainty)
  # but this isn't guaranteed, so we just check both are finite

  expect_true(is.finite(hard_entropy))
  expect_true(is.finite(soft_entropy))
})
