# Tests for fit_two_step_baseline function
# Updated to use Dynamic v2 API (simulate_dynamic_synergy_data)

test_that("fit_two_step_baseline works with PCA method", {
  skip_if_not_installed("mclust")
  
  # Generate test data using Dynamic v2 simulator
  set.seed(123)
  N <- 10
  K_true <- 2
  r_true <- 2
  M <- 4
  T_each <- 20
  
  sim_data <- simulate_dynamic_synergy_data(
    N = N, K = K_true, r = r_true, M = M, T_each = T_each,
    cluster_sep_spatial = 2.0,
    cluster_sep_temporal = 2.0,
    cluster_sep_stability = 2.0,
    seed = 123
  )
  
  # Test PCA method
  result <- fit_two_step_baseline(
    sim_data$X_list,
    r = r_true,
    K = K_true,
    method = "PCA",
    seed = 123,
    z_true = sim_data$z_true
  )
  
  # Check return structure
  expect_type(result, "list")
  expect_true("z_est" %in% names(result))
  expect_true("W_list" %in% names(result))
  expect_true("W_centroids" %in% names(result))
  expect_true("mu_list" %in% names(result))
  expect_true("sd_list" %in% names(result))
  expect_true("method" %in% names(result))
  expect_true("BIC" %in% names(result))
  expect_true("ARI" %in% names(result))
  
  # Check mu_list and sd_list
  expect_type(result$mu_list, "list")
  expect_type(result$sd_list, "list")
  expect_length(result$mu_list, N)
  expect_length(result$sd_list, N)
  expect_equal(length(result$mu_list[[1]]), M)
  expect_equal(length(result$sd_list[[1]]), M)
  
  # Check z_est
  expect_type(result$z_est, "integer")
  expect_length(result$z_est, N)
  expect_true(all(result$z_est >= 1 & result$z_est <= K_true))
  
  # Check W_list
  expect_type(result$W_list, "list")
  expect_length(result$W_list, N)
  expect_equal(dim(result$W_list[[1]]), c(M, r_true))
  
  # Check W_centroids
  expect_true(is.matrix(result$W_centroids))
  expect_equal(dim(result$W_centroids), c(M, K_true))
  
  # Check method
  expect_equal(result$method, "TwoStep_PCA")
  
  # Check BIC is NA (k-means doesn't have BIC)
  expect_true(is.na(result$BIC))
  
  # Check ARI is computed
  expect_type(result$ARI, "double")
  expect_true(!is.na(result$ARI))
  expect_true(result$ARI >= -1 && result$ARI <= 1)
})

test_that("fit_two_step_baseline works with FA method when psych is available", {
  skip_if_not_installed("mclust")
  skip_if_not_installed("psych")
  
  # Generate test data using Dynamic v2 simulator
  set.seed(456)
  N <- 10
  K_true <- 2
  r_true <- 2
  M <- 4
  T_each <- 20
  
  sim_data <- simulate_dynamic_synergy_data(
    N = N, K = K_true, r = r_true, M = M, T_each = T_each,
    cluster_sep_spatial = 2.0,
    cluster_sep_temporal = 2.0,
    cluster_sep_stability = 2.0,
    seed = 456
  )
  
  # Test FA method
  result <- fit_two_step_baseline(
    sim_data$X_list,
    r = r_true,
    K = K_true,
    method = "FA",
    seed = 456,
    z_true = sim_data$z_true
  )
  
  # Check return structure
  expect_type(result, "list")
  expect_equal(result$method, "TwoStep_FA")
  expect_type(result$z_est, "integer")
  expect_length(result$z_est, N)
})

test_that("fit_two_step_baseline works without z_true", {
  skip_if_not_installed("mclust")
  
  # Generate test data using Dynamic v2 simulator
  set.seed(789)
  N <- 10
  K_true <- 2
  r_true <- 2
  M <- 4
  T_each <- 20
  
  sim_data <- simulate_dynamic_synergy_data(
    N = N, K = K_true, r = r_true, M = M, T_each = T_each,
    cluster_sep_spatial = 2.0,
    cluster_sep_temporal = 2.0,
    cluster_sep_stability = 2.0,
    seed = 789
  )
  
  # Test without z_true
  result <- fit_two_step_baseline(
    sim_data$X_list,
    r = r_true,
    K = K_true,
    method = "PCA",
    seed = 789,
    z_true = NULL
  )
  
  # Check ARI is NA when z_true is not provided
  expect_true(is.na(result$ARI))
  
  # But z_est should still be computed
  expect_type(result$z_est, "integer")
  expect_length(result$z_est, N)
})

test_that("plot_fig1_baseline_comparison creates valid plot", {
  skip_if_not_installed("ggplot2")
  
  # Create sample data
  df_test <- data.frame(
    Method = c("SingleFA", "SinglePCA", "TwoStep_FA", "TwoStep_PCA", 
               "MixtureFA", "MixturePCA"),
    BIC = c(100, 110, 90, 95, 80, 85),
    ARI = c(NA, NA, 0.5, 0.6, 0.8, 0.85),
    SSE = c(50, 55, 40, 42, 30, 32),
    VAF = c(0.5, 0.45, 0.6, 0.58, 0.7, 0.68)
  )
  
  # Test ARI plot
  p_ari <- plot_fig1_baseline_comparison(df_test, metric = "ARI")
  expect_s3_class(p_ari, "ggplot")
  
  # Test VAF plot
  p_vaf <- plot_fig1_baseline_comparison(df_test, metric = "VAF")
  expect_s3_class(p_vaf, "ggplot")
})

test_that("calc_reconstruction_error_twostep computes SSE correctly", {
  # Generate test data using Dynamic v2 simulator
  set.seed(999)
  N <- 5
  K_true <- 2
  r_true <- 2
  M <- 4
  T_each <- 20
  
  sim_data <- simulate_dynamic_synergy_data(
    N = N, K = K_true, r = r_true, M = M, T_each = T_each,
    cluster_sep_spatial = 2.0,
    cluster_sep_temporal = 2.0,
    cluster_sep_stability = 2.0,
    seed = 999
  )
  
  # Fit two-step baseline
  result <- fit_two_step_baseline(
    sim_data$X_list,
    r = r_true,
    K = K_true,
    method = "PCA",
    seed = 999,
    z_true = sim_data$z_true
  )
  
  # Compute SSE
  sse <- calc_reconstruction_error_twostep(sim_data$X_list, result)
  
  # Check SSE is positive
  expect_type(sse, "double")
  expect_true(sse > 0)
  expect_true(is.finite(sse))
})
