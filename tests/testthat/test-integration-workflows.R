# Integration Tests for Complete Workflows
# Tests for end-to-end analysis pipelines

# Load test fixtures
small_data <- readRDS("fixtures/small_test_data.rds")
medium_data <- readRDS("fixtures/medium_test_data.rds")

test_that("Complete MFA workflow executes successfully", {
  # Complete workflow: data -> fit -> evaluate -> visualize
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Step 1: Fit MFA model
  fit_mfa <- mfa_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 10,
    nIterFA = 5,
    tol = 1e-3,
    n_init = 2,
    use_kmeans_init = TRUE,
    mc_cores = 1,
    seed = 123
  )
  
  expect_type(fit_mfa, "list")
  expect_true(is.finite(fit_mfa$logLik))
  
  # Step 2: Compute factor scores
  fit_mfa <- compute_factor_scores_mfa(fit_mfa, list_of_data)
  
  expect_true("factor_scores" %in% names(fit_mfa))
  expect_equal(length(fit_mfa$factor_scores), true_params$N)
  
  # Step 3: Evaluate model performance
  vaf <- compute_global_vaf_mfa(list_of_data, fit_mfa)
  
  expect_true(vaf >= 0)
  expect_true(vaf <= 1)
  
  # Step 4: Get cluster sizes
  sizes <- compute_cluster_sizes(fit_mfa)
  
  expect_equal(sum(sizes), true_params$N)
  
  # Workflow completed successfully
  expect_true(TRUE)
})

test_that("Complete MPCA workflow executes successfully", {
  # Complete workflow: data -> fit -> evaluate
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Step 1: Fit MPCA model
  fit_mpca <- mixture_pca_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 10,
    nIterPCA = 5,
    tol = 1e-3,
    method = "EM",
    n_init = 2,
    use_kmeans_init = TRUE,
    mc_cores = 1,
    seed = 456
  )
  
  expect_type(fit_mpca, "list")
  expect_true(is.finite(fit_mpca$logLik))
  
  # Step 2: Compute factor scores
  fit_mpca <- compute_factor_scores_mpca(fit_mpca, list_of_data)
  
  expect_true("factor_scores" %in% names(fit_mpca))
  expect_equal(length(fit_mpca$factor_scores), true_params$N)
  
  # Step 3: Evaluate model performance
  vaf <- compute_global_vaf_mpca(list_of_data, fit_mpca)
  
  expect_true(vaf >= 0)
  expect_true(vaf <= 1)
  
  # Step 4: Get cluster sizes
  sizes <- compute_cluster_sizes_mpca(fit_mpca)
  
  expect_equal(sum(sizes), true_params$N)
  
  # Workflow completed successfully
  expect_true(TRUE)
})

test_that("Model selection workflow executes successfully", {
  # Complete model selection workflow
  list_of_data <- small_data$data
  
  # Step 1: Grid search for optimal K and r
  result_mfa <- select_optimal_K_r_mfa(
    list_of_data = list_of_data,
    Kvec = 1:2,
    rvec = 1:2,
    max_iter = 5,
    nIterFA = 3,
    tol = 1e-3,
    n_init = 1,
    use_kmeans_init = FALSE
  )
  
  expect_type(result_mfa, "list")
  expect_true("best_model_info" %in% names(result_mfa))
  
  # Step 2: Extract best model
  best_model <- result_mfa$best_model_info$model
  
  expect_type(best_model, "list")
  expect_true(is.finite(best_model$logLik))
  
  # Step 3: Post-hoc evaluation
  df_posthoc <- posthoc_mfa_evaluation(list_of_data, result_mfa)
  
  expect_s3_class(df_posthoc, "data.frame")
  expect_true("GlobalVAF" %in% names(df_posthoc))
  
  # Step 4: Retrieve specific model
  retrieved_model <- get_model_by_K_r(
    result_mfa, 
    K_target = result_mfa$best_model_info$K,
    r_target = result_mfa$best_model_info$r
  )
  
  expect_type(retrieved_model, "list")
  
  # Workflow completed successfully
  expect_true(TRUE)
})

test_that("MFA vs MPCA comparison workflow executes successfully", {
  skip_if_not_installed("mclust")

  # Compare MFA and MPCA on same data
  list_of_data <- small_data$data
  true_params <- small_data$true_params

  # Fit both models with same K and r
  fit_mfa <- mfa_em_fit(
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
  
  fit_mpca <- mixture_pca_em_fit(
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
  
  # Compare log-likelihoods
  expect_true(is.finite(fit_mfa$logLik))
  expect_true(is.finite(fit_mpca$logLik))
  
  # Compare VAF
  vaf_mfa <- compute_global_vaf_mfa(list_of_data, fit_mfa)
  vaf_mpca <- compute_global_vaf_mpca(list_of_data, fit_mpca)
  
  expect_true(vaf_mfa >= 0 && vaf_mfa <= 1)
  expect_true(vaf_mpca >= 0 && vaf_mpca <= 1)
  
  # Compare cluster assignments (using ARI)
  ari <- mclust::adjustedRandIndex(fit_mfa$z, fit_mpca$z)
  
  expect_true(ari >= -1 && ari <= 1)
  
  # Both methods should produce valid results
  expect_true(TRUE)
})

test_that("Rotation and alignment workflow executes successfully", {
  # Workflow with rotation and sign alignment
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Step 1: Fit MFA model
  fit_mfa <- mfa_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 10,
    nIterFA = 5,
    tol = 1e-3,
    n_init = 1,
    use_kmeans_init = FALSE,
    seed = 789
  )
  
  # Step 2: Compute factor scores
  fit_mfa <- compute_factor_scores_mfa(fit_mfa, list_of_data)
  
  # Step 3: Rotate loadings
  rotated_mfa <- rotate_mfa_model(
    mfa_fit = fit_mfa,
    rotation = "varimax",
    rotate_scores = TRUE
  )
  
  expect_type(rotated_mfa, "list")
  expect_true("Lambda" %in% names(rotated_mfa))
  
  # Step 4: Align factor signs
  aligned_mfa <- align_factor_signs_mfa(rotated_mfa)
  
  expect_type(aligned_mfa, "list")
  expect_true("Lambda" %in% names(aligned_mfa))
  
  # Check that dimensions are preserved
  for (k in 1:true_params$K) {
    expect_equal(dim(aligned_mfa$Lambda[[k]]), c(true_params$M, true_params$r))
  }
  
  # Workflow completed successfully
  expect_true(TRUE)
})

test_that("Multi-initialization workflow improves results", {
  # Test that multi-initialization finds better solutions
  list_of_data <- medium_data$data
  true_params <- medium_data$true_params
  
  # Single initialization
  fit_single <- mfa_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 10,
    nIterFA = 5,
    tol = 1e-3,
    n_init = 1,
    use_kmeans_init = FALSE,
    seed = 111
  )
  
  # Multiple initializations
  fit_multi <- mfa_em_fit(
    list_of_data = list_of_data,
    K = true_params$K,
    r = true_params$r,
    max_iter = 10,
    nIterFA = 5,
    tol = 1e-3,
    n_init = 5,
    use_kmeans_init = TRUE,
    mc_cores = 1,
    seed = 111
  )
  
  # Multi-init should have equal or better log-likelihood
  expect_true(fit_multi$logLik >= fit_single$logLik - 1e-6)
  
  # Both should produce valid results
  expect_true(is.finite(fit_single$logLik))
  expect_true(is.finite(fit_multi$logLik))
})

test_that("Complete workflow handles edge cases", {
  # Test workflow with edge case data
  edge_cases <- readRDS("fixtures/edge_case_data.rds")
  
  # Single cluster case
  single_cluster_data <- edge_cases$single_cluster$data
  
  fit_single_cluster <- mfa_em_fit(
    list_of_data = single_cluster_data,
    K = 1,
    r = 2,
    max_iter = 5,
    nIterFA = 3,
    tol = 1e-3,
    n_init = 1,
    use_kmeans_init = FALSE
  )
  
  expect_equal(length(unique(fit_single_cluster$z)), 1)
  expect_true(is.finite(fit_single_cluster$logLik))
  
  # Single factor case
  single_factor_data <- edge_cases$single_factor$data
  
  fit_single_factor <- mfa_em_fit(
    list_of_data = single_factor_data,
    K = 2,
    r = 1,
    max_iter = 5,
    nIterFA = 3,
    tol = 1e-3,
    n_init = 1,
    use_kmeans_init = FALSE
  )
  
  expect_equal(ncol(fit_single_factor$Lambda[[1]]), 1)
  expect_true(is.finite(fit_single_factor$logLik))
  
  # Workflows completed successfully
  expect_true(TRUE)
})

test_that("Workflow produces consistent results across runs", {
  # Test reproducibility of complete workflow
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Run complete workflow twice
  run_workflow <- function(seed_val) {
    fit <- mfa_em_fit(
      list_of_data = list_of_data,
      K = true_params$K,
      r = true_params$r,
      max_iter = 10,
      nIterFA = 5,
      tol = 1e-3,
      n_init = 2,
      use_kmeans_init = TRUE,
      mc_cores = 1,
      seed = seed_val
    )
    
    fit <- compute_factor_scores_mfa(fit, list_of_data)
    vaf <- compute_global_vaf_mfa(list_of_data, fit)
    
    list(fit = fit, vaf = vaf)
  }
  
  result1 <- run_workflow(999)
  result2 <- run_workflow(999)
  
  # Results should be identical with same seed
  expect_equal(result1$fit$z, result2$fit$z)
  expect_equal(result1$fit$logLik, result2$fit$logLik, tolerance = 1e-10)
  expect_equal(result1$vaf, result2$vaf, tolerance = 1e-10)
})

test_that("Workflow handles missing or invalid data gracefully", {
  # Test error handling in workflow
  
  # Empty data list
  expect_error(
    mfa_em_fit(list_of_data = list(), K = 2, r = 2)
  )
  
  # Invalid K (will fail in kmeans with specific message)
  expect_error(
    mfa_em_fit(list_of_data = small_data$data, K = 0, r = 2)
  )
  
  # Invalid r (will fail in C++ code)
  expect_error(
    mfa_em_fit(list_of_data = small_data$data, K = 2, r = 0)
  )
})
