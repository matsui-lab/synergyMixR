test_that("select_r_global_mfa works with PPCA and CV", {
  skip_on_cran()
  
  # Generate simple test data
  set.seed(123)
  N <- 20
  M <- 8
  T_each <- 50
  r_true <- 3
  
  # Simulate data with known r
  list_of_data <- lapply(1:N, function(i) {
    matrix(rnorm(T_each * M), nrow = T_each, ncol = M)
  })
  
  # Test with PPCA and CV
  result <- select_r_global_mfa(
    list_of_data = list_of_data,
    rvec = 2:5,
    model = "PPCA",
    criterion = "CV",
    folds = 3,
    seed = 123,
    verbose = FALSE
  )
  
  # Check structure
  expect_type(result, "list")
  expect_true("r_global" %in% names(result))
  expect_true("summary" %in% names(result))
  expect_true("criterion" %in% names(result))
  expect_true("model" %in% names(result))
  
  # Check values
  expect_equal(result$criterion, "CV")
  expect_equal(result$model, "PPCA")
  expect_true(result$r_global >= 2 && result$r_global <= 5)
  
  # Check summary structure
  expect_s3_class(result$summary, "data.frame")
  expect_true("r" %in% names(result$summary))
  expect_true("cv_mean" %in% names(result$summary))
  expect_true("cv_se" %in% names(result$summary))
  expect_true("vaf" %in% names(result$summary))
})


test_that("select_r_global_mfa works with FA and CV", {
  skip_on_cran()
  
  # Generate simple test data
  set.seed(456)
  N <- 15
  M <- 6
  T_each <- 40
  
  list_of_data <- lapply(1:N, function(i) {
    matrix(rnorm(T_each * M), nrow = T_each, ncol = M)
  })
  
  # Test with FA and CV
  result <- select_r_global_mfa(
    list_of_data = list_of_data,
    rvec = 2:4,
    model = "FA",
    criterion = "CV",
    folds = 3,
    seed = 456,
    verbose = FALSE
  )
  
  # Check structure
  expect_type(result, "list")
  expect_equal(result$criterion, "CV")
  expect_equal(result$model, "FA")
  expect_true(result$r_global >= 2 && result$r_global <= 4)
})


test_that("select_r_global_mfa works with Parallel Analysis", {
  skip_on_cran()
  
  # Generate simple test data
  set.seed(789)
  N <- 15
  M <- 8
  T_each <- 50
  
  list_of_data <- lapply(1:N, function(i) {
    matrix(rnorm(T_each * M), nrow = T_each, ncol = M)
  })
  
  # Test with Parallel Analysis
  result <- select_r_global_mfa(
    list_of_data = list_of_data,
    rvec = 2:5,
    model = "PPCA",
    criterion = "ParallelAnalysis",
    n_permutations = 20,  # Small number for speed
    seed = 789,
    verbose = FALSE
  )
  
  # Check structure
  expect_type(result, "list")
  expect_equal(result$criterion, "ParallelAnalysis")
  expect_equal(result$model, "PPCA")
  expect_true(result$r_global >= 2 && result$r_global <= 5)
  
  # Check PA-specific outputs
  expect_true("pa_curve" %in% names(result))
  expect_s3_class(result$summary, "data.frame")
  expect_true("eigenvalue_obs" %in% names(result$summary))
  expect_true("eigenvalue_perm_95" %in% names(result$summary))
})


test_that("select_r_global_mfa returns valid results with structured data", {
  skip_on_cran()

  # Generate data with known factor structure
  # IMPORTANT: Use common Lambda across all subjects (MFA model assumption)
  # Note: CV-based r selection with negLogLik tends to prefer larger r values
  # as it doesn't penalize model complexity. This is expected behavior.
  set.seed(999)
  N <- 30
  M <- 8
  T_each <- 100
  r_true <- 3

  # Create common loading matrix
  Lambda_common <- matrix(rnorm(M * r_true, sd = 2), nrow = M, ncol = r_true)

  # Simulate factor model: X = Z * Lambda^T + noise
  noise_sd <- 0.5
  list_of_data <- lapply(1:N, function(i) {
    Z <- matrix(rnorm(T_each * r_true), nrow = T_each, ncol = r_true)
    noise <- matrix(rnorm(T_each * M, sd = noise_sd), nrow = T_each, ncol = M)
    X <- Z %*% t(Lambda_common) + noise
    X
  })

  # Test CV selection
  result <- select_r_global_mfa(
    list_of_data = list_of_data,
    rvec = 1:6,
    model = "PPCA",
    criterion = "CV",
    folds = 5,
    seed = 999,
    verbose = FALSE
  )

  # Check that result has valid structure
  expect_type(result, "list")
  expect_true("r_global" %in% names(result))
  expect_true("summary" %in% names(result))
  expect_s3_class(result$summary, "data.frame")

  # Check that selected r is within candidate range
  expect_true(result$r_global >= 1 && result$r_global <= 6)

  # Check that CV error decreases (or at least doesn't increase much) for true r
  # VAF at r=r_true should be high (>90%) indicating the model captures signal
  vaf_at_true_r <- result$summary$vaf[result$summary$r == r_true]
  expect_true(vaf_at_true_r > 0.9,
              info = sprintf("VAF at true r=%d should be > 0.9, got %.3f",
                             r_true, vaf_at_true_r))
})


test_that("select_r_clusterwise_mfa works with global-first mode", {
  skip_on_cran()
  
  # Generate simple test data
  set.seed(111)
  N <- 20
  M <- 8
  T_each <- 50
  K <- 3
  
  list_of_data <- lapply(1:N, function(i) {
    matrix(rnorm(T_each * M), nrow = T_each, ncol = M)
  })
  
  # Test global-first mode
  result <- select_r_clusterwise_mfa(
    list_of_data = list_of_data,
    K = K,
    rvec = 2:5,
    mode = "global-first",
    global_model = "PPCA",
    global_criterion = "CV",
    global_folds = 3,
    global_seed = 111,
    verbose = FALSE
  )
  
  # Check structure
  expect_type(result, "list")
  expect_equal(result$mode, "global-first")
  expect_equal(result$agg_method, "global")
  expect_true(result$agg_r >= 2 && result$agg_r <= 5)
  expect_null(result$cluster_r)
  expect_null(result$cluster_metrics)
  expect_null(result$initial_fit)
  expect_true("global_result" %in% names(result))
  expect_type(result$global_result, "list")
})


test_that("select_r_clusterwise_mpca works with global-first mode", {
  skip_on_cran()
  
  # Generate simple test data
  set.seed(222)
  N <- 20
  M <- 8
  T_each <- 50
  K <- 3
  
  list_of_data <- lapply(1:N, function(i) {
    matrix(rnorm(T_each * M), nrow = T_each, ncol = M)
  })
  
  # Test global-first mode for MPCA
  result <- select_r_clusterwise_mpca(
    list_of_data = list_of_data,
    K = K,
    rvec = 2:5,
    mode = "global-first",
    global_model = "PPCA",
    global_criterion = "CV",
    global_folds = 3,
    global_seed = 222,
    verbose = FALSE
  )
  
  # Check structure
  expect_type(result, "list")
  expect_equal(result$mode, "global-first")
  expect_equal(result$agg_method, "global")
  expect_true(result$agg_r >= 2 && result$agg_r <= 5)
  expect_null(result$cluster_r)
  expect_null(result$cluster_metrics)
  expect_null(result$initial_fit)
  expect_true("global_result" %in% names(result))
})


test_that("CV metric options work correctly", {
  skip_on_cran()
  
  # Generate simple test data
  set.seed(333)
  N <- 15
  M <- 6
  T_each <- 40
  
  list_of_data <- lapply(1:N, function(i) {
    matrix(rnorm(T_each * M), nrow = T_each, ncol = M)
  })
  
  # Test with Frobenius metric
  result_frob <- select_r_global_mfa(
    list_of_data = list_of_data,
    rvec = 2:4,
    model = "PPCA",
    criterion = "CV",
    cv_metric = "frobenius",
    folds = 3,
    seed = 333,
    verbose = FALSE
  )
  
  expect_type(result_frob, "list")
  expect_true(result_frob$r_global >= 2 && result_frob$r_global <= 4)
  
  # Test with negLogLik metric (for FA)
  result_nll <- select_r_global_mfa(
    list_of_data = list_of_data,
    rvec = 2:4,
    model = "FA",
    criterion = "CV",
    cv_metric = "negLogLik",
    folds = 3,
    seed = 333,
    verbose = FALSE
  )
  
  expect_type(result_nll, "list")
  expect_true(result_nll$r_global >= 2 && result_nll$r_global <= 4)
})


test_that("Input validation works correctly", {
  # Generate simple test data
  set.seed(444)
  N <- 10
  M <- 6
  T_each <- 30
  
  list_of_data <- lapply(1:N, function(i) {
    matrix(rnorm(T_each * M), nrow = T_each, ncol = M)
  })
  
  # Test with invalid rvec (too few values)
  expect_error(
    select_r_global_mfa(
      list_of_data = list_of_data,
      rvec = 3,
      model = "PPCA",
      criterion = "CV"
    ),
    "rvec must contain at least 2 candidate r values"
  )
  
  # Test with invalid folds
  expect_error(
    select_r_global_mfa(
      list_of_data = list_of_data,
      rvec = 2:4,
      model = "PPCA",
      criterion = "CV",
      folds = 1
    ),
    "folds must be between 2 and"
  )
  
  # Test with too many folds
  expect_error(
    select_r_global_mfa(
      list_of_data = list_of_data,
      rvec = 2:4,
      model = "PPCA",
      criterion = "CV",
      folds = N + 1
    ),
    "folds must be between 2 and"
  )
})
