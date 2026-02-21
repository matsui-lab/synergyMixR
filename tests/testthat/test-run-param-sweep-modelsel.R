# Tests for run_param_sweep_modelsel
# test-run-param-sweep-modelsel.R

test_that("run_param_sweep_modelsel completes in test mode", {
  skip_on_cran()
  skip_on_ci()
  
  temp_dir <- tempdir()
  output_dir <- file.path(temp_dir, "test_sweep_modelsel")
  
  result <- run_param_sweep_modelsel(
    test = TRUE,
    cores = 1,
    output_dir = output_dir,
    strategies = c("grid"),
    methods = c("MFA"),
    max_iter = 10,
    n_init = 1,
    use_kmeans_init = FALSE,
    save_selection_summaries = FALSE,
    verbose = FALSE
  )
  
  expect_type(result, "list")
  expect_true("raw_results" %in% names(result))
  expect_true("aggregated_results" %in% names(result))
  expect_true("win_rates" %in% names(result))
  expect_true("output_dir" %in% names(result))
  
  expect_true(nrow(result$raw_results) > 0)
  
  unlink(output_dir, recursive = TRUE)
})

test_that("run_param_sweep_modelsel raw_results has required columns", {
  skip_on_cran()
  skip_on_ci()
  
  temp_dir <- tempdir()
  output_dir <- file.path(temp_dir, "test_sweep_modelsel_cols")
  
  result <- run_param_sweep_modelsel(
    test = TRUE,
    cores = 1,
    output_dir = output_dir,
    strategies = c("grid"),
    methods = c("MFA"),
    max_iter = 10,
    n_init = 1,
    use_kmeans_init = FALSE,
    save_selection_summaries = FALSE,
    verbose = FALSE
  )
  
  required_cols <- c(
    "run_id", "seed",
    "N", "K_true", "r_true", "M", "T_each",
    "cluster_sep_spatial", "cluster_sep_temporal", "cluster_sep_stability",
    "method", "selection_strategy",
    "K_hat", "r_hat", "K_correct", "r_correct", "Kr_correct",
    "ARI", "VAF", "logLik", "BIC", "runtime_sec", "error_msg"
  )
  
  for (col in required_cols) {
    expect_true(col %in% names(result$raw_results),
                info = paste("Missing column:", col))
  }
  
  unlink(output_dir, recursive = TRUE)
})

test_that("run_param_sweep_modelsel creates output files", {
  skip_on_cran()
  skip_on_ci()
  
  temp_dir <- tempdir()
  output_dir <- file.path(temp_dir, "test_sweep_modelsel_files")
  
  result <- run_param_sweep_modelsel(
    test = TRUE,
    cores = 1,
    output_dir = output_dir,
    strategies = c("grid"),
    methods = c("MFA"),
    max_iter = 10,
    n_init = 1,
    use_kmeans_init = FALSE,
    save_selection_summaries = FALSE,
    verbose = FALSE
  )
  
  expect_true(dir.exists(output_dir))
  
  rds_files <- list.files(output_dir, pattern = "sweep_modelsel_raw_.*\\.rds$")
  expect_true(length(rds_files) > 0)
  
  expect_true(file.exists(file.path(output_dir, "sweep_modelsel_aggregated.csv")))
  
  expect_true(file.exists(file.path(output_dir, "win_rates.csv")))
  
  unlink(output_dir, recursive = TRUE)
})

test_that("run_param_sweep_modelsel handles multiple methods and strategies", {
  skip_on_cran()
  skip_on_ci()
  
  temp_dir <- tempdir()
  output_dir <- file.path(temp_dir, "test_sweep_modelsel_multi")
  
  result <- run_param_sweep_modelsel(
    test = TRUE,
    cores = 1,
    output_dir = output_dir,
    strategies = c("grid", "twostage"),
    methods = c("MFA", "MPCA"),
    max_iter = 10,
    n_init = 1,
    use_kmeans_init = FALSE,
    save_selection_summaries = FALSE,
    verbose = FALSE
  )
  
  expect_true("MFA" %in% result$raw_results$method)
  expect_true("MPCA" %in% result$raw_results$method)
  expect_true("grid" %in% result$raw_results$selection_strategy)
  expect_true("twostage" %in% result$raw_results$selection_strategy)
  
  unlink(output_dir, recursive = TRUE)
})

test_that("make_param_grid_modelsel_test creates valid grid", {
  grid <- make_param_grid_modelsel_test()
  
  expect_s3_class(grid, "data.frame")
  expect_true(nrow(grid) > 0)
  
  required_cols <- c("N", "K", "r", "M", "T_each",
                     "cluster_sep_spatial", "cluster_sep_temporal",
                     "cluster_sep_stability", "seed", "run_id")
  
  for (col in required_cols) {
    expect_true(col %in% names(grid),
                info = paste("Missing column:", col))
  }
})

test_that("make_param_grid_modelsel_full creates valid grid", {
  grid <- make_param_grid_modelsel_full()
  
  expect_s3_class(grid, "data.frame")
  expect_true(nrow(grid) > 0)
  
  expect_true(nrow(grid) > nrow(make_param_grid_modelsel_test()))
})

test_that("aggregate_modelsel_results produces valid output", {
  skip_on_cran()
  
  raw_df <- data.frame(
    N = c(20, 20, 20, 20),
    K_true = c(2, 2, 2, 2),
    r_true = c(2, 2, 2, 2),
    M = c(8, 8, 8, 8),
    T_each = c(100, 100, 100, 100),
    cluster_sep_spatial = c(1.0, 1.0, 1.0, 1.0),
    cluster_sep_temporal = c(1.0, 1.0, 1.0, 1.0),
    cluster_sep_stability = c(1.0, 1.0, 1.0, 1.0),
    method = c("MFA", "MFA", "MPCA", "MPCA"),
    selection_strategy = c("grid", "twostage", "grid", "twostage"),
    ARI = c(0.8, 0.85, 0.75, 0.78),
    VAF = c(0.9, 0.92, 0.88, 0.89),
    BIC = c(1000, 950, 1100, 1050),
    subspace_similarity = c(0.85, 0.88, 0.82, 0.84),
    runtime_sec = c(10, 5, 12, 6),
    K_correct = c(TRUE, TRUE, FALSE, TRUE),
    r_correct = c(TRUE, FALSE, TRUE, TRUE),
    Kr_correct = c(TRUE, FALSE, FALSE, TRUE),
    error_msg = c("", "", "", ""),
    stringsAsFactors = FALSE
  )

  agg <- aggregate_modelsel_results(raw_df)

  expect_s3_class(agg, "data.frame")
  expect_true(nrow(agg) > 0)

  expect_true("ARI_mean" %in% names(agg))
  expect_true("VAF_mean" %in% names(agg))
  expect_true("subspace_similarity_mean" %in% names(agg))
  expect_true("K_correct_rate" %in% names(agg))
  expect_true("r_correct_rate" %in% names(agg))
  expect_true("Kr_correct_rate" %in% names(agg))
})

test_that("compute_modelsel_win_rates produces valid output", {
  skip_on_cran()
  
  agg_df <- data.frame(
    N = c(20, 20, 20, 20),
    K_true = c(2, 2, 2, 2),
    r_true = c(2, 2, 2, 2),
    M = c(8, 8, 8, 8),
    T_each = c(100, 100, 100, 100),
    cluster_sep_spatial = c(1.0, 1.0, 1.0, 1.0),
    cluster_sep_temporal = c(1.0, 1.0, 1.0, 1.0),
    cluster_sep_stability = c(1.0, 1.0, 1.0, 1.0),
    method = c("MFA", "MFA", "MPCA", "MPCA"),
    selection_strategy = c("grid", "twostage", "grid", "twostage"),
    ARI_mean = c(0.8, 0.85, 0.75, 0.78),
    VAF_mean = c(0.9, 0.92, 0.88, 0.89),
    BIC_mean = c(1000, 950, 1100, 1050),
    runtime_mean = c(10, 5, 12, 6),
    Kr_correct_rate = c(0.8, 0.7, 0.6, 0.75),
    stringsAsFactors = FALSE
  )
  
  win_rates <- compute_modelsel_win_rates(agg_df)
  
  expect_s3_class(win_rates, "data.frame")
  expect_true(nrow(win_rates) > 0)
  
  expect_true("metric" %in% names(win_rates))
  expect_true("method" %in% names(win_rates))
  expect_true("selection_strategy" %in% names(win_rates))
  expect_true("win_rate" %in% names(win_rates))
})

test_that("extract_selection_summary works for grid strategy", {
  skip_on_cran()
  
  data <- generate_synergy_data(
    N = 10, K = 2, r = 2, M = 6, T_each = 50, seed = 123
  )
  
  result <- select_K_r(
    list_of_data = data$list_of_data,
    method = "MFA",
    strategy = "grid",
    Kvec = 1:3,
    rvec = 1:3,
    max_iter = 10,
    n_init = 1,
    use_kmeans_init = FALSE,
    verbose = FALSE
  )
  
  summary <- extract_selection_summary(result)
  
  expect_s3_class(summary, "data.frame")
  expect_true("K" %in% names(summary))
  expect_true("r" %in% names(summary))
  expect_true("BIC" %in% names(summary))
})

test_that("extract_selection_summary works for twostage strategy", {
  skip_on_cran()
  
  data <- generate_synergy_data(
    N = 10, K = 2, r = 2, M = 6, T_each = 50, seed = 123
  )
  
  result <- select_K_r(
    list_of_data = data$list_of_data,
    method = "MFA",
    strategy = "twostage",
    Kvec = 1:3,
    rvec = 1:3,
    r_fixed = 2,
    max_iter = 10,
    n_init = 1,
    use_kmeans_init = FALSE,
    verbose = FALSE
  )
  
  summary <- extract_selection_summary(result)
  
  expect_s3_class(summary, "data.frame")
  expect_true("stage" %in% names(summary))
})
