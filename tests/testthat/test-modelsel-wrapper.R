# Tests for model selection wrappers
# test-modelsel-wrapper.R

test_that("generate_synergy_data returns expected structure", {
  skip_on_cran()
  
  data <- generate_synergy_data(
    N = 10, K = 2, r = 2, M = 6, T_each = 50,
    cluster_sep_spatial = 1.0,
    cluster_sep_temporal = 1.0,
    cluster_sep_stability = 1.0,
    seed = 123
  )
  
  expect_type(data, "list")
  expect_true("list_of_data" %in% names(data))
  expect_true("true_cluster" %in% names(data))
  expect_true("params" %in% names(data))
  
  expect_length(data$list_of_data, 10)
  expect_length(data$true_cluster, 10)
  
  expect_true(all(data$true_cluster %in% 1:2))
  
  expect_equal(nrow(data$list_of_data[[1]]), 50)
  expect_equal(ncol(data$list_of_data[[1]]), 6)
  
  expect_equal(data$params$N, 10)
  expect_equal(data$params$K, 2)
  expect_equal(data$params$r, 2)
  expect_equal(data$params$seed, 123)
})

test_that("generate_synergy_data requires seed", {
  expect_error(
    generate_synergy_data(N = 10, K = 2, r = 2),
    "seed is required"
  )
})

test_that("generate_synergy_data is reproducible with same seed", {
  skip_on_cran()
  
  data1 <- generate_synergy_data(
    N = 5, K = 2, r = 2, M = 6, T_each = 30, seed = 456
  )
  data2 <- generate_synergy_data(
    N = 5, K = 2, r = 2, M = 6, T_each = 30, seed = 456
  )
  
  expect_equal(data1$true_cluster, data2$true_cluster)
  expect_equal(data1$list_of_data[[1]], data2$list_of_data[[1]])
})

test_that("select_K_r returns expected structure for MFA grid", {
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
  
  expect_type(result, "list")
  expect_true("best_K" %in% names(result))
  expect_true("best_r" %in% names(result))
  expect_true("best_model" %in% names(result))
  expect_true("summary_df" %in% names(result))
  expect_true("meta" %in% names(result))
  
  expect_true(is.numeric(result$best_K))
  expect_true(is.numeric(result$best_r))
  expect_true(result$best_K >= 1 && result$best_K <= 3)
  expect_true(result$best_r >= 1 && result$best_r <= 3)
  
  expect_equal(result$meta$method, "MFA")
  expect_equal(result$meta$strategy, "grid")
})

test_that("select_K_r returns expected structure for MPCA grid", {
  skip_on_cran()
  
  data <- generate_synergy_data(
    N = 10, K = 2, r = 2, M = 6, T_each = 50, seed = 123
  )
  
  result <- select_K_r(
    list_of_data = data$list_of_data,
    method = "MPCA",
    strategy = "grid",
    Kvec = 1:3,
    rvec = 1:3,
    max_iter = 10,
    n_init = 1,
    use_kmeans_init = FALSE,
    verbose = FALSE
  )
  
  expect_type(result, "list")
  expect_true("best_K" %in% names(result))
  expect_true("best_r" %in% names(result))
  expect_true("best_model" %in% names(result))
  
  expect_equal(result$meta$method, "MPCA")
  expect_equal(result$meta$strategy, "grid")
})

test_that("select_K_r returns expected structure for twostage", {
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
  
  expect_type(result, "list")
  expect_true("best_K" %in% names(result))
  expect_true("best_r" %in% names(result))
  expect_true("best_model" %in% names(result))
  
  expect_equal(result$meta$method, "MFA")
  expect_equal(result$meta$strategy, "twostage")
  
  expect_type(result$summary_df, "list")
  expect_true("stage1" %in% names(result$summary_df))
  expect_true("stage2" %in% names(result$summary_df))
})

test_that("extract_cluster_labels works for MFA fit", {
  skip_on_cran()
  
  data <- generate_synergy_data(
    N = 10, K = 2, r = 2, M = 6, T_each = 50, seed = 123
  )
  
  fit <- mfa_em_fit(
    data$list_of_data, K = 2, r = 2,
    max_iter = 10, n_init = 1, use_kmeans_init = FALSE, verbose = FALSE
  )
  
  labels <- extract_cluster_labels(fit)
  
  expect_length(labels, 10)
  expect_true(all(labels %in% 1:2))
})

test_that("extract_cluster_labels works for MPCA fit", {
  skip_on_cran()
  
  data <- generate_synergy_data(
    N = 10, K = 2, r = 2, M = 6, T_each = 50, seed = 123
  )
  
  fit <- mixture_pca_em_fit(
    data$list_of_data, K = 2, r = 2,
    max_iter = 10, n_init = 1, use_kmeans_init = FALSE, verbose = FALSE
  )
  
  labels <- extract_cluster_labels(fit)
  
  expect_length(labels, 10)
  expect_true(all(labels %in% 1:2))
})

test_that("extract_cluster_labels returns NA for NULL input", {
  expect_true(is.na(extract_cluster_labels(NULL)))
})

test_that("compute_ari returns valid values", {
  skip_if_not_installed("mclust")
  true_labels <- c(1, 1, 1, 2, 2, 2)
  pred_labels <- c(1, 1, 2, 2, 2, 2)

  ari <- compute_ari(true_labels, pred_labels)

  expect_true(is.numeric(ari))
  expect_true(ari >= -1 && ari <= 1)
})

test_that("compute_ari returns 1 for perfect match", {
  skip_if_not_installed("mclust")
  labels <- c(1, 1, 1, 2, 2, 2)

  ari <- compute_ari(labels, labels)

  expect_equal(ari, 1)
})

test_that("compute_ari handles NA inputs", {
  expect_true(is.na(compute_ari(c(1, NA, 2), c(1, 2, 2))))
  expect_true(is.na(compute_ari(c(1, 2, 2), c(1, NA, 2))))
})

test_that("compute_ari handles length mismatch", {
  expect_warning(
    result <- compute_ari(c(1, 2), c(1, 2, 3)),
    "Length mismatch"
  )
  expect_true(is.na(result))
})

test_that("compute_vaf_wrapper works for MFA", {
  skip_on_cran()
  
  data <- generate_synergy_data(
    N = 10, K = 2, r = 2, M = 6, T_each = 50, seed = 123
  )
  
  fit <- mfa_em_fit(
    data$list_of_data, K = 2, r = 2,
    max_iter = 10, n_init = 1, use_kmeans_init = FALSE, verbose = FALSE
  )
  
  vaf <- compute_vaf_wrapper("MFA", data$list_of_data, fit)
  
  expect_true(is.numeric(vaf))
  expect_true(vaf >= 0 && vaf <= 1)
})

test_that("compute_vaf_wrapper works for MPCA", {
  skip_on_cran()
  
  data <- generate_synergy_data(
    N = 10, K = 2, r = 2, M = 6, T_each = 50, seed = 123
  )
  
  fit <- mixture_pca_em_fit(
    data$list_of_data, K = 2, r = 2,
    max_iter = 10, n_init = 1, use_kmeans_init = FALSE, verbose = FALSE
  )
  
  vaf <- compute_vaf_wrapper("MPCA", data$list_of_data, fit)
  
  expect_true(is.numeric(vaf))
  expect_true(vaf >= 0 && vaf <= 1)
})

test_that("compute_vaf_wrapper returns NA for NULL fit", {
  expect_true(is.na(compute_vaf_wrapper("MFA", list(), NULL)))
})

test_that("compute_bic_wrapper works for MFA", {
  skip_on_cran()
  
  data <- generate_synergy_data(
    N = 10, K = 2, r = 2, M = 6, T_each = 50, seed = 123
  )
  
  fit <- mfa_em_fit(
    data$list_of_data, K = 2, r = 2,
    max_iter = 10, n_init = 1, use_kmeans_init = FALSE, verbose = FALSE
  )
  
  N_obs <- sum(sapply(data$list_of_data, nrow))
  bic <- compute_bic_wrapper("MFA", fit, K = 2, r = 2, M = 6, N_obs = N_obs)
  
  expect_true(is.numeric(bic))
  expect_true(is.finite(bic))
})

test_that("evaluate_model_selection returns complete results", {
  skip_on_cran()
  
  data <- generate_synergy_data(
    N = 10, K = 2, r = 2, M = 6, T_each = 50, seed = 123
  )
  
  select_result <- select_K_r(
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
  
  eval_result <- evaluate_model_selection(
    list_of_data = data$list_of_data,
    true_cluster = data$true_cluster,
    K_true = 2,
    r_true = 2,
    select_result = select_result,
    runtime_sec = 1.5
  )
  
  expect_type(eval_result, "list")
  expect_true("K_hat" %in% names(eval_result))
  expect_true("r_hat" %in% names(eval_result))
  expect_true("K_correct" %in% names(eval_result))
  expect_true("r_correct" %in% names(eval_result))
  expect_true("Kr_correct" %in% names(eval_result))
  expect_true("ARI" %in% names(eval_result))
  expect_true("VAF" %in% names(eval_result))
  expect_true("logLik" %in% names(eval_result))
  expect_true("BIC" %in% names(eval_result))
  expect_true("runtime_sec" %in% names(eval_result))
  
  expect_true(is.logical(eval_result$K_correct))
  expect_true(is.logical(eval_result$r_correct))
  expect_true(is.logical(eval_result$Kr_correct))
  expect_equal(eval_result$runtime_sec, 1.5)
})
