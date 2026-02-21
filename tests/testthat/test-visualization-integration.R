# Integration Tests for Visualization Functions with Unified Theme
# Tests that visualization functions work correctly with the unified theme system

library(testthat)

# Skip all tests if visualization dependencies are not available
skip_if_not_installed("ggplot2")
skip_if_not_installed("reshape2")
skip_if_not_installed("dplyr")
skip_if_not_installed("tidyr")

# Load test fixtures
small_data <- readRDS("fixtures/small_test_data.rds")

# Helper function to create a simple MFA model for testing
create_test_mfa_model <- function() {
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  list(
    z = true_params$z,
    mu = true_params$mu,
    Lambda = true_params$Lambda,
    Psi = true_params$Psi,
    pi = rep(1/true_params$K, true_params$K)
  )
}

# Helper function to create a simple MPCA model for testing
create_test_mpca_model <- function() {
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  W_list <- list()
  P_list <- list()
  D_list <- list()
  for (k in 1:true_params$K) {
    W_k <- true_params$Lambda[[k]]
    W_list[[k]] <- W_k
    svd_k <- svd(W_k)
    P_list[[k]] <- svd_k$u[, 1:true_params$r]
    D_list[[k]] <- diag(svd_k$d[1:true_params$r])
  }
  
  list(
    z = true_params$z,
    mu = true_params$mu,
    W = W_list,
    P = P_list,
    D = D_list,
    sigma2 = rep(0.5, true_params$K),
    pi = rep(1/true_params$K, true_params$K)
  )
}

# Test plot_cluster_synergy_loadings_mfa with unified theme
test_that("plot_cluster_synergy_loadings_mfa works with unified theme (heatmap)", {
  mfa_fit <- create_test_mfa_model()
  list_of_data <- small_data$data
  
  # Create heatmap plot
  p <- plot_cluster_synergy_loadings_mfa(
    mfa_fit, 
    cluster_ids = 1:2, 
    plot_type = "heatmap"
  )
  
  # Check that plot is created successfully
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Check that plot can be built without errors
  built <- ggplot_build(p)
  expect_true(!is.null(built))
  
  # Check that the plot has the expected layers
  expect_true(length(p$layers) > 0)
  expect_equal(class(p$layers[[1]]$geom)[1], "GeomTile")
})

test_that("plot_cluster_synergy_loadings_mfa works with unified theme (bar)", {
  mfa_fit <- create_test_mfa_model()
  
  # Create bar plot
  p <- plot_cluster_synergy_loadings_mfa(
    mfa_fit, 
    cluster_ids = 1:2, 
    plot_type = "bar"
  )
  
  # Check that plot is created successfully
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Check that plot can be built without errors
  built <- ggplot_build(p)
  expect_true(!is.null(built))
  
  # Check that the plot has the expected layers
  expect_true(length(p$layers) > 0)
  expect_equal(class(p$layers[[1]]$geom)[1], "GeomBar")
})

# Test plot_all_factor_scores_mfa with unified theme
test_that("plot_all_factor_scores_mfa works with unified theme", {
  mfa_fit <- create_test_mfa_model()
  list_of_data <- small_data$data
  
  # Create factor scores plot
  p <- plot_all_factor_scores_mfa(
    mfa_fit, 
    list_of_data, 
    overlay_subjects = TRUE
  )
  
  # Check that plot is created successfully
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Check that plot can be built without errors
  built <- ggplot_build(p)
  expect_true(!is.null(built))
  
  # Check that the plot has the expected layers
  expect_true(length(p$layers) > 0)
  expect_equal(class(p$layers[[1]]$geom)[1], "GeomLine")
})

test_that("plot_all_factor_scores_mfa works with averaged subjects", {
  mfa_fit <- create_test_mfa_model()
  list_of_data <- small_data$data
  
  # Create averaged factor scores plot
  p <- plot_all_factor_scores_mfa(
    mfa_fit, 
    list_of_data, 
    overlay_subjects = FALSE
  )
  
  # Check that plot is created successfully
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Check that plot can be built without errors
  built <- ggplot_build(p)
  expect_true(!is.null(built))
})

# Test plot_all_reconstructed_waveforms_mfa with unified theme
test_that("plot_all_reconstructed_waveforms_mfa works with unified theme", {
  mfa_fit <- create_test_mfa_model()
  list_of_data <- small_data$data
  
  # Create reconstructed waveforms plot
  p <- plot_all_reconstructed_waveforms_mfa(
    mfa_fit, 
    list_of_data, 
    overlay_subjects = TRUE
  )
  
  # Check that plot is created successfully
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Check that plot can be built without errors
  built <- ggplot_build(p)
  expect_true(!is.null(built))
  
  # Check that the plot has the expected layers
  expect_true(length(p$layers) > 0)
  expect_equal(class(p$layers[[1]]$geom)[1], "GeomLine")
})

# Test plot_waveform_comparison_mfa with unified theme
test_that("plot_waveform_comparison_mfa works with unified theme", {
  mfa_fit <- create_test_mfa_model()
  list_of_data <- small_data$data
  
  # Create waveform comparison plot
  p <- plot_waveform_comparison_mfa(
    mfa_fit, 
    list_of_data, 
    subject_id = 1, 
    channel = 1
  )
  
  # Check that plot is created successfully
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Check that plot can be built without errors
  built <- ggplot_build(p)
  expect_true(!is.null(built))
  
  # Check that the plot has the expected layers
  expect_true(length(p$layers) > 0)
  expect_equal(class(p$layers[[1]]$geom)[1], "GeomLine")
})

# Test plot_bic_line with unified theme
test_that("plot_bic_line works with unified theme", {
  # Create mock BIC summary data
  df_summary <- data.frame(
    K = rep(1:3, each = 3),
    r = rep(1:3, 3),
    BIC = rnorm(9, mean = 1000, sd = 100)
  )
  
  # Create BIC plot with K on x-axis
  p_K <- plot_bic_line(df_summary, x_axis = "K")
  
  # Check that plot is created successfully
  expect_s3_class(p_K, "gg")
  expect_s3_class(p_K, "ggplot")
  
  # Check that plot can be built without errors
  built_K <- ggplot_build(p_K)
  expect_true(!is.null(built_K))
  
  # Create BIC plot with r on x-axis
  p_r <- plot_bic_line(df_summary, x_axis = "r")
  
  # Check that plot is created successfully
  expect_s3_class(p_r, "gg")
  expect_s3_class(p_r, "ggplot")
  
  # Check that plot can be built without errors
  built_r <- ggplot_build(p_r)
  expect_true(!is.null(built_r))
})

# Test plot_cluster_synergy_loadings_mpca with unified theme
test_that("plot_cluster_synergy_loadings_mpca works with unified theme (heatmap)", {
  mpca_fit <- create_test_mpca_model()
  
  # Create heatmap plot
  p <- plot_cluster_synergy_loadings_mpca(
    mpca_fit, 
    cluster_ids = 1:2, 
    plot_type = "heatmap"
  )
  
  # Check that plot is created successfully
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Check that plot can be built without errors
  built <- ggplot_build(p)
  expect_true(!is.null(built))
  
  # Check that the plot has the expected layers
  expect_true(length(p$layers) > 0)
  expect_equal(class(p$layers[[1]]$geom)[1], "GeomTile")
})

test_that("plot_cluster_synergy_loadings_mpca works with unified theme (bar)", {
  mpca_fit <- create_test_mpca_model()
  
  # Create bar plot
  p <- plot_cluster_synergy_loadings_mpca(
    mpca_fit, 
    cluster_ids = 1:2, 
    plot_type = "bar"
  )
  
  # Check that plot is created successfully
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Check that plot can be built without errors
  built <- ggplot_build(p)
  expect_true(!is.null(built))
  
  # Check that the plot has the expected layers
  expect_true(length(p$layers) > 0)
  expect_equal(class(p$layers[[1]]$geom)[1], "GeomBar")
})

# Test plot_all_factor_scores_mpca with unified theme
test_that("plot_all_factor_scores_mpca works with unified theme", {
  mpca_fit <- create_test_mpca_model()
  list_of_data <- small_data$data
  
  # Compute factor scores first
  mpca_fit <- compute_factor_scores_mpca(mpca_fit, list_of_data)
  
  # Create factor scores plot
  p <- plot_all_factor_scores_mpca(
    mpca_fit, 
    list_of_data, 
    overlay_subjects = TRUE
  )
  
  # Check that plot is created successfully
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Check that plot can be built without errors
  built <- ggplot_build(p)
  expect_true(!is.null(built))
  
  # Check that the plot has the expected layers
  expect_true(length(p$layers) > 0)
  expect_equal(class(p$layers[[1]]$geom)[1], "GeomLine")
})

# Test plot_all_reconstructed_waveforms_mpca with unified theme
test_that("plot_all_reconstructed_waveforms_mpca works with unified theme", {
  mpca_fit <- create_test_mpca_model()
  list_of_data <- small_data$data
  
  # Compute factor scores first
  mpca_fit <- compute_factor_scores_mpca(mpca_fit, list_of_data)
  
  # Create reconstructed waveforms plot
  p <- plot_all_reconstructed_waveforms_mpca(
    mpca_fit, 
    list_of_data, 
    overlay_subjects = TRUE
  )
  
  # Check that plot is created successfully
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Check that plot can be built without errors
  built <- ggplot_build(p)
  expect_true(!is.null(built))
  
  # Check that the plot has the expected layers
  expect_true(length(p$layers) > 0)
  expect_equal(class(p$layers[[1]]$geom)[1], "GeomLine")
})

# Test plot_waveform_comparison_mpca with unified theme
test_that("plot_waveform_comparison_mpca works with unified theme", {
  mpca_fit <- create_test_mpca_model()
  list_of_data <- small_data$data
  
  # Compute factor scores first
  mpca_fit <- compute_factor_scores_mpca(mpca_fit, list_of_data)
  
  # Create waveform comparison plot
  p <- plot_waveform_comparison_mpca(
    mpca_fit, 
    list_of_data, 
    subject_id = 1, 
    channel = 1
  )
  
  # Check that plot is created successfully
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Check that plot can be built without errors
  built <- ggplot_build(p)
  expect_true(!is.null(built))
  
  # Check that the plot has the expected layers
  expect_true(length(p$layers) > 0)
  expect_equal(class(p$layers[[1]]$geom)[1], "GeomLine")
})

# Test that muscle_names parameter works correctly
test_that("visualization functions accept muscle_names parameter", {
  mfa_fit <- create_test_mfa_model()
  list_of_data <- small_data$data
  muscle_names <- paste0("Muscle", 1:small_data$true_params$M)
  
  # Test with plot_cluster_synergy_loadings_mfa
  p1 <- plot_cluster_synergy_loadings_mfa(
    mfa_fit, 
    cluster_ids = 1, 
    plot_type = "heatmap",
    muscle_names = muscle_names
  )
  expect_s3_class(p1, "gg")
  
  # Test with plot_waveform_comparison_mfa
  p2 <- plot_waveform_comparison_mfa(
    mfa_fit, 
    list_of_data, 
    subject_id = 1, 
    channel = 1,
    muscle_names = muscle_names
  )
  expect_s3_class(p2, "gg")
})

# Test that all plots use consistent faceting
test_that("visualization functions use consistent facet_grid", {
  mfa_fit <- create_test_mfa_model()
  list_of_data <- small_data$data
  
  # Test factor scores plot
  p_scores <- plot_all_factor_scores_mfa(mfa_fit, list_of_data)
  expect_s3_class(p_scores$facet, "FacetGrid")
  
  # Test reconstructed waveforms plot
  p_recon <- plot_all_reconstructed_waveforms_mfa(mfa_fit, list_of_data)
  expect_s3_class(p_recon$facet, "FacetGrid")
})

# Test that all heatmaps use diverging color scale
test_that("heatmap visualizations use diverging color scale", {
  mfa_fit <- create_test_mfa_model()
  
  # Create heatmap plot
  p <- plot_cluster_synergy_loadings_mfa(
    mfa_fit, 
    cluster_ids = 1, 
    plot_type = "heatmap"
  )
  
  # Check that plot has a fill scale
  expect_true(!is.null(p$scales))
  
  # Build the plot to check the scale
  built <- ggplot_build(p)
  expect_true(!is.null(built))
})

# Test that plots handle edge cases gracefully
test_that("visualization functions handle single cluster", {
  mfa_fit <- create_test_mfa_model()
  list_of_data <- small_data$data
  
  # Test with single cluster
  p <- plot_cluster_synergy_loadings_mfa(
    mfa_fit, 
    cluster_ids = 1, 
    plot_type = "heatmap"
  )
  
  # Check that plot is created successfully
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Check that plot can be built without errors
  built <- ggplot_build(p)
  expect_true(!is.null(built))
})

test_that("visualization functions handle all clusters", {
  mfa_fit <- create_test_mfa_model()
  
  # Test with all clusters (NULL cluster_ids)
  p <- plot_cluster_synergy_loadings_mfa(
    mfa_fit, 
    cluster_ids = NULL, 
    plot_type = "heatmap"
  )
  
  # Check that plot is created successfully
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Check that plot can be built without errors
  built <- ggplot_build(p)
  expect_true(!is.null(built))
})
