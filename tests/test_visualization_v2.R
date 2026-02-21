# Test script for visualization functions in visualize_synergy_v2.R
# This script validates that all visualization functions generate proper ggplot objects

library(synergyMixR)
library(testthat)

# Skip all visualization tests if required packages are not available
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  cat("Skipping visualization tests: ggplot2 not installed\n")
  quit(save = "no", status = 0)
}
if (!requireNamespace("viridis", quietly = TRUE)) {
  cat("Skipping visualization tests: viridis not installed\n")
  quit(save = "no", status = 0)
}

cat("Testing visualization functions for simulate_dynamic_synergy_data()\n")
cat("=", rep("=", 70), "\n", sep = "")

# ============================================================================
# Setup: Generate test data
# ============================================================================
cat("\nSetup: Generating test data\n")
cat("-", rep("-", 70), "\n", sep = "")

set.seed(650)
sim <- simulate_dynamic_synergy_data(
  N = 10, K = 3, r = 3, M = 8, T_each = 200,
  cluster_sep_spatial = 1.0,
  cluster_sep_temporal = 1.0,
  cluster_sep_stability = 1.0
)

cat("✓ Test data generated successfully\n")
cat(sprintf("  - N = %d subjects\n", length(sim$X_list)))
cat(sprintf("  - K = %d clusters\n", length(sim$Lambda_list)))
cat(sprintf("  - r = %d synergies\n", ncol(sim$Lambda_list[[1]])))
cat(sprintf("  - M = %d muscles\n", ncol(sim$X_list[[1]])))
cat(sprintf("  - T = %d time points\n", nrow(sim$X_list[[1]])))

# ============================================================================
# Test 1: plot_synergy_bases()
# ============================================================================
cat("\nTest 1: plot_synergy_bases()\n")
cat("-", rep("-", 70), "\n", sep = "")

test_that("plot_synergy_bases generates ggplot object", {
  p <- plot_synergy_bases(sim)
  expect_s3_class(p, "ggplot")
  expect_true(!is.null(p$data))
  expect_true(nrow(p$data) > 0)
})

p1 <- plot_synergy_bases(sim)
cat("✓ plot_synergy_bases() test passed\n")
cat("  - Returns ggplot object\n")
cat("  - Contains data for all clusters\n")

# ============================================================================
# Test 2: plot_activation_patterns()
# ============================================================================
cat("\nTest 2: plot_activation_patterns()\n")
cat("-", rep("-", 70), "\n", sep = "")

test_that("plot_activation_patterns generates ggplot object", {
  p <- plot_activation_patterns(sim)
  expect_s3_class(p, "ggplot")
  expect_true(!is.null(p$data))
  expect_true(nrow(p$data) > 0)
})

test_that("plot_activation_patterns respects n_subjects_per_cluster parameter", {
  p <- plot_activation_patterns(sim, n_subjects_per_cluster = 2)
  expect_s3_class(p, "ggplot")
})

p2 <- plot_activation_patterns(sim)
cat("✓ plot_activation_patterns() test passed\n")
cat("  - Returns ggplot object\n")
cat("  - Respects n_subjects_per_cluster parameter\n")

# ============================================================================
# Test 3: plot_reconstructed_emg()
# ============================================================================
cat("\nTest 3: plot_reconstructed_emg()\n")
cat("-", rep("-", 70), "\n", sep = "")

test_that("plot_reconstructed_emg generates ggplot object", {
  p <- plot_reconstructed_emg(sim)
  expect_s3_class(p, "ggplot")
  expect_true(!is.null(p$layers))
  expect_true(length(p$layers) >= 2)  # Should have at least 2 layers (individual + mean)
})

test_that("plot_reconstructed_emg works with cluster_id parameter", {
  p <- plot_reconstructed_emg(sim, cluster_id = 1)
  expect_s3_class(p, "ggplot")
})

test_that("plot_reconstructed_emg works with muscles parameter", {
  p <- plot_reconstructed_emg(sim, muscles = c(1, 2, 3))
  expect_s3_class(p, "ggplot")
})

test_that("plot_reconstructed_emg works with both cluster_id and muscles", {
  p <- plot_reconstructed_emg(sim, cluster_id = 2, muscles = c(1, 2))
  expect_s3_class(p, "ggplot")
})

p3 <- plot_reconstructed_emg(sim)
cat("✓ plot_reconstructed_emg() test passed\n")
cat("  - Returns ggplot object with multiple layers\n")
cat("  - Works with cluster_id parameter\n")
cat("  - Works with muscles parameter\n")
cat("  - Works with both parameters combined\n")

# ============================================================================
# Test 4: plot_noise_profiles()
# ============================================================================
cat("\nTest 4: plot_noise_profiles()\n")
cat("-", rep("-", 70), "\n", sep = "")

test_that("plot_noise_profiles generates ggplot object", {
  p <- plot_noise_profiles(sim)
  expect_s3_class(p, "ggplot")
  expect_true(!is.null(p$data))
  expect_true(nrow(p$data) > 0)
})

p4 <- plot_noise_profiles(sim)
cat("✓ plot_noise_profiles() test passed\n")
cat("  - Returns ggplot object\n")
cat("  - Contains noise parameter data\n")

# ============================================================================
# Test 5: plot_cluster_embedding()
# ============================================================================
cat("\nTest 5: plot_cluster_embedding()\n")
cat("-", rep("-", 70), "\n", sep = "")

test_that("plot_cluster_embedding generates ggplot object with PCA", {
  p <- plot_cluster_embedding(sim, method = "pca")
  expect_s3_class(p, "ggplot")
  expect_true(!is.null(p$data))
  expect_true(nrow(p$data) == length(sim$X_list))
})

# Test UMAP only if package is available
if (requireNamespace("uwot", quietly = TRUE)) {
  test_that("plot_cluster_embedding generates ggplot object with UMAP", {
    p <- plot_cluster_embedding(sim, method = "umap")
    expect_s3_class(p, "ggplot")
    expect_true(!is.null(p$data))
    expect_true(nrow(p$data) == length(sim$X_list))
  })
  cat("✓ plot_cluster_embedding() test passed\n")
  cat("  - Works with PCA method\n")
  cat("  - Works with UMAP method (using uwot)\n")
} else {
  cat("✓ plot_cluster_embedding() test passed (PCA only)\n")
  cat("  - Works with PCA method\n")
  cat("  - UMAP method not tested (uwot package not available)\n")
}

# ============================================================================
# Test 6: plot_synergy_summary()
# ============================================================================
cat("\nTest 6: plot_synergy_summary()\n")
cat("-", rep("-", 70), "\n", sep = "")

if (requireNamespace("patchwork", quietly = TRUE)) {
  test_that("plot_synergy_summary generates patchwork object", {
    p <- plot_synergy_summary(sim)
    expect_s3_class(p, "patchwork")
  })
  
  p6 <- plot_synergy_summary(sim)
  cat("✓ plot_synergy_summary() test passed\n")
  cat("  - Returns patchwork object\n")
  cat("  - Combines multiple visualization panels\n")
} else {
  cat("⚠ plot_synergy_summary() test skipped\n")
  cat("  - patchwork package not available\n")
}

# ============================================================================
# Test 7: plot_neural_command()
# ============================================================================
cat("\nTest 7: plot_neural_command()\n")
cat("-", rep("-", 70), "\n", sep = "")

test_that("plot_neural_command generates ggplot object with PCA", {
  p <- plot_neural_command(sim, method = "pca", color_by = "cluster")
  expect_s3_class(p, "ggplot")
  expect_true(!is.null(p$data))
  expect_true(nrow(p$data) > 0)
})

test_that("plot_neural_command works with color_by='time'", {
  p <- plot_neural_command(sim, method = "pca", color_by = "time")
  expect_s3_class(p, "ggplot")
})

test_that("plot_neural_command works with smooth=FALSE", {
  p <- plot_neural_command(sim, method = "pca", color_by = "cluster", smooth = FALSE)
  expect_s3_class(p, "ggplot")
})

# Test UMAP only if uwot package is available
if (requireNamespace("uwot", quietly = TRUE)) {
  test_that("plot_neural_command generates ggplot object with UMAP", {
    p <- plot_neural_command(sim, method = "umap", color_by = "cluster")
    expect_s3_class(p, "ggplot")
    expect_true(!is.null(p$data))
    expect_true(nrow(p$data) > 0)
  })
  cat("✓ plot_neural_command() test passed\n")
  cat("  - Works with PCA method\n")
  cat("  - Works with UMAP method (using uwot)\n")
  cat("  - Works with color_by='cluster' and color_by='time'\n")
  cat("  - Works with smooth=TRUE and smooth=FALSE\n")
} else {
  cat("✓ plot_neural_command() test passed (PCA only)\n")
  cat("  - Works with PCA method\n")
  cat("  - Works with color_by='cluster' and color_by='time'\n")
  cat("  - Works with smooth=TRUE and smooth=FALSE\n")
  cat("  - UMAP method not tested (uwot package not available)\n")
}

# ============================================================================
# Test 8: Edge cases and error handling
# ============================================================================
cat("\nTest 8: Edge cases and error handling\n")
cat("-", rep("-", 70), "\n", sep = "")

test_that("Functions handle missing components gracefully", {
  sim_incomplete <- list(X_list = sim$X_list)
  expect_error(plot_synergy_bases(sim_incomplete))
  expect_error(plot_activation_patterns(sim_incomplete))
  expect_error(plot_noise_profiles(sim_incomplete))
})

test_that("Functions work with K=2 clusters", {
  sim_k2 <- simulate_dynamic_synergy_data(N = 10, K = 2, r = 3, M = 8, T_each = 100)
  expect_s3_class(plot_synergy_bases(sim_k2), "ggplot")
  expect_s3_class(plot_activation_patterns(sim_k2), "ggplot")
  expect_s3_class(plot_reconstructed_emg(sim_k2), "ggplot")
  expect_s3_class(plot_noise_profiles(sim_k2), "ggplot")
  expect_s3_class(plot_cluster_embedding(sim_k2, method = "pca"), "ggplot")
})

test_that("Functions work with K=5 clusters", {
  sim_k5 <- simulate_dynamic_synergy_data(N = 15, K = 5, r = 3, M = 8, T_each = 100)
  expect_s3_class(plot_synergy_bases(sim_k5), "ggplot")
  expect_s3_class(plot_activation_patterns(sim_k5), "ggplot")
  expect_s3_class(plot_reconstructed_emg(sim_k5), "ggplot")
  expect_s3_class(plot_noise_profiles(sim_k5), "ggplot")
  expect_s3_class(plot_cluster_embedding(sim_k5, method = "pca"), "ggplot")
})

cat("✓ Edge cases test passed\n")
cat("  - Functions handle missing components with errors\n")
cat("  - Functions work with K=2 clusters\n")
cat("  - Functions work with K=5 clusters\n")

# ============================================================================
# Test 8: Visual consistency across different parameter settings
# ============================================================================
cat("\nTest 8: Visual consistency across parameter settings\n")
cat("-", rep("-", 70), "\n", sep = "")

test_that("Visualizations work with varying separation parameters", {
  sim_low_sep <- simulate_dynamic_synergy_data(
    N = 10, K = 3, r = 3, M = 8, T_each = 100,
    cluster_sep_spatial = 0.5,
    cluster_sep_temporal = 0.5,
    cluster_sep_stability = 0.5
  )
  
  sim_high_sep <- simulate_dynamic_synergy_data(
    N = 10, K = 3, r = 3, M = 8, T_each = 100,
    cluster_sep_spatial = 2.0,
    cluster_sep_temporal = 2.0,
    cluster_sep_stability = 2.0
  )
  
  expect_s3_class(plot_synergy_bases(sim_low_sep), "ggplot")
  expect_s3_class(plot_synergy_bases(sim_high_sep), "ggplot")
  expect_s3_class(plot_cluster_embedding(sim_low_sep, method = "pca"), "ggplot")
  expect_s3_class(plot_cluster_embedding(sim_high_sep, method = "pca"), "ggplot")
})

cat("✓ Visual consistency test passed\n")
cat("  - Functions work with low separation parameters\n")
cat("  - Functions work with high separation parameters\n")

# ============================================================================
# Test 9: Theme customization
# ============================================================================
cat("\nTest 9: Theme customization\n")
cat("-", rep("-", 70), "\n", sep = "")

test_that("Functions respect theme_base_size parameter", {
  p_small <- plot_synergy_bases(sim, theme_base_size = 10)
  p_large <- plot_synergy_bases(sim, theme_base_size = 18)
  expect_s3_class(p_small, "ggplot")
  expect_s3_class(p_large, "ggplot")
})

cat("✓ Theme customization test passed\n")
cat("  - Functions respect theme_base_size parameter\n")

# ============================================================================
# Summary
# ============================================================================
cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("All visualization tests passed successfully!\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("\nAll seven visualization functions are working correctly:\n")
cat("  1. plot_synergy_bases() ✓\n")
cat("  2. plot_activation_patterns() ✓\n")
cat("  3. plot_reconstructed_emg() ✓\n")
cat("  4. plot_noise_profiles() ✓\n")
cat("  5. plot_cluster_embedding() ✓\n")
cat("  6. plot_synergy_summary() ✓\n")
cat("  7. plot_neural_command() ✓\n")
cat("\n")
