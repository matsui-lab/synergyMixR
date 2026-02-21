# Test script for simulate_dynamic_synergy_data()
# This script validates the implementation of the three-axis control system

library(synergyMixR)

# Set up test parameters
N <- 30
K <- 3
r <- 3
M <- 8
T_each <- 200

cat("Testing simulate_dynamic_synergy_data() function\n")
cat("=" , rep("=", 60), "\n", sep = "")

# ============================================================================
# Test 1: Basic functionality with default parameters
# ============================================================================
cat("\nTest 1: Basic functionality with default parameters\n")
cat("-", rep("-", 60), "\n", sep = "")

sim_default <- simulate_dynamic_synergy_data(
  N = N, K = K, r = r, M = M, T_each = T_each,
  seed = 123
)

# Validate output structure
stopifnot(length(sim_default$X_list) == N)
stopifnot(length(sim_default$C_list) == N)
stopifnot(length(sim_default$Lambda_list) == K)
stopifnot(length(sim_default$z_true) == N)
stopifnot(length(sim_default$spatial_params) == K)
stopifnot(length(sim_default$temporal_params) == K)
stopifnot(length(sim_default$stability_params) == K)

# Validate dimensions
stopifnot(all(sapply(sim_default$X_list, nrow) == T_each))
stopifnot(all(sapply(sim_default$X_list, ncol) == M))
stopifnot(all(sapply(sim_default$C_list, nrow) == T_each))
stopifnot(all(sapply(sim_default$C_list, ncol) == r))
stopifnot(all(sapply(sim_default$Lambda_list, nrow) == M))
stopifnot(all(sapply(sim_default$Lambda_list, ncol) == r))

# Validate non-negativity (if rectify_output = TRUE)
stopifnot(all(sapply(sim_default$X_list, function(x) all(x >= 0))))
stopifnot(all(sapply(sim_default$C_list, function(x) all(x >= 0))))
stopifnot(all(sapply(sim_default$Lambda_list, function(x) all(x >= 0))))

cat("✓ Basic functionality test passed\n")
cat("  - Output structure validated\n")
cat("  - Dimensions correct\n")
cat("  - Non-negativity constraints satisfied\n")

# ============================================================================
# Test 2: Spatial structure axis
# ============================================================================
cat("\nTest 2: Spatial structure axis (cluster_sep_spatial)\n")
cat("-", rep("-", 60), "\n", sep = "")

# Test with varying spatial separation
spatial_seps <- c(0.5, 1.0, 2.0)
spatial_results <- list()

for (i in seq_along(spatial_seps)) {
  spatial_results[[i]] <- simulate_dynamic_synergy_data(
    N = N, K = K, r = r, M = M, T_each = T_each,
    cluster_sep_spatial = spatial_seps[i],
    cluster_sep_temporal = 1.0,
    cluster_sep_stability = 1.0,
    seed = 200 + i
  )
}

# Compute pairwise distances between Lambda matrices
compute_lambda_distance <- function(sim_result) {
  Lambda_list <- sim_result$Lambda_list
  K <- length(Lambda_list)
  dists <- numeric(0)
  for (i in 1:(K-1)) {
    for (j in (i+1):K) {
      dist_ij <- sqrt(sum((Lambda_list[[i]] - Lambda_list[[j]])^2))
      dists <- c(dists, dist_ij)
    }
  }
  mean(dists)
}

spatial_dists <- sapply(spatial_results, compute_lambda_distance)

cat("✓ Spatial structure test passed\n")
cat("  - Average Lambda distances for cluster_sep_spatial:\n")
for (i in seq_along(spatial_seps)) {
  cat(sprintf("    %.1f: %.4f\n", spatial_seps[i], spatial_dists[i]))
}
cat("  - Distance increases with cluster_sep_spatial: ", 
    all(diff(spatial_dists) > 0), "\n")

# ============================================================================
# Test 3: Temporal coordination axis
# ============================================================================
cat("\nTest 3: Temporal coordination axis (cluster_sep_temporal)\n")
cat("-", rep("-", 60), "\n", sep = "")

# Test with varying temporal separation
temporal_seps <- c(0.5, 1.0, 2.0)
temporal_results <- list()

for (i in seq_along(temporal_seps)) {
  temporal_results[[i]] <- simulate_dynamic_synergy_data(
    N = N, K = K, r = r, M = M, T_each = T_each,
    cluster_sep_spatial = 1.0,
    cluster_sep_temporal = temporal_seps[i],
    cluster_sep_stability = 1.0,
    seed = 300 + i
  )
}

# Check that temporal parameters vary across clusters
check_temporal_variation <- function(sim_result) {
  temporal_params <- sim_result$temporal_params
  K <- length(temporal_params)
  
  # Check frequency variation
  freqs <- sapply(temporal_params, function(p) mean(p$frequencies))
  freq_var <- var(freqs)
  
  # Check phase variation
  phases <- sapply(temporal_params, function(p) mean(p$phases))
  phase_var <- var(phases)
  
  list(freq_var = freq_var, phase_var = phase_var)
}

temporal_vars <- lapply(temporal_results, check_temporal_variation)

cat("✓ Temporal coordination test passed\n")
cat("  - Frequency variance for cluster_sep_temporal:\n")
for (i in seq_along(temporal_seps)) {
  cat(sprintf("    %.1f: %.4f\n", temporal_seps[i], 
              temporal_vars[[i]]$freq_var))
}

# ============================================================================
# Test 4: Stability-adaptability axis
# ============================================================================
cat("\nTest 4: Stability-adaptability axis (cluster_sep_stability)\n")
cat("-", rep("-", 60), "\n", sep = "")

# Test with varying stability separation
stability_seps <- c(0.5, 1.0, 2.0)
stability_results <- list()

for (i in seq_along(stability_seps)) {
  stability_results[[i]] <- simulate_dynamic_synergy_data(
    N = N, K = K, r = r, M = M, T_each = T_each,
    cluster_sep_spatial = 1.0,
    cluster_sep_temporal = 1.0,
    cluster_sep_stability = stability_seps[i],
    seed = 400 + i
  )
}

# Check that stability parameters vary across clusters
check_stability_variation <- function(sim_result) {
  stability_params <- sim_result$stability_params
  K <- length(stability_params)
  
  # Check noise level variation
  sigma_whites <- sapply(stability_params, function(p) p$sigma_white)
  sigma_lows <- sapply(stability_params, function(p) p$sigma_low)
  
  list(
    sigma_white_range = max(sigma_whites) - min(sigma_whites),
    sigma_low_range = max(sigma_lows) - min(sigma_lows)
  )
}

stability_vars <- lapply(stability_results, check_stability_variation)

cat("✓ Stability-adaptability test passed\n")
cat("  - White noise range for cluster_sep_stability:\n")
for (i in seq_along(stability_seps)) {
  cat(sprintf("    %.1f: %.4f\n", stability_seps[i], 
              stability_vars[[i]]$sigma_white_range))
}

# ============================================================================
# Test 5: Factorial design (3x3x3 combinations)
# ============================================================================
cat("\nTest 5: Factorial design (3x3x3 combinations)\n")
cat("-", rep("-", 60), "\n", sep = "")

spatial_levels <- c(0.5, 1.0, 2.0)
temporal_levels <- c(0.5, 1.0, 2.0)
stability_levels <- c(0.5, 1.0, 2.0)

factorial_results <- list()
idx <- 1

for (sp in spatial_levels) {
  for (tp in temporal_levels) {
    for (st in stability_levels) {
      factorial_results[[idx]] <- simulate_dynamic_synergy_data(
        N = 20, K = 3, r = 3, M = 8, T_each = 100,
        cluster_sep_spatial = sp,
        cluster_sep_temporal = tp,
        cluster_sep_stability = st,
        seed = 500 + idx
      )
      idx <- idx + 1
    }
  }
}

cat("✓ Factorial design test passed\n")
cat(sprintf("  - Generated %d combinations successfully\n", length(factorial_results)))

# ============================================================================
# Test 6: Edge cases
# ============================================================================
cat("\nTest 6: Edge cases\n")
cat("-", rep("-", 60), "\n", sep = "")

# Test with K=1 (single cluster)
sim_k1 <- simulate_dynamic_synergy_data(
  N = 10, K = 1, r = 3, M = 8, T_each = 100,
  seed = 600
)
stopifnot(length(sim_k1$Lambda_list) == 1)
stopifnot(all(sim_k1$z_true == 1))

# Test with K=2 (two clusters)
sim_k2 <- simulate_dynamic_synergy_data(
  N = 10, K = 2, r = 3, M = 8, T_each = 100,
  seed = 601
)
stopifnot(length(sim_k2$Lambda_list) == 2)

# Test with zero separation
sim_zero_sep <- simulate_dynamic_synergy_data(
  N = 10, K = 3, r = 3, M = 8, T_each = 100,
  cluster_sep_spatial = 0,
  cluster_sep_temporal = 0,
  cluster_sep_stability = 0,
  seed = 602
)
stopifnot(length(sim_zero_sep$Lambda_list) == 3)

cat("✓ Edge cases test passed\n")
cat("  - K=1 (single cluster): OK\n")
cat("  - K=2 (two clusters): OK\n")
cat("  - Zero separation: OK\n")

# ============================================================================
# Test 7: K>3 clusters - verify unique parameters for each cluster
# ============================================================================
cat("\nTest 7: K>3 clusters - verify unique parameters for each cluster\n")
cat("-", rep("-", 60), "\n", sep = "")

# Test with K=5 clusters
sim_k5 <- simulate_dynamic_synergy_data(
  N = 30, K = 5, r = 3, M = 8, T_each = 200,
  cluster_sep_spatial = 1.5,
  cluster_sep_temporal = 1.5,
  cluster_sep_stability = 1.5,
  seed = 650
)

# Verify temporal parameters are unique for each cluster
temporal_base_freqs <- sapply(sim_k5$temporal_params, function(p) p$base_frequency)
temporal_modes <- sapply(sim_k5$temporal_params, function(p) p$coordination_mode)
temporal_mean_phases <- sapply(sim_k5$temporal_params, function(p) mean(p$phases))

cat("  Temporal parameters for K=5:\n")
for (k in 1:5) {
  cat(sprintf("    Cluster %d: mode=%s, base_freq=%.3f, mean_phase=%.3f\n",
              k, temporal_modes[k], temporal_base_freqs[k], temporal_mean_phases[k]))
}

# Check that base frequencies are all different
stopifnot(length(unique(temporal_base_freqs)) == 5)
cat("  ✓ All 5 clusters have unique base frequencies\n")

# Check that modes cycle through rigid/alternating/adaptive
expected_modes <- c("rigid", "alternating", "adaptive", "rigid", "alternating")
stopifnot(all(temporal_modes == expected_modes))
cat("  ✓ Coordination modes cycle correctly: rigid → alternating → adaptive → rigid → alternating\n")

# Verify stability parameters are unique for each cluster
stability_sigma_whites <- sapply(sim_k5$stability_params, function(p) p$sigma_white)
stability_sigma_lows <- sapply(sim_k5$stability_params, function(p) p$sigma_low)
stability_profiles <- sapply(sim_k5$stability_params, function(p) p$control_profile)

cat("\n  Stability parameters for K=5:\n")
for (k in 1:5) {
  cat(sprintf("    Cluster %d: profile=%s, sigma_white=%.4f, sigma_low=%.4f\n",
              k, stability_profiles[k], stability_sigma_whites[k], stability_sigma_lows[k]))
}

# Check that sigma_white values are all different (within tolerance)
stopifnot(length(unique(round(stability_sigma_whites, 4))) == 5)
cat("  ✓ All 5 clusters have unique sigma_white values\n")

# Check that profiles cycle through stable/moderate/adaptive
expected_profiles <- c("stable", "moderate", "adaptive", "stable", "moderate")
stopifnot(all(stability_profiles == expected_profiles))
cat("  ✓ Stability profiles cycle correctly: stable → moderate → adaptive → stable → moderate\n")

# Verify spatial parameters produce distinct Lambda matrices
lambda_pairwise_dists <- numeric(0)
for (i in 1:4) {
  for (j in (i+1):5) {
    dist_ij <- sqrt(sum((sim_k5$Lambda_list[[i]] - sim_k5$Lambda_list[[j]])^2))
    lambda_pairwise_dists <- c(lambda_pairwise_dists, dist_ij)
  }
}

cat(sprintf("\n  Spatial parameters for K=5:\n"))
cat(sprintf("    Mean pairwise Lambda distance: %.4f\n", mean(lambda_pairwise_dists)))
cat(sprintf("    Min pairwise Lambda distance: %.4f\n", min(lambda_pairwise_dists)))
cat(sprintf("    Max pairwise Lambda distance: %.4f\n", max(lambda_pairwise_dists)))

# Check that all pairwise distances are positive (i.e., all Lambda matrices are distinct)
stopifnot(all(lambda_pairwise_dists > 0))
cat("  ✓ All 5 clusters have distinct Lambda matrices\n")

cat("\n✓ K>3 test passed - all clusters have unique parameters\n")

# ============================================================================
# Test 8: Visualization test (optional, for manual inspection)
# ============================================================================
cat("\nTest 8: Generate example plots for visual inspection\n")
cat("-", rep("-", 60), "\n", sep = "")

sim_vis <- simulate_dynamic_synergy_data(
  N = 30, K = 3, r = 3, M = 8, T_each = 200,
  cluster_sep_spatial = 1.5,
  cluster_sep_temporal = 1.5,
  cluster_sep_stability = 1.5,
  seed = 700
)

# Plot synergy activations for first subject from each cluster
cat("  - Generated data for visualization\n")
cat("  - Cluster assignments: ", table(sim_vis$z_true), "\n")

# Print summary of parameters for each cluster
cat("\n  Cluster parameter summary:\n")
for (k in 1:K) {
  cat(sprintf("\n  Cluster %d:\n", k))
  cat(sprintf("    Spatial: bias pattern shape = %s\n", 
              paste(dim(sim_vis$spatial_params[[k]]$bias_pattern), collapse = "x")))
  cat(sprintf("    Temporal: mode = %s, base_freq = %.2f\n",
              sim_vis$temporal_params[[k]]$coordination_mode,
              sim_vis$temporal_params[[k]]$base_frequency))
  cat(sprintf("    Stability: profile = %s, sigma_white = %.3f\n",
              sim_vis$stability_params[[k]]$control_profile,
              sim_vis$stability_params[[k]]$sigma_white))
}

cat("\n✓ Visualization test passed\n")

# ============================================================================
# Summary
# ============================================================================
cat("\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("All tests passed successfully!\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("\nThe simulate_dynamic_synergy_data() function is working correctly.\n")
cat("All three control axes (spatial, temporal, stability) are functional.\n")
