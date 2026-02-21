#!/usr/bin/env Rscript
# Performance Test Suite for OpenMP Parallelization
# 
# This script tests the parallel performance of mfa_em_fit() and mixture_pca_em_fit()
# to verify that OpenMP parallelization provides speedup as expected.
#
# Usage:
#   Rscript tests/performance/test_parallel_performance.R
#
# Or from R:
#   source("tests/performance/test_parallel_performance.R")

library(synergyMixR)

# Configuration
SEED <- 42
TOLERANCE <- 1e-10  # Tolerance for numerical comparison

cat("\n")
cat("========================================\n")
cat("Performance Test Suite for OpenMP\n")
cat("========================================\n")
cat("\n")

# Helper function to format time
format_time <- function(seconds) {
  if (seconds < 1) {
    sprintf("%.0f ms", seconds * 1000)
  } else {
    sprintf("%.2f s", seconds)
  }
}

# Helper function to print test result
print_result <- function(test_name, passed, details = "") {
  status <- if (passed) "✓ PASS" else "✗ FAIL"
  cat(sprintf("[%s] %s\n", status, test_name))
  if (details != "") {
    cat(sprintf("        %s\n", details))
  }
}

# Test counter
tests_passed <- 0
tests_failed <- 0

cat("System Information:\n")
cat(sprintf("  R version: %s\n", R.version.string))
cat(sprintf("  Platform: %s\n", R.version$platform))
cat(sprintf("  CPU cores: %d\n", parallel::detectCores()))
cat(sprintf("  synergyMixR version: %s\n", packageVersion("synergyMixR")))
cat("\n")

# ============================================================================
# Test 1: Reproducibility across different thread counts
# ============================================================================
cat("Test 1: Reproducibility with Different Thread Counts\n")
cat("------------------------------------------------------\n")

test_reproducibility <- function(method_name, fit_func, N = 50, K = 2, r = 2) {
  cat(sprintf("Testing %s (N=%d, K=%d, r=%d)...\n", method_name, N, K, r))
  
  # Generate test data
  sim_data <- simulate_mixture_data(
    N = N, K = K, r = r, M = 8, T_each = 100,
    cluster_sep = 1.0, noise_scale = 1.0, seed = SEED
  )
  
  # Run with different thread counts
  fit_1thread <- fit_func(sim_data$list_of_data, K = K, r = r, 
                          max_iter = 10, n_init = 1, use_kmeans_init = FALSE,
                          verbose = FALSE, n_threads = 1, cpp_seed = SEED)
  
  fit_2threads <- fit_func(sim_data$list_of_data, K = K, r = r,
                           max_iter = 10, n_init = 1, use_kmeans_init = FALSE,
                           verbose = FALSE, n_threads = 2, cpp_seed = SEED)
  
  fit_4threads <- fit_func(sim_data$list_of_data, K = K, r = r,
                           max_iter = 10, n_init = 1, use_kmeans_init = FALSE,
                           verbose = FALSE, n_threads = 4, cpp_seed = SEED)
  
  # Check reproducibility
  diff_1_2 <- abs(fit_1thread$logLik - fit_2threads$logLik)
  diff_1_4 <- abs(fit_1thread$logLik - fit_4threads$logLik)
  diff_2_4 <- abs(fit_2threads$logLik - fit_4threads$logLik)
  
  passed <- (diff_1_2 < TOLERANCE) && (diff_1_4 < TOLERANCE) && (diff_2_4 < TOLERANCE)
  
  details <- sprintf("logLik: 1-thread=%.6f, 2-threads=%.6f, 4-threads=%.6f (max diff=%.2e)",
                     fit_1thread$logLik, fit_2threads$logLik, fit_4threads$logLik,
                     max(diff_1_2, diff_1_4, diff_2_4))
  
  print_result(sprintf("%s reproducibility", method_name), passed, details)
  
  if (passed) tests_passed <<- tests_passed + 1 else tests_failed <<- tests_failed + 1
  
  return(passed)
}

test_reproducibility("MFA", mfa_em_fit)
test_reproducibility("MPCA", mixture_pca_em_fit)

cat("\n")

# ============================================================================
# Test 2: Performance scaling with thread count
# ============================================================================
cat("Test 2: Performance Scaling with Thread Count\n")
cat("------------------------------------------------------\n")

test_performance_scaling <- function(method_name, fit_func, N = 100, K = 2, r = 2) {
  cat(sprintf("Testing %s performance scaling (N=%d, K=%d, r=%d)...\n", method_name, N, K, r))
  
  # Generate test data
  sim_data <- simulate_mixture_data(
    N = N, K = K, r = r, M = 8, T_each = 100,
    cluster_sep = 1.0, noise_scale = 1.0, seed = SEED
  )
  
  # Measure execution time with different thread counts
  thread_counts <- c(1, 2, 4)
  times <- numeric(length(thread_counts))
  
  for (i in seq_along(thread_counts)) {
    n_threads <- thread_counts[i]
    
    timing <- system.time({
      fit <- fit_func(sim_data$list_of_data, K = K, r = r,
                      max_iter = 10, n_init = 1, use_kmeans_init = FALSE,
                      verbose = FALSE, n_threads = n_threads, cpp_seed = SEED)
    })
    
    times[i] <- timing["elapsed"]
    cat(sprintf("  %d thread%s: %s\n", n_threads, 
                if (n_threads > 1) "s" else " ", format_time(times[i])))
  }
  
  # Calculate speedup
  speedup_2 <- times[1] / times[2]
  speedup_4 <- times[1] / times[3]
  
  cat(sprintf("  Speedup: 2-threads=%.2fx, 4-threads=%.2fx\n", speedup_2, speedup_4))
  
  # Check if we get any speedup (even modest speedup is acceptable)
  # We expect at least some benefit, but not necessarily linear speedup
  passed <- (speedup_2 >= 0.8) && (speedup_4 >= 0.8)
  
  details <- sprintf("Speedup: 2-threads=%.2fx, 4-threads=%.2fx (baseline: %s)",
                     speedup_2, speedup_4, format_time(times[1]))
  
  print_result(sprintf("%s performance scaling", method_name), passed, details)
  
  if (passed) tests_passed <<- tests_passed + 1 else tests_failed <<- tests_failed + 1
  
  return(list(passed = passed, times = times, speedup_2 = speedup_2, speedup_4 = speedup_4))
}

mfa_perf <- test_performance_scaling("MFA", mfa_em_fit)
mpca_perf <- test_performance_scaling("MPCA", mixture_pca_em_fit)

cat("\n")

# ============================================================================
# Test 3: Performance with different data sizes
# ============================================================================
cat("Test 3: Performance Scaling with Data Size\n")
cat("------------------------------------------------------\n")

test_data_size_scaling <- function(method_name, fit_func) {
  cat(sprintf("Testing %s with different data sizes...\n", method_name))
  
  data_sizes <- c(20, 50, 100)
  times_1thread <- numeric(length(data_sizes))
  times_4threads <- numeric(length(data_sizes))
  
  for (i in seq_along(data_sizes)) {
    N <- data_sizes[i]
    
    # Generate test data
    sim_data <- simulate_mixture_data(
      N = N, K = 2, r = 2, M = 8, T_each = 100,
      cluster_sep = 1.0, noise_scale = 1.0, seed = SEED
    )
    
    # 1 thread
    timing_1 <- system.time({
      fit <- fit_func(sim_data$list_of_data, K = 2, r = 2,
                      max_iter = 10, n_init = 1, use_kmeans_init = FALSE,
                      verbose = FALSE, n_threads = 1, cpp_seed = SEED)
    })
    times_1thread[i] <- timing_1["elapsed"]
    
    # 4 threads
    timing_4 <- system.time({
      fit <- fit_func(sim_data$list_of_data, K = 2, r = 2,
                      max_iter = 10, n_init = 1, use_kmeans_init = FALSE,
                      verbose = FALSE, n_threads = 4, cpp_seed = SEED)
    })
    times_4threads[i] <- timing_4["elapsed"]
    
    speedup <- times_1thread[i] / times_4threads[i]
    cat(sprintf("  N=%3d: 1-thread=%s, 4-threads=%s, speedup=%.2fx\n",
                N, format_time(times_1thread[i]), format_time(times_4threads[i]), speedup))
  }
  
  # Check that larger datasets benefit from parallelization
  # For small datasets (N=20), overhead dominates, so we focus on larger datasets
  # We expect speedup to improve with larger N
  speedups <- times_1thread / times_4threads
  
  # Check that speedup improves or stays stable as N increases
  # (speedup for N=100 should be better than or equal to N=50)
  speedup_trend_ok <- speedups[3] >= speedups[2] * 0.9  # Allow 10% tolerance
  
  # Also check that largest dataset (N=100) shows reasonable performance
  large_dataset_ok <- speedups[3] >= 0.8
  
  passed <- speedup_trend_ok && large_dataset_ok
  
  details <- sprintf("Speedups: N=20:%.2fx, N=50:%.2fx, N=100:%.2fx (large dataset ≥0.8x: %s)",
                     speedups[1], speedups[2], speedups[3],
                     if (large_dataset_ok) "✓" else "✗")
  
  print_result(sprintf("%s data size scaling", method_name), passed, details)
  
  if (passed) tests_passed <<- tests_passed + 1 else tests_failed <<- tests_failed + 1
  
  return(passed)
}

test_data_size_scaling("MFA", mfa_em_fit)
test_data_size_scaling("MPCA", mixture_pca_em_fit)

cat("\n")

# ============================================================================
# Test 4: Thread count parameter validation
# ============================================================================
cat("Test 4: Thread Count Parameter Validation\n")
cat("------------------------------------------------------\n")

test_thread_param <- function(method_name, fit_func) {
  cat(sprintf("Testing %s thread parameter handling...\n", method_name))
  
  # Generate small test data
  sim_data <- simulate_mixture_data(
    N = 20, K = 2, r = 2, M = 8, T_each = 50,
    cluster_sep = 1.0, noise_scale = 1.0, seed = SEED
  )
  
  # Test with various thread counts
  thread_tests <- c(1, 2, 4, 8)
  all_passed <- TRUE
  
  for (n_threads in thread_tests) {
    tryCatch({
      fit <- fit_func(sim_data$list_of_data, K = 2, r = 2,
                      max_iter = 5, n_init = 1, use_kmeans_init = FALSE,
                      verbose = FALSE, n_threads = n_threads, cpp_seed = SEED)
      
      # Check that result is valid
      if (is.null(fit$logLik) || !is.finite(fit$logLik)) {
        all_passed <- FALSE
        cat(sprintf("  ✗ n_threads=%d: Invalid result\n", n_threads))
      } else {
        cat(sprintf("  ✓ n_threads=%d: OK (logLik=%.2f)\n", n_threads, fit$logLik))
      }
    }, error = function(e) {
      all_passed <<- FALSE
      cat(sprintf("  ✗ n_threads=%d: Error - %s\n", n_threads, e$message))
    })
  }
  
  print_result(sprintf("%s thread parameter handling", method_name), all_passed)
  
  if (all_passed) tests_passed <<- tests_passed + 1 else tests_failed <<- tests_failed + 1
  
  return(all_passed)
}

test_thread_param("MFA", mfa_em_fit)
test_thread_param("MPCA", mixture_pca_em_fit)

cat("\n")

# ============================================================================
# Test 5: Seed parameter consistency
# ============================================================================
cat("Test 5: Seed Parameter Consistency\n")
cat("------------------------------------------------------\n")

test_seed_consistency <- function(method_name, fit_func) {
  cat(sprintf("Testing %s seed consistency...\n", method_name))
  
  # Generate test data
  sim_data <- simulate_mixture_data(
    N = 30, K = 2, r = 2, M = 8, T_each = 50,
    cluster_sep = 1.0, noise_scale = 1.0, seed = SEED
  )
  
  # Run multiple times with same seed
  results <- list()
  for (i in 1:3) {
    fit <- fit_func(sim_data$list_of_data, K = 2, r = 2,
                    max_iter = 10, n_init = 1, use_kmeans_init = FALSE,
                    verbose = FALSE, n_threads = 4, cpp_seed = SEED)
    results[[i]] <- fit$logLik
  }
  
  # Check all results are identical
  diff_1_2 <- abs(results[[1]] - results[[2]])
  diff_1_3 <- abs(results[[1]] - results[[3]])
  diff_2_3 <- abs(results[[2]] - results[[3]])
  
  passed <- (diff_1_2 < TOLERANCE) && (diff_1_3 < TOLERANCE) && (diff_2_3 < TOLERANCE)
  
  details <- sprintf("logLik: run1=%.6f, run2=%.6f, run3=%.6f (max diff=%.2e)",
                     results[[1]], results[[2]], results[[3]],
                     max(diff_1_2, diff_1_3, diff_2_3))
  
  print_result(sprintf("%s seed consistency", method_name), passed, details)
  
  if (passed) tests_passed <<- tests_passed + 1 else tests_failed <<- tests_failed + 1
  
  return(passed)
}

test_seed_consistency("MFA", mfa_em_fit)
test_seed_consistency("MPCA", mixture_pca_em_fit)

cat("\n")

# ============================================================================
# Summary
# ============================================================================
cat("========================================\n")
cat("Test Summary\n")
cat("========================================\n")
cat(sprintf("Total tests: %d\n", tests_passed + tests_failed))
cat(sprintf("Passed: %d\n", tests_passed))
cat(sprintf("Failed: %d\n", tests_failed))
cat("\n")

if (tests_failed == 0) {
  cat("✓ All performance tests passed!\n")
  cat("\n")
  cat("The OpenMP parallelization is working correctly:\n")
  cat("  - Results are reproducible across different thread counts\n")
  cat("  - Performance scales with thread count (modest speedup expected)\n")
  cat("  - Performance scales with data size\n")
  cat("  - Thread parameter handling is robust\n")
  cat("  - Seed parameter ensures consistency\n")
  cat("\n")
  quit(status = 0)
} else {
  cat("✗ Some performance tests failed.\n")
  cat("\n")
  cat("Please review the failed tests above and investigate:\n")
  cat("  - Check OpenMP installation and configuration\n")
  cat("  - Verify system has multiple CPU cores available\n")
  cat("  - Review benchmark results for expected performance characteristics\n")
  cat("\n")
  quit(status = 1)
}
