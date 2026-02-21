# Internal Helper Functions for Parameter Sweep
# These functions support the run_param_sweep() function

# Aggregate Results Across Seeds
aggregate_results <- function(results_df) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Please install it.")
  }
  
  # Safety filter: Remove rows with non-finite metric values
  # This prevents NaN/Inf from propagating into aggregated results
  # Note: Only filter on SSE and VAF, as BIC can be NA for TwoStep methods (k-means)
  # and ARI can be NA for Single methods (by design)
  metric_cols <- c("VAF", "SSE")
  metric_cols <- metric_cols[metric_cols %in% names(results_df)]
  
  if (length(metric_cols) > 0) {
    n_before <- nrow(results_df)
    results_df <- results_df %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(metric_cols), is.finite))
    n_after <- nrow(results_df)
    
    if (n_before > n_after) {
      message("Note: aggregate_results() filtered ", n_before - n_after, 
              " rows with non-finite SSE/VAF values")
    }
  }
  
  # Group by all parameters except seed and run_id
  grouping_vars <- c("model_type", "Method", "N", "K", "r", "M", "T_each",
                     "cluster_sep_spatial", "cluster_sep_temporal", 
                     "cluster_sep_stability")
  
  # Filter to only include grouping variables that exist in the data
  grouping_vars <- grouping_vars[grouping_vars %in% names(results_df)]
  
  # Aggregate metrics
  aggregated <- results_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) %>%
    dplyr::summarize(
      ARI_mean = mean(ARI, na.rm = TRUE),
      ARI_sd = sd(ARI, na.rm = TRUE),
      BIC_mean = mean(BIC, na.rm = TRUE),
      BIC_sd = sd(BIC, na.rm = TRUE),
      VAF_mean = mean(VAF, na.rm = TRUE),
      VAF_sd = sd(VAF, na.rm = TRUE),
      SSE_mean = mean(SSE, na.rm = TRUE),
      SSE_sd = sd(SSE, na.rm = TRUE),
      n_reps = dplyr::n(),
      .groups = "drop"
    )
  
  return(aggregated)
}

# Save Results to Multiple Formats
save_results <- function(results_df, output_dir, prefix = "sweep_results", 
                        timestamp = TRUE) {
  # Create output directory if it doesn't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Generate filename
  if (timestamp) {
    timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
    base_name <- paste0(prefix, "_", timestamp_str)
  } else {
    base_name <- prefix
  }
  
  # Save as RDS
  rds_file <- file.path(output_dir, paste0(base_name, ".rds"))
  saveRDS(results_df, rds_file)
  message("Results saved to: ", rds_file)
  
  # Save as CSV
  csv_file <- file.path(output_dir, paste0(base_name, ".csv"))
  write.csv(results_df, csv_file, row.names = FALSE)
  message("Results saved to: ", csv_file)
  
  return(invisible(list(rds = rds_file, csv = csv_file)))
}

# Create Summary Table
create_summary_table <- function(results_df, metric = "ARI") {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Please install it.")
  }
  
  metric_col <- paste0(metric, "_mean")
  metric_sd_col <- paste0(metric, "_sd")
  
  if (!metric_col %in% names(results_df)) {
    stop("Metric '", metric, "' not found in results. ",
         "Make sure to aggregate results first using aggregate_results().")
  }
  
  # Safety filter: Remove rows with non-finite values in the metric column
  n_before <- nrow(results_df)
  results_df <- results_df %>%
    dplyr::filter(is.finite(!!rlang::sym(metric_col)))
  n_after <- nrow(results_df)
  
  if (n_before > n_after) {
    message("Note: create_summary_table() filtered ", n_before - n_after, 
            " rows with non-finite ", metric, " values")
  }
  
  summary_table <- results_df %>%
    dplyr::group_by(model_type, Method) %>%
    dplyr::summarize(
      mean = mean(!!rlang::sym(metric_col), na.rm = TRUE),
      sd = sd(!!rlang::sym(metric_col), na.rm = TRUE),
      min = min(!!rlang::sym(metric_col), na.rm = TRUE),
      max = max(!!rlang::sym(metric_col), na.rm = TRUE),
      .groups = "drop"
    )
  
  return(summary_table)
}

# Extract Best Performing Method
extract_best_method <- function(results_df, metric = "ARI", maximize = TRUE) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Please install it.")
  }
  
  metric_col <- paste0(metric, "_mean")
  
  if (!metric_col %in% names(results_df)) {
    stop("Metric '", metric, "' not found in results. ",
         "Make sure to aggregate results first using aggregate_results().")
  }
  
  # Group by parameter combination
  grouping_vars <- c("model_type", "N", "K", "r", "M", "T_each",
                     "cluster_sep_spatial", "cluster_sep_temporal", 
                     "cluster_sep_stability")
  grouping_vars <- grouping_vars[grouping_vars %in% names(results_df)]
  
  if (maximize) {
    best_methods <- results_df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) %>%
      dplyr::filter(!!rlang::sym(metric_col) == max(!!rlang::sym(metric_col), na.rm = TRUE)) %>%
      dplyr::ungroup()
  } else {
    best_methods <- results_df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) %>%
      dplyr::filter(!!rlang::sym(metric_col) == min(!!rlang::sym(metric_col), na.rm = TRUE)) %>%
      dplyr::ungroup()
  }
  
  return(best_methods)
}

# Compute Win Rates
compute_win_rates <- function(results_df, metric = "ARI", maximize = TRUE) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Please install it.")
  }
  
  # Get best methods for each parameter combination
  best_methods <- extract_best_method(results_df, metric = metric, maximize = maximize)
  
  # Count wins for each method
  win_rates <- best_methods %>%
    dplyr::group_by(model_type, Method) %>%
    dplyr::summarize(
      n_wins = dplyr::n(),
      .groups = "drop"
    )
  
  # Calculate total number of parameter combinations
  grouping_vars <- c("model_type", "N", "K", "r", "M", "T_each",
                     "cluster_sep_spatial", "cluster_sep_temporal", 
                     "cluster_sep_stability")
  grouping_vars <- grouping_vars[grouping_vars %in% names(results_df)]
  
  total_combinations <- results_df %>%
    dplyr::select(dplyr::all_of(grouping_vars)) %>%
    dplyr::distinct() %>%
    nrow()
  
  # Add win rate percentage
  win_rates$win_rate <- win_rates$n_wins / total_combinations * 100
  
  # Sort by win rate
  win_rates <- win_rates %>%
    dplyr::arrange(dplyr::desc(win_rate))
  
  return(win_rates)
}

# Progress Tracker for Parameter Sweep
print_progress <- function(current, total, start_time, extra_info = NULL) {
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  pct_complete <- current / total * 100
  
  # Estimate time remaining
  if (current > 0) {
    avg_time_per_iter <- elapsed / current
    time_remaining <- avg_time_per_iter * (total - current)
    eta <- Sys.time() + time_remaining
    
    cat(sprintf("\r[%d/%d] %.1f%% complete | Elapsed: %.1f min | ETA: %s",
                current, total, pct_complete, elapsed / 60,
                format(eta, "%H:%M:%S")))
  } else {
    cat(sprintf("\r[%d/%d] %.1f%% complete | Elapsed: %.1f min",
                current, total, pct_complete, elapsed / 60))
  }
  
  if (!is.null(extra_info)) {
    cat(" |", extra_info)
  }
  
  if (current == total) {
    cat("\n")
  }
  
  flush.console()
}

# Create Parameter Grid (Full)
make_param_grid_full <- function() {
  # Define parameter ranges
  N_values <- c(20, 30, 50)
  K_values <- c(2, 3)
  r_values <- c(2, 3)
  M_fixed <- 8
  T_each_fixed <- 200
  cluster_sep_spatial_values <- c(0.3, 0.6, 1.0)
  cluster_sep_temporal_values <- c(0.3, 0.6, 1.0)
  cluster_sep_stability_values <- c(0.3, 0.6, 1.0)
  seed_values <- 1:5
  
  # Create the full parameter grid
  param_grid <- expand.grid(
    N = N_values,
    K = K_values,
    r = r_values,
    M = M_fixed,
    T_each = T_each_fixed,
    cluster_sep_spatial = cluster_sep_spatial_values,
    cluster_sep_temporal = cluster_sep_temporal_values,
    cluster_sep_stability = cluster_sep_stability_values,
    seed = seed_values,
    stringsAsFactors = FALSE
  )
  
  # Add a unique run ID for each parameter combination
  param_grid$run_id <- seq_len(nrow(param_grid))
  
  return(param_grid)
}

# Create Parameter Grid (Test)
make_param_grid_test <- function() {
  # Smaller test grid for quick validation
  param_grid_test <- expand.grid(
    N = c(20),
    K = c(2),
    r = c(2),
    M = 8,
    T_each = 100,
    cluster_sep_spatial = c(0.5, 1.0),
    cluster_sep_temporal = c(0.5, 1.0),
    cluster_sep_stability = c(0.5, 1.0),
    seed = 1:2,
    stringsAsFactors = FALSE
  )
  param_grid_test$run_id <- seq_len(nrow(param_grid_test))
  
  return(param_grid_test)
}
