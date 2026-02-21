#' Run Parameter Sweep Experiment
#'
#' Orchestrates a comprehensive parameter sweep experiment to evaluate the
#' performance of MFA, MPCA, and baseline methods using the v2 dynamic synergy
#' simulator with three orthogonal control axes.
#'
#' @param test Logical. If TRUE, uses a smaller test grid for quick validation.
#'   Default is FALSE (full parameter grid).
#' @param cores Integer. Number of CPU cores to use for parallel processing.
#'   Default is 4. Set to 1 for sequential processing.
#' @param output_dir Character. Directory to save results. Default is
#'   "inst/paper/results_sweep" if it exists (repo use), otherwise
#'   "results_sweep" in the current working directory.
#' @param param_grid Optional. Custom parameter grid data frame. If NULL,
#'   uses the built-in full or test grid based on the \code{test} parameter.
#' @param verbose Logical. If TRUE, prints progress information. Default is TRUE.
#'
#' @return Invisibly returns a list containing:
#'   \itemize{
#'     \item \code{raw_results}: Data frame with raw results from all runs
#'     \item \code{aggregated_results}: Data frame with results aggregated across seeds
#'     \item \code{summary_file}: Path to the summary CSV file
#'     \item \code{win_rates_file}: Path to the win rates CSV file
#'     \item \code{output_dir}: Directory where results were saved
#'   }
#'
#' @details
#' The parameter sweep evaluates:
#' \itemize{
#'   \item 6 methods: SingleFA, SinglePCA, TwoStep_FA, TwoStep_PCA, MixtureFA, MixturePCA
#'   \item 2 model types: Mixture (static) vs Dynamic_v2 (with 3 orthogonal axes)
#'   \item Metrics: ARI (clustering accuracy), VAF (reconstruction quality), BIC, SSE
#'   \item Variable parameters: N, K, r, cluster_sep_spatial, cluster_sep_temporal, cluster_sep_stability
#' }
#'
#' For parallel processing (cores > 1), the function uses \code{parLapplyLB} for
#' load balancing. If the \code{pbapply} package is installed, a progress bar
#' will be displayed.
#'
#' @examples
#' \dontrun{
#' # Run test mode with 2 cores
#' results <- run_param_sweep(test = TRUE, cores = 2)
#'
#' # Run full sweep with 8 cores
#' results <- run_param_sweep(test = FALSE, cores = 8)
#'
#' # Run with custom output directory
#' results <- run_param_sweep(test = TRUE, cores = 2, output_dir = "my_results")
#'
#' # Run with custom parameter grid
#' library(dplyr)
#' custom_grid <- expand.grid(
#'   N = c(20, 50),
#'   K = c(2, 3),
#'   r = 2,
#'   M = 8,
#'   T_each = 200,
#'   cluster_sep_spatial = c(0.3, 0.6),
#'   cluster_sep_temporal = 0.6,
#'   cluster_sep_stability = 1.0,
#'   seed = 1:3,
#'   stringsAsFactors = FALSE
#' )
#' custom_grid$run_id <- seq_len(nrow(custom_grid))
#' 
#' results <- run_param_sweep(
#'   param_grid = custom_grid,
#'   cores = 8,
#'   output_dir = "inst/paper/results_custom"
#' )
#' 
#' # Access results directly
#' head(results$raw_results)
#' head(results$aggregated_results)
#' }
#'
#' @export
run_param_sweep <- function(test = FALSE, cores = 4, output_dir = NULL,
                            param_grid = NULL, verbose = TRUE) {
  
  # Load required libraries
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Please install it.")
  }
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required. Please install it.")
  }
  
  # Determine output directory
  if (is.null(output_dir)) {
    if (dir.exists("inst/paper")) {
      output_dir <- file.path("inst", "paper", "results_sweep")
      if (verbose) message("Using output directory: ", output_dir)
    } else {
      output_dir <- file.path(getwd(), "results_sweep")
      if (verbose) message("inst/paper not found. Using output directory: ", output_dir)
    }
  }
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Select or create parameter grid
  if (is.null(param_grid)) {
    if (test) {
      if (verbose) {
        cat("\n=== RUNNING TEST MODE ===\n")
        cat("Using smaller parameter grid for quick validation\n\n")
      }
      param_grid_to_use <- make_param_grid_test()
    } else {
      if (verbose) {
        cat("\n=== RUNNING FULL PARAMETER SWEEP ===\n")
        cat("This may take several hours to complete\n\n")
      }
      param_grid_to_use <- make_param_grid_full()
    }
  } else {
    param_grid_to_use <- param_grid
    if (verbose) {
      cat("\n=== RUNNING CUSTOM PARAMETER GRID ===\n\n")
    }
  }
  
  # Print configuration
  if (verbose) {
    cat("Configuration:\n")
    cat("  Parameter combinations:", nrow(param_grid_to_use), "\n")
    cat("  CPU cores:", cores, "\n")
    cat("  Output directory:", output_dir, "\n\n")
  }
  
  # Define the function to run a single parameter combination
  run_single_config <- function(i, param_grid) {
    p <- param_grid[i, ]
    
    tryCatch({
      # Run comparison using Dynamic_v2 model
      result <- compare_static_vs_dynamic(
        model_types = "Dynamic_v2",
        N = p$N,
        K = p$K,
        r = p$r,
        M = p$M,
        T_each = p$T_each,
        cluster_sep_spatial = p$cluster_sep_spatial,
        cluster_sep_temporal = p$cluster_sep_temporal,
        cluster_sep_stability = p$cluster_sep_stability,
        seed = p$seed,
        output_dir = tempdir(),  # Don't save individual files
        n_init = 1,
        use_kmeans_init = FALSE,
        mc_cores = 1  # Each worker uses 1 core
      )
      
      # Add run_id to track which parameter combination this is
      result$run_id <- p$run_id
      
      return(result)
      
    }, error = function(e) {
      warning("Error in run ", i, ": ", e$message)
      return(NULL)
    })
  }
  
  # Run parameter sweep
  if (verbose) cat("Starting parameter sweep...\n")
  start_time <- Sys.time()
  
  if (cores > 1) {
    if (verbose) cat("Using parallel processing with", cores, "cores\n\n")
    
    # Set up parallel cluster
    cl <- parallel::makeCluster(cores)
    
    # Ensure cluster is stopped even on error
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    # Export necessary objects and functions to workers
    parallel::clusterExport(cl, c("param_grid_to_use", "run_single_config"),
                           envir = environment())
    
    # Load synergyMixR on each worker
    parallel::clusterEvalQ(cl, {
      library(synergyMixR)
    })
    
    # Run parameter sweep in parallel with load balancing
    # Use pbapply for progress bar if available, otherwise use parLapplyLB
    idx <- seq_len(nrow(param_grid_to_use))
    
    if (requireNamespace("pbapply", quietly = TRUE)) {
      if (verbose) cat("Using pbapply for progress tracking\n\n")
      results_list <- pbapply::pblapply(
        idx,
        function(i) run_single_config(i, param_grid_to_use),
        cl = cl
      )
    } else {
      if (verbose) {
        cat("Note: Install 'pbapply' package for progress tracking\n")
        cat("Running without progress bar...\n\n")
      }
      results_list <- parallel::parLapplyLB(
        cl, idx, run_single_config,
        param_grid = param_grid_to_use
      )
    }
    
  } else {
    if (verbose) cat("Using sequential processing (single core)\n\n")
    
    # Run sequentially with progress tracking
    results_list <- vector("list", nrow(param_grid_to_use))
    
    for (i in seq_len(nrow(param_grid_to_use))) {
      results_list[[i]] <- run_single_config(i, param_grid_to_use)
      if (verbose) print_progress(i, nrow(param_grid_to_use), start_time)
    }
  }
  
  if (verbose) {
    cat("\nParameter sweep completed!\n")
    elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    cat("Total time:", round(elapsed_time, 1), "minutes\n\n")
  }
  
  # Combine results
  if (verbose) cat("Combining results...\n")
  results_list <- results_list[!sapply(results_list, is.null)]  # Remove NULL entries
  all_results <- dplyr::bind_rows(results_list)
  
  if (verbose) {
    cat("Total rows:", nrow(all_results), "\n")
    cat("Successful runs:", length(unique(all_results$run_id)), "/", 
        nrow(param_grid_to_use), "\n\n")
  }
  
  # Filter out non-finite values (NaN, Inf, -Inf) from numeric columns
  # This prevents warnings and errors in downstream aggregation and plotting
  if (verbose) cat("Filtering non-finite values...\n")
  n_rows_before <- nrow(all_results)
  
  # Identify numeric columns that should be finite
  numeric_cols <- c("ARI", "BIC", "VAF", "SSE")
  numeric_cols <- numeric_cols[numeric_cols %in% names(all_results)]
  
  if (length(numeric_cols) > 0) {
    # Filter to rows where all numeric metrics are finite
    all_results <- all_results %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(numeric_cols), is.finite))
    
    n_rows_after <- nrow(all_results)
    n_removed <- n_rows_before - n_rows_after
    
    if (verbose && n_removed > 0) {
      cat("Warning: Removed", n_removed, "rows with non-finite values (NaN/Inf/-Inf)\n")
      cat("   This typically occurs when SingleFA/SinglePCA fail to converge\n\n")
    } else if (verbose) {
      cat("All values are finite\n\n")
    }
  }
  
  # Save raw results
  if (verbose) cat("Saving raw results...\n")
  raw_files <- save_results(all_results, output_dir, prefix = "sweep_results_raw", 
                           timestamp = TRUE)
  
  # Aggregate results across seeds
  if (verbose) cat("Aggregating results across seeds...\n")
  aggregated_results <- aggregate_results(all_results)
  
  if (verbose) cat("Aggregated rows:", nrow(aggregated_results), "\n\n")
  
  # Save aggregated results
  if (verbose) cat("Saving aggregated results...\n")
  agg_files <- save_results(aggregated_results, output_dir, 
                           prefix = "sweep_results_aggregated", 
                           timestamp = TRUE)
  
  # Create summary statistics
  if (verbose) cat("Creating summary statistics...\n")
  
  # Summary by method and model type
  summary_ari <- create_summary_table(aggregated_results, metric = "ARI")
  summary_vaf <- create_summary_table(aggregated_results, metric = "VAF")
  summary_bic <- create_summary_table(aggregated_results, metric = "BIC")
  
  if (verbose) {
    cat("\nSummary - ARI (Clustering Accuracy):\n")
    print(summary_ari, n = Inf)
    
    cat("\nSummary - VAF (Variance Accounted For):\n")
    print(summary_vaf, n = Inf)
  }
  
  # Compute win rates
  if (verbose) cat("\nComputing win rates...\n")
  win_rates_ari <- compute_win_rates(aggregated_results, metric = "ARI", maximize = TRUE)
  win_rates_vaf <- compute_win_rates(aggregated_results, metric = "VAF", maximize = TRUE)
  win_rates_bic <- compute_win_rates(aggregated_results, metric = "BIC", maximize = FALSE)
  
  if (verbose) {
    cat("\nWin Rates - ARI (% of parameter combinations where method performed best):\n")
    print(win_rates_ari, n = Inf)
    
    cat("\nWin Rates - VAF:\n")
    print(win_rates_vaf, n = Inf)
    
    cat("\nWin Rates - BIC (lower is better):\n")
    print(win_rates_bic, n = Inf)
  }
  
  # Save summary tables
  summary_file <- file.path(output_dir, "sweep_summary.csv")
  write.csv(aggregated_results, summary_file, row.names = FALSE)
  if (verbose) cat("\nSummary saved to:", summary_file, "\n")
  
  # Save win rates
  win_rates_file <- file.path(output_dir, "win_rates.csv")
  win_rates_combined <- dplyr::bind_rows(
    win_rates_ari %>% dplyr::mutate(metric = "ARI"),
    win_rates_vaf %>% dplyr::mutate(metric = "VAF"),
    win_rates_bic %>% dplyr::mutate(metric = "BIC")
  )
  write.csv(win_rates_combined, win_rates_file, row.names = FALSE)
  if (verbose) cat("Win rates saved to:", win_rates_file, "\n")
  
  if (verbose) {
    cat("\n=== PARAMETER SWEEP COMPLETE ===\n")
    cat("Results saved to:", output_dir, "\n")
    cat("Next steps:\n")
    cat("  1. Run plot_param_sweep_results.R to generate figures\n")
    cat("  2. Review results in", output_dir, "\n\n")
  }
  
  # Return results invisibly for interactive use
  invisible(list(
    raw_results = all_results,
    aggregated_results = aggregated_results,
    summary_file = summary_file,
    win_rates_file = win_rates_file,
    output_dir = output_dir
  ))
}
