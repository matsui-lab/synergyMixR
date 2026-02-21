#' Run Parameter Sweep for Baseline Method Comparisons (v2)
#'
#' Orchestrates a parameter sweep experiment to evaluate baseline method performance
#' (SingleFA, SinglePCA, TwoStep_FA, TwoStep_PCA, MixtureFA, MixturePCA) using the
#' v2 dynamic synergy simulator with three orthogonal control axes. This function
#' focuses on method-level benchmarking rather than model comparison.
#'
#' @param test Logical. If TRUE, uses a smaller test grid for quick validation.
#'   Default is FALSE (full parameter grid).
#' @param grid_type Character. Parameter grid type: "test", "paper_core", or
#'   "paper_stress". If NULL (default), uses the \code{test} parameter to
#'   determine the grid. When specified, overrides the \code{test} parameter.
#' @param cores Integer. Number of CPU cores to use for parallel processing.
#'   Default is 4. Set to 1 for sequential processing.
#' @param output_dir Character. Directory to save results. Default is
#'   "inst/paper/results_baselines_v2" if it exists (repo use), otherwise
#'   "results_baselines_v2" in the current working directory.
#' @param param_grid Optional. Custom parameter grid data frame. If NULL,
#'   uses the built-in full or test grid based on the \code{test} parameter.
#'   Required columns: N, K, r, M, T_each, cluster_sep_spatial,
#'   cluster_sep_temporal, cluster_sep_stability, seed, run_id.
#' @param verbose Logical. If TRUE, prints progress information. Default is TRUE.
#' @param align_basis Logical. If \code{TRUE} (default), applies Hungarian algorithm
#'   alignment in Step 1.5 of TwoStep methods to ensure corresponding components across
#'   subjects are in the same order before clustering. Set to \code{FALSE} to disable
#'   alignment and compare performance with/without this step in the parameter sweep.
#' @param refine_assignments Logical. If \code{TRUE} (default), performs a single
#'   reassignment pass after fitting cluster-specific models in Step 4 of TwoStep methods.
#'   Each subject is reassigned to the cluster that minimizes reconstruction error. Set to
#'   \code{FALSE} to disable reassignment for ablation studies.
#' @param parallel_strategy Character. Parallelization strategy:
#'   \itemize{
#'     \item "grid": Local multicore parallelization using parLapplyLB
#'     \item "sequential": No parallelization (single core)
#'     \item "pbs": PBS/OpenPBS array job submission for HPC clusters
#'     \item "future": Use future framework (supports local and HPC backends)
#'   }
#'   Default is "grid".
#' @param pbs_config List. PBS configuration options (only used when
#'   parallel_strategy = "pbs"). See \code{\link{submit_pbs_sweep}} for details.
#' @param pbs_submit Logical. If TRUE and parallel_strategy = "pbs", submit the
#'   job to PBS. If FALSE, only generate scripts (dry run). Default is FALSE.
#' @param pbs_wait Logical. If TRUE, wait for PBS job completion and aggregate
#'   results. Default is FALSE.
#' @param future_backend Character. Backend for future parallelization (only used
#'   when parallel_strategy = "future"):
#'   \itemize{
#'     \item "multisession": Local parallel R sessions (default)
#'     \item "multicore": Local forked processes (Unix only)
#'     \item "batchtools_torque": PBS/TORQUE scheduler
#'     \item "batchtools_slurm": SLURM scheduler
#'   }
#'   See \code{\link{setup_future_backend}} for all options.
#' @param future_resources List. HPC resource specifications for future.batchtools
#'   backends (e.g., list(walltime = "04:00:00", memory = "4gb")).
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
#' This function complements \code{\link{run_param_sweep}} by focusing on
#' method-level benchmarking rather than model comparison. It evaluates:
#' \itemize{
#'   \item 6 methods: SingleFA, SinglePCA, TwoStep_FA, TwoStep_PCA, MixtureFA, MixturePCA
#'   \item 1 model type: Dynamic_v2 (with 3 orthogonal axes)
#'   \item Metrics: ARI (clustering accuracy), VAF (reconstruction quality), BIC, SSE
#'   \item Variable parameters: N, K, r, cluster_sep_spatial, cluster_sep_temporal, cluster_sep_stability
#' }
#'
#' For parallel processing (cores > 1), the function uses \code{parLapplyLB} for
#' load balancing. If the \code{pbapply} package is installed, a progress bar
#' will be displayed.
#'
#' @section PBS Mode:
#' When \code{parallel_strategy = "pbs"}, this function generates PBS array job
#' scripts and optionally submits them. Each array task runs one parameter
#' combination. Results are saved as individual RDS files and can be aggregated
#' after job completion using \code{\link{aggregate_partial_results}}.
#'
#' @section Future Mode:
#' When \code{parallel_strategy = "future"}, this function uses the future
#' framework for flexible parallelization. Local backends (multisession,
#' multicore) work out of the box. HPC backends (batchtools_torque,
#' batchtools_slurm) require the future.batchtools package and appropriate
#' cluster configuration.
#'
#' @examples
#' \dontrun{
#' # Run test mode with 2 cores
#' results <- run_param_sweep_baselines(test = TRUE, cores = 2)
#'
#' # Run full sweep with 8 cores
#' results <- run_param_sweep_baselines(test = FALSE, cores = 8)
#'
#' # Run with custom output directory
#' results <- run_param_sweep_baselines(
#'   test = TRUE,
#'   cores = 2,
#'   output_dir = "my_baselines_results"
#' )
#'
#' # Use future with local multisession backend
#' results <- run_param_sweep_baselines(
#'   test = TRUE,
#'   parallel_strategy = "future",
#'   future_backend = "multisession",
#'   cores = 4
#' )
#'
#' # Use PBS for HPC cluster
#' results <- run_param_sweep_baselines(
#'   test = FALSE,
#'   parallel_strategy = "pbs",
#'   pbs_config = list(queue = "batch", walltime = "04:00:00", mem = "4gb"),
#'   pbs_submit = TRUE,
#'   pbs_wait = TRUE
#' )
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
#' results <- run_param_sweep_baselines(
#'   param_grid = custom_grid,
#'   cores = 8,
#'   output_dir = "inst/paper/results_baselines_custom"
#' )
#'
#' # Access results directly
#' head(results$raw_results)
#' head(results$aggregated_results)
#' }
#'
#' @seealso
#' \code{\link{compare_baselines}} for single-run baseline comparison,
#' \code{\link{run_param_sweep}} for model comparison parameter sweep,
#' \code{\link{plot_sweep_summary}} for visualization
#'
#' @export
run_param_sweep_baselines <- function(test = FALSE,
                                         grid_type = NULL,
                                         cores = 4,
                                         output_dir = NULL,
                                         param_grid = NULL,
                                         verbose = TRUE,
                                         align_basis = TRUE,
                                         refine_assignments = TRUE,
                                         parallel_strategy = c("grid", "sequential", "pbs", "future"),
                                         pbs_config = list(),
                                         pbs_submit = FALSE,
                                         pbs_wait = FALSE,
                                         future_backend = "multisession",
                                         future_resources = list()) {
  
  # Load required libraries
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Please install it.")
  }
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required. Please install it.")
  }

  parallel_strategy <- match.arg(parallel_strategy)

  # Validate grid_type if provided
  if (!is.null(grid_type)) {
    valid_grids <- c("test", "paper_core", "paper_stress")
    if (!grid_type %in% valid_grids) {
      stop("Invalid grid_type: ", grid_type,
           ". Valid options: ", paste(valid_grids, collapse = ", "))
    }
  }

  # Determine output directory
  if (is.null(output_dir)) {
    if (dir.exists("inst/paper")) {
      output_dir <- file.path("inst", "paper", "results_baselines_v2")
      if (verbose) message("Using output directory: ", output_dir)
    } else {
      output_dir <- file.path(getwd(), "results_baselines_v2")
      if (verbose) message("inst/paper not found. Using output directory: ", output_dir)
    }
  }

  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Determine effective grid type
  effective_grid <- if (!is.null(grid_type)) {
    grid_type
  } else if (test) {
    "test"
  } else {
    "paper_core"
  }

  # Get parameter grid
  if (is.null(param_grid)) {
    # Source config_param_grid.R to get grid functions
    config_file <- system.file("paper", "v2", "config_param_grid.R",
                               package = "synergyMixR")
    if (config_file == "" || !file.exists(config_file)) {
      config_file <- "inst/paper/v2/config_param_grid.R"
    }
    if (file.exists(config_file)) {
      source(config_file, local = TRUE)
    }

    # Use grid functions based on effective_grid
    param_grid <- switch(effective_grid,
      "test" = if (exists("make_param_grid_test", mode = "function"))
                 make_param_grid_test() else make_param_grid_baselines_test(),
      "paper_core" = if (exists("make_param_grid_paper_core", mode = "function"))
                       make_param_grid_paper_core() else make_param_grid_full(),
      "paper_stress" = if (exists("make_param_grid_paper_stress", mode = "function"))
                         make_param_grid_paper_stress() else make_param_grid_full(),
      make_param_grid_baselines_test()  # fallback
    )

    if (verbose) {
      message("\n=== RUNNING ", toupper(effective_grid), " GRID (BASELINES V2) ===")
      if (effective_grid == "test") {
        message("Using smaller parameter grid for quick validation\n")
      } else {
        message("This may take several hours to complete\n")
      }
    }
  } else {
    # Validate custom parameter grid
    required_cols <- c("N", "K", "r", "M", "T_each",
                      "cluster_sep_spatial", "cluster_sep_temporal", "cluster_sep_stability",
                      "seed", "run_id")
    missing_cols <- setdiff(required_cols, colnames(param_grid))
    if (length(missing_cols) > 0) {
      stop("Custom parameter grid is missing required columns: ",
           paste(missing_cols, collapse = ", "))
    }
    if (verbose) {
      message("\n=== RUNNING CUSTOM PARAMETER GRID (BASELINES V2) ===\n")
    }
  }

  # Print configuration
  if (verbose) {
    message("Configuration:")
    message("  Parameter combinations: ", nrow(param_grid))
    message("  Parallel strategy: ", parallel_strategy)
    if (parallel_strategy %in% c("grid", "sequential")) {
      message("  CPU cores: ", cores)
    }
    message("  Align basis: ", align_basis)
    message("  Refine assignments: ", refine_assignments)
    message("  Output directory: ", output_dir)
    if (parallel_strategy == "future") {
      message("  Future backend: ", future_backend)
    }
    if (parallel_strategy == "pbs") {
      message("  PBS submit: ", pbs_submit)
      message("  PBS wait: ", pbs_wait)
      if (!is.null(pbs_config$queue)) message("  PBS queue: ", pbs_config$queue)
      if (!is.null(pbs_config$walltime)) message("  PBS walltime: ", pbs_config$walltime)
      if (!is.null(pbs_config$mem)) message("  PBS memory: ", pbs_config$mem)
    }
    message("")
  }

  # ============================================================================
  # PBS MODE
  # ============================================================================
  if (parallel_strategy == "pbs") {
    if (verbose) {
      cat("\n=== PBS MODE ===\n")
      cat("Preparing PBS array job submission...\n\n")
    }

    # Prepare sweep configuration for baselines
    sweep_config <- list(
      align_basis = align_basis,
      refine_assignments = refine_assignments
    )

    # Submit PBS sweep
    result <- submit_pbs_sweep(
      param_grid = param_grid,
      output_dir = output_dir,
      pbs_config = pbs_config,
      sweep_type = "baselines",
      sweep_config = sweep_config,
      submit = pbs_submit,
      verbose = verbose
    )

    # If waiting for completion, aggregate results
    if (pbs_submit && pbs_wait && !is.null(result$job_id)) {
      if (verbose) cat("\nWaiting for PBS job completion...\n")
      wait_for_pbs_job(
        result$job_id,
        output_dir = output_dir,
        n_tasks = result$n_tasks,
        verbose = verbose
      )

      if (verbose) cat("\nAggregating results...\n")
      all_results <- aggregate_partial_results(output_dir, verbose = verbose)

      if (!is.null(all_results) && nrow(all_results) > 0) {
        result$raw_results <- all_results
        result$aggregated_results <- aggregate_results(all_results)
      } else {
        warning("No results found after job completion. Check PBS logs in: ",
                file.path(output_dir, "logs"))
      }
    }

    return(invisible(result))
  }

  # ============================================================================
  # FUTURE MODE
  # ============================================================================
  if (parallel_strategy == "future") {
    if (verbose) {
      cat("\n=== FUTURE MODE ===\n")
      cat("Using future framework for parallelization\n")
      cat("Backend:", future_backend, "\n\n")
    }

    # Check for required packages
    if (!check_future_available(batchtools = grepl("^batchtools_", future_backend))) {
      stop("Required future packages not available. Install with:\n",
           "  install.packages(c('future', 'future.apply'))\n",
           "For HPC backends also install:\n",
           "  install.packages('future.batchtools')")
    }

    # Define run function for future
    run_future_baselines_config <- function(p, align_basis, refine_assignments) {
      # Get K and r from param_grid (handle both naming conventions)
      K <- if (!is.null(p$K_true)) p$K_true else p$K
      r <- if (!is.null(p$r_true)) p$r_true else p$r

      tryCatch({
        result <- compare_baselines(
          N = p$N,
          K = K,
          r = r,
          M = p$M,
          T_each = p$T_each,
          cluster_sep_spatial = p$cluster_sep_spatial,
          cluster_sep_temporal = p$cluster_sep_temporal,
          cluster_sep_stability = p$cluster_sep_stability,
          seed = p$seed,
          n_init = 1,
          use_kmeans_init = FALSE,
          mc_cores = 1,
          align_basis = align_basis,
          refine_assignments = refine_assignments
        )
        result$run_id <- p$run_id
        return(result)
      }, error = function(e) {
        warning("Error in run ", p$run_id, ": ", e$message)
        return(NULL)
      })
    }

    # Run with future
    all_results <- run_sweep_future(
      param_grid = param_grid,
      run_func = run_future_baselines_config,
      align_basis = align_basis,
      refine_assignments = refine_assignments,
      future_backend = future_backend,
      workers = cores,
      resources = future_resources,
      verbose = verbose
    )

    # Filter out NULL results
    all_results <- all_results[!sapply(all_results, is.null)]

    if (length(all_results) == 0) {
      stop("No valid results were produced. Check error messages above.")
    }

    all_results <- dplyr::bind_rows(all_results)

    # Save and aggregate results
    return(finalize_baselines_results(
      all_results = all_results,
      output_dir = output_dir,
      verbose = verbose
    ))
  }

  # ============================================================================
  # LOCAL MODE (grid or sequential)
  # ============================================================================
  
  # Define the function to run a single parameter combination
  run_single_config <- function(i, param_grid) {
    p <- param_grid[i, ]
    
    tryCatch({
      # Run baseline comparison using compare_baselines
      result <- compare_baselines(
        N = p$N,
        K = p$K,
        r = p$r,
        M = p$M,
        T_each = p$T_each,
        cluster_sep_spatial = p$cluster_sep_spatial,
        cluster_sep_temporal = p$cluster_sep_temporal,
        cluster_sep_stability = p$cluster_sep_stability,
        seed = p$seed,
        n_init = 1,
        use_kmeans_init = FALSE,
        mc_cores = 1,  # Each worker uses 1 core
        align_basis = align_basis,
        refine_assignments = refine_assignments
      )
      
      # Add run_id to track which parameter combination this is
      result$run_id <- p$run_id
      
      return(result)
      
    }, error = function(e) {
      warning("Error in run ", i, ": ", e$message)
      return(NULL)
    })
  }
  
  # Run parameter sweep with parallel processing
  if (verbose) message("Starting baseline parameter sweep...")
  start_time <- Sys.time()
  
  if (cores > 1) {
    if (verbose) message("Using parallel processing with ", cores, " cores\n")
    
    # Set up parallel cluster
    cl <- parallel::makeCluster(cores)
    
    # Ensure cluster is stopped on exit
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    # Export necessary objects and functions to workers
    parallel::clusterExport(cl, c("param_grid", "run_single_config",
                                  "align_basis", "refine_assignments"),
                           envir = environment())
    
    # Load synergyMixR on each worker
    parallel::clusterEvalQ(cl, {
      library(synergyMixR)
    })
    
    # Run parameter sweep in parallel with load balancing
    if (requireNamespace("pbapply", quietly = TRUE)) {
      # Use pbapply for progress bar if available
      results_list <- pbapply::pblapply(
        seq_len(nrow(param_grid)),
        run_single_config,
        param_grid = param_grid,
        cl = cl
      )
    } else {
      # Use parLapplyLB for load balancing without progress bar
      results_list <- parallel::parLapplyLB(
        cl,
        seq_len(nrow(param_grid)),
        run_single_config,
        param_grid = param_grid
      )
    }
    
  } else {
    if (verbose) message("Using sequential processing\n")
    
    # Sequential processing
    results_list <- lapply(seq_len(nrow(param_grid)), run_single_config,
                          param_grid = param_grid)
  }
  
  # Combine results
  if (verbose) message("\nCombining results...")
  results_list <- Filter(Negate(is.null), results_list)

  if (length(results_list) == 0) {
    stop("No valid results were produced. Check error messages above.")
  }

  all_results <- dplyr::bind_rows(results_list)

  # Finalize and save results
  finalize_baselines_results(
    all_results = all_results,
    output_dir = output_dir,
    start_time = start_time,
    verbose = verbose
  )
}


#' Finalize and Save Baseline Results
#'
#' Helper function to save raw results, aggregate, and compute win rates.
#'
#' @param all_results Data frame with all results.
#' @param output_dir Output directory.
#' @param start_time Optional start time for timing calculation.
#' @param verbose Logical for verbose output.
#' @return Invisibly returns list with results.
#' @keywords internal
finalize_baselines_results <- function(all_results, output_dir, start_time = NULL,
                                       verbose = TRUE) {
  # Filter out NaN/Inf values (only check SSE and VAF)
  # Note: BIC can be NA for TwoStep methods (k-means), ARI can be NA for Single methods
  n_before <- nrow(all_results)
  all_results <- all_results[
    is.finite(all_results$SSE) &
    is.finite(all_results$VAF),
  ]
  n_after <- nrow(all_results)

  if (verbose && n_before > n_after) {
    message("Filtered out ", n_before - n_after, " rows with NaN/Inf values in SSE/VAF")
  }

  # Save raw results
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  raw_file <- file.path(output_dir, paste0("baseline_results_raw_", timestamp, ".rds"))
  saveRDS(all_results, raw_file)
  if (verbose) message("Raw results saved to: ", raw_file)

  # Also save as CSV
  raw_csv <- file.path(output_dir, paste0("baseline_results_raw_", timestamp, ".csv"))
  write.csv(all_results, raw_csv, row.names = FALSE)

  # Aggregate results across seeds
  if (verbose) message("Aggregating results across seeds...")
  aggregated_results <- aggregate_results(all_results)

  agg_file <- file.path(output_dir, paste0("baseline_results_aggregated_", timestamp, ".rds"))
  saveRDS(aggregated_results, agg_file)
  if (verbose) message("Aggregated results saved to: ", agg_file)

  # Also save aggregated as CSV
  agg_csv <- file.path(output_dir, "baseline_results_aggregated.csv")
  write.csv(aggregated_results, agg_csv, row.names = FALSE)

  # Create summary table
  if (verbose) message("Creating summary table...")
  summary_table <- create_summary_table(aggregated_results)

  summary_file <- file.path(output_dir, "baseline_summary.csv")
  write.csv(summary_table, summary_file, row.names = FALSE)
  if (verbose) message("Summary table saved to: ", summary_file)

  # Calculate win rates for multiple metrics
  if (verbose) message("Calculating win rates...")
  win_rates_ari <- compute_win_rates(aggregated_results, metric = "ARI", maximize = TRUE)
  win_rates_vaf <- compute_win_rates(aggregated_results, metric = "VAF", maximize = TRUE)
  win_rates_bic <- compute_win_rates(aggregated_results, metric = "BIC", maximize = FALSE)

  # Combine win rates with metric column
  win_rates_combined <- dplyr::bind_rows(
    win_rates_ari %>% dplyr::mutate(metric = "ARI"),
    win_rates_vaf %>% dplyr::mutate(metric = "VAF"),
    win_rates_bic %>% dplyr::mutate(metric = "BIC")
  )

  win_rates_file <- file.path(output_dir, "baseline_win_rates.csv")
  write.csv(win_rates_combined, win_rates_file, row.names = FALSE)
  if (verbose) message("Win rates saved to: ", win_rates_file)

  # Print timing information
  if (!is.null(start_time)) {
    end_time <- Sys.time()
    elapsed <- difftime(end_time, start_time, units = "mins")
    if (verbose) {
      message("\n=== BASELINE PARAMETER SWEEP COMPLETE ===")
      message("Total time: ", round(elapsed, 2), " minutes")
      message("Total runs: ", nrow(all_results))
      message("Output directory: ", output_dir)
    }
  } else {
    if (verbose) {
      message("\n=== BASELINE PARAMETER SWEEP COMPLETE ===")
      message("Total runs: ", nrow(all_results))
      message("Output directory: ", output_dir)
    }
  }

  # Return results invisibly
  invisible(list(
    raw_results = all_results,
    aggregated_results = aggregated_results,
    summary_file = summary_file,
    win_rates_file = win_rates_file,
    output_dir = output_dir
  ))
}


#' Create Parameter Grid for Baseline Sweep (Test)
#'
#' @return A data frame with parameter combinations for test sweep.
#' @keywords internal
make_param_grid_baselines_test <- function() {
  param_grid <- expand.grid(
    N = 20,
    K = 2,
    r = 2,
    M = 8,
    T_each = 100,
    cluster_sep_spatial = c(0.5, 1.0),
    cluster_sep_temporal = c(0.5, 1.0),
    cluster_sep_stability = c(0.5, 1.0),
    seed = 1:2,
    stringsAsFactors = FALSE
  )
  param_grid$run_id <- seq_len(nrow(param_grid))
  return(param_grid)
}


#' Run Single Baseline Task
#'
#' Runs a single baseline comparison task. Used by PBS workers and future backends.
#'
#' @param task_id Task ID (1-indexed).
#' @param param_grid_file Path to RDS file containing parameter grid.
#' @param output_dir Output directory for results.
#' @param align_basis Logical for basis alignment.
#' @param refine_assignments Logical for assignment refinement.
#' @param verbose Logical for verbose output.
#' @return Invisibly returns the result data frame.
#' @keywords internal
#' @export
run_single_task_baselines <- function(task_id, param_grid_file, output_dir,
                                      align_basis = TRUE, refine_assignments = TRUE,
                                      verbose = TRUE) {
  # Load parameter grid
  param_grid <- readRDS(param_grid_file)

  if (task_id < 1 || task_id > nrow(param_grid)) {
    stop("Invalid task_id: ", task_id, ". Must be between 1 and ", nrow(param_grid))
  }

  p <- param_grid[task_id, ]

  if (verbose) {
    cat("Running task", task_id, "of", nrow(param_grid), "\n")
    cat("Parameters: N=", p$N, ", K=", p$K, ", r=", p$r,
        ", seed=", p$seed, "\n")
  }

  # Get K and r (handle both naming conventions)
  K <- if (!is.null(p$K_true)) p$K_true else p$K
  r <- if (!is.null(p$r_true)) p$r_true else p$r

  # Run baseline comparison
  result <- tryCatch({
    res <- compare_baselines(
      N = p$N,
      K = K,
      r = r,
      M = p$M,
      T_each = p$T_each,
      cluster_sep_spatial = p$cluster_sep_spatial,
      cluster_sep_temporal = p$cluster_sep_temporal,
      cluster_sep_stability = p$cluster_sep_stability,
      seed = p$seed,
      n_init = 1,
      use_kmeans_init = FALSE,
      mc_cores = 1,
      align_basis = align_basis,
      refine_assignments = refine_assignments
    )
    res$run_id <- p$run_id
    res
  }, error = function(e) {
    warning("Error in task ", task_id, ": ", e$message)
    NULL
  })

  # Save result
  if (!is.null(result)) {
    result_file <- file.path(output_dir, "partial_results",
                             sprintf("result_task_%05d.rds", task_id))
    dir.create(dirname(result_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(result, result_file)
    if (verbose) cat("Result saved to:", result_file, "\n")
  }

  invisible(result)
}
