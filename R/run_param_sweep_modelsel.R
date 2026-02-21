#' Run Parameter Sweep with Model Selection
#'
#' Orchestrates a comprehensive parameter sweep experiment that includes
#' model selection (K, r) for each parameter combination. This function
#' generates v2 dynamic synergy data, performs model selection using
#' specified strategies, and evaluates performance metrics.
#'
#' @param test Logical. If TRUE, uses a smaller test grid for quick validation.
#'   Default is FALSE (full parameter grid).
#' @param grid_type Character. Parameter grid type: "test", "paper_core", or
#'   "paper_stress". If NULL (default), uses the \code{test} parameter to
#'   determine the grid. When specified, overrides the \code{test} parameter.
#' @param mode Character. Selection mode controlling the search behavior:
#'   \itemize{
#'     \item "modelsel": Normal model selection (K,r grid search)
#'     \item "oracle_r": Use true r, only select K
#'     \item "oracle_Kr": Use true K and r (no selection)
#'   }
#'   Default is "modelsel".
#' @param cores Integer. Number of CPU cores for parallel processing.
#'   Default is 4. Set to 1 for sequential processing.
#' @param output_dir Character. Directory to save results. Default is
#'   "inst/paper/results_sweep_modelsel" if it exists, otherwise
#'   "results_sweep_modelsel" in the current working directory.
#' @param param_grid Optional. Custom parameter grid data frame. If NULL,
#'   uses the built-in full or test grid based on the \code{test} parameter.
#' @param K_candidates Integer vector. Candidate K values for model selection.
#'   Default is 1:5.
#' @param r_candidates Integer vector. Candidate r values for model selection.
#'   Default is 1:5.
#' @param strategies Character vector. Model selection strategies to evaluate.
#'   Options: "grid", "twostage". Default is c("grid", "twostage").
#' @param methods Character vector. Model types to evaluate.
#'   Options: "MFA", "MPCA". Default is c("MFA", "MPCA").
#' @param criterion Character. Model selection criterion: "ICL" or "BIC".
#'   Default is "ICL".
#' @param icl_temperature Numeric. Temperature parameter for ICL entropy scaling.
#'   Higher values favor more clusters. Default is 2.0.
#' @param parallel_strategy Character. Parallelization strategy:
#'   \itemize{
#'     \item "grid": Local multicore parallelization using mclapply
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
#' @param max_iter Maximum EM iterations. Default is 50.
#' @param nIterFA Sub-iterations for FA update. Default is 5.
#' @param nIterPCA Sub-iterations for PCA update. Default is 5.
#' @param n_init Number of random initializations. Default is 1.
#' @param tol Convergence tolerance. Default is 1e-3.
#' @param use_kmeans_init Whether to use k-means initialization. Default is TRUE.
#' @param n_threads Number of OpenMP threads. Default is 1.
#' @param save_models Logical. If TRUE, save fitted models. Default is FALSE.
#' @param save_selection_summaries Logical. If TRUE, save BIC summaries for
#'   each run. Default is TRUE.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return Invisibly returns a list containing:
#'   \itemize{
#'     \item \code{raw_results}: Data frame with raw results from all runs
#'     \item \code{aggregated_results}: Data frame with results aggregated across seeds
#'     \item \code{win_rates}: Data frame with win rates by method/strategy
#'     \item \code{output_dir}: Directory where results were saved
#'   }
#'
#' @details
#' The parameter sweep evaluates model selection performance across:
#' \itemize{
#'   \item Methods: MFA, MPCA
#'   \item Strategies: grid (full BIC grid search), twostage (K first, then r)
#'   \item Metrics: ARI (clustering accuracy), VAF (reconstruction quality),
#'     BIC, selection accuracy (K_correct, r_correct, Kr_correct)
#'   \item Variable parameters: N, K, r, cluster_sep_spatial, cluster_sep_temporal,
#'     cluster_sep_stability
#' }
#'
#' @section Output Schema:
#' The raw_results data frame contains the following columns:
#' \itemize{
#'   \item \code{run_id}, \code{seed}: Identifiers
#'   \item \code{N}, \code{K_true}, \code{r_true}, \code{M}, \code{T_each}: True parameters
#'   \item \code{cluster_sep_spatial}, \code{cluster_sep_temporal}, \code{cluster_sep_stability}: Separation parameters
#'   \item \code{method}: "MFA" or "MPCA"
#'   \item \code{selection_strategy}: "grid" or "twostage"
#'   \item \code{mode}: Selection mode ("modelsel", "oracle_r", or "oracle_Kr")
#'   \item \code{r_mode}: Always "global" (for future cluster-wise r extension)
#'   \item \code{K_hat}, \code{r_hat}: Estimated parameters
#'   \item \code{K_correct}, \code{r_correct}, \code{Kr_correct}: Selection accuracy
#'   \item \code{ARI}, \code{VAF}, \code{logLik}, \code{BIC}: Performance metrics
#'   \item \code{subspace_similarity}: Subspace recovery score (mean across clusters)
#'   \item \code{runtime_sec}: Runtime in seconds
#'   \item \code{error_msg}: Error message (empty string if successful)
#' }
#'
#' @examples
#' \dontrun{
#' # Run test mode with 2 cores
#' results <- run_param_sweep_modelsel(test = TRUE, cores = 2)
#'
#' # Run full sweep with 8 cores
#' results <- run_param_sweep_modelsel(test = FALSE, cores = 8)
#'
#' # Run with custom strategies
#' results <- run_param_sweep_modelsel(
#'   test = TRUE, cores = 4,
#'   strategies = c("twostage"),
#'   methods = c("MFA")
#' )
#' }
#'
#' @export
run_param_sweep_modelsel <- function(
    test = FALSE,
    grid_type = NULL,
    mode = c("modelsel", "oracle_r", "oracle_Kr"),
    cores = 4,
    output_dir = NULL,
    param_grid = NULL,
    K_candidates = 1:5,
    r_candidates = 1:5,
    strategies = c("grid", "twostage"),
    methods = c("MFA", "MPCA"),
    criterion = c("ICL", "BIC"),
    icl_temperature = 2.0,
    parallel_strategy = c("grid", "sequential", "pbs", "future"),
    pbs_config = list(),
    pbs_submit = FALSE,
    pbs_wait = FALSE,
    future_backend = "multisession",
    future_resources = list(),
    max_iter = 50,
    nIterFA = 5,
    nIterPCA = 5,
    n_init = 1,
    tol = 1e-3,
    use_kmeans_init = TRUE,
    n_threads = 1,
    save_models = FALSE,
    save_selection_summaries = TRUE,
    verbose = TRUE
) {
  # Load required packages
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Please install it.")
  }
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required. Please install it.")
  }
  
  parallel_strategy <- match.arg(parallel_strategy)
  mode <- match.arg(mode)
  criterion <- match.arg(criterion)

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
      output_dir <- file.path("inst", "paper", "results_sweep_modelsel")
    } else {
      output_dir <- file.path(getwd(), "results_sweep_modelsel")
    }
  }
  
  # Create output directories
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (save_selection_summaries) {
    dir.create(file.path(output_dir, "selection_summaries"), 
               recursive = TRUE, showWarnings = FALSE)
  }
  dir.create(file.path(output_dir, "figures"), 
             recursive = TRUE, showWarnings = FALSE)
  
  if (verbose) {
    message("Output directory: ", output_dir)
  }
  
  # Select or create parameter grid
  if (is.null(param_grid)) {
    # Determine effective grid type
    effective_grid <- if (!is.null(grid_type)) {
      grid_type
    } else if (test) {
      "test"
    } else {
      "paper_core"
    }
    
    if (verbose) {
      cat("\n=== RUNNING", toupper(effective_grid), "GRID ===\n")
      cat("Selection mode:", mode, "\n\n")
    }
    
    # Source config_param_grid.R to get grid functions
    # Try system.file first (for installed package), then fall back to relative path
    config_file <- system.file("paper", "v2", "config_param_grid.R",
                               package = "synergyMixR")
    if (config_file == "" || !file.exists(config_file)) {
      # Fall back to relative path (for running from repo checkout)
      config_file <- "inst/paper/v2/config_param_grid.R"
    }
    if (file.exists(config_file)) {
      source(config_file, local = TRUE)
    }
    
    # Use grid functions from config_param_grid.R
    param_grid <- switch(effective_grid,
      "test" = if (exists("make_param_grid_test", mode = "function"))
                 make_param_grid_test() else make_param_grid_modelsel_test(),
      "paper_core" = if (exists("make_param_grid_paper_core", mode = "function"))
                       make_param_grid_paper_core() else make_param_grid_modelsel_full(),
      "paper_stress" = if (exists("make_param_grid_paper_stress", mode = "function"))
                         make_param_grid_paper_stress() else make_param_grid_modelsel_full(),
      make_param_grid_modelsel_test()  # fallback
    )
  }

  # Handle PBS mode - submit jobs and return early
  if (parallel_strategy == "pbs") {
    if (verbose) {
      cat("\n=== PBS MODE ===\n")
      cat("Preparing PBS array job submission...\n\n")
    }

    # Prepare sweep configuration
    sweep_config <- list(
      methods = methods,
      strategies = strategies,
      K_candidates = K_candidates,
      r_candidates = r_candidates,
      mode = mode,
      criterion = criterion,
      icl_temperature = icl_temperature
    )

    # Submit PBS sweep
    result <- submit_pbs_sweep(
      param_grid = param_grid,
      output_dir = output_dir,
      pbs_config = pbs_config,
      sweep_type = "modelsel",
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
        result$aggregated_results <- aggregate_modelsel_results(all_results)
      } else {
        warning("No results found after job completion. Check PBS logs in: ",
                file.path(output_dir, "logs"))
      }
    }

    return(invisible(result))
  }

  # Handle future mode - use future framework for parallelization
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
    run_future_config <- function(p, methods, strategies,
                                  K_candidates, r_candidates, selection_mode,
                                  criterion, icl_temperature,
                                  max_iter, nIterFA, nIterPCA, n_init, tol,
                                  use_kmeans_init, n_threads,
                                  save_selection_summaries, output_dir) {

      # Get true K and r from param_grid (handle both naming conventions)
      K_true <- if (!is.null(p$K_true)) p$K_true else p$K
      r_true <- if (!is.null(p$r_true)) p$r_true else p$r

      run_single_task_modelsel(
        p = p,
        K_true = K_true,
        r_true = r_true,
        methods = methods,
        strategies = strategies,
        K_candidates = K_candidates,
        r_candidates = r_candidates,
        mode = selection_mode,
        criterion = criterion,
        icl_temperature = icl_temperature,
        max_iter = max_iter,
        tol = tol,
        n_init = n_init,
        use_kmeans_init = use_kmeans_init,
        verbose = FALSE
      )
    }

    # Run with future
    all_results <- run_sweep_future(
      param_grid = param_grid,
      run_func = run_future_config,
      methods = methods,
      strategies = strategies,
      K_candidates = K_candidates,
      r_candidates = r_candidates,
      selection_mode = mode,
      criterion = criterion,
      icl_temperature = icl_temperature,
      max_iter = max_iter,
      nIterFA = nIterFA,
      nIterPCA = nIterPCA,
      n_init = n_init,
      tol = tol,
      use_kmeans_init = use_kmeans_init,
      n_threads = n_threads,
      save_selection_summaries = save_selection_summaries,
      output_dir = output_dir,
      future_backend = future_backend,
      workers = cores,
      resources = future_resources,
      verbose = verbose
    )

    # Save and aggregate results (same as other modes)
    if (verbose) cat("Saving raw results...\n")
    timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
    raw_file <- file.path(output_dir, paste0("sweep_modelsel_raw_", timestamp_str, ".rds"))
    saveRDS(all_results, raw_file)

    raw_csv <- file.path(output_dir, paste0("sweep_modelsel_raw_", timestamp_str, ".csv"))
    write.csv(all_results, raw_csv, row.names = FALSE)

    if (verbose) cat("Raw results saved to:", raw_file, "\n\n")

    # Aggregate results
    if (verbose) cat("Aggregating results...\n")
    aggregated_results <- aggregate_modelsel_results(all_results)

    agg_csv <- file.path(output_dir, "sweep_modelsel_aggregated.csv")
    write.csv(aggregated_results, agg_csv, row.names = FALSE)

    # Compute win rates
    if (verbose) cat("Computing win rates...\n")
    win_rates <- compute_modelsel_win_rates(aggregated_results)

    win_rates_file <- file.path(output_dir, "win_rates.csv")
    write.csv(win_rates, win_rates_file, row.names = FALSE)

    if (verbose) {
      cat("\n=== PARAMETER SWEEP COMPLETE ===\n")
      cat("Results saved to:", output_dir, "\n")
    }

    return(invisible(list(
      raw_results = all_results,
      aggregated_results = aggregated_results,
      win_rates = win_rates,
      output_dir = output_dir
    )))
  }

  # Windows parallel fallback (for grid mode)
  if (.Platform$OS.type == "windows" && cores > 1 && parallel_strategy == "grid") {
    warning("Windows detected: forcing sequential execution (cores = 1)")
    cores <- 1
    parallel_strategy <- "sequential"
  }
  
  # Parallel safety: if using grid-level parallelization, force n_threads = 1
  if (parallel_strategy == "grid" && cores > 1 && n_threads > 1) {
    if (verbose) {
      message("Grid-level parallelization enabled: forcing n_threads = 1")
    }
    n_threads <- 1
  }
  
  # Print configuration
  if (verbose) {
    cat("Configuration:\n")
    cat("  Parameter combinations:", nrow(param_grid), "\n")
    cat("  Selection mode:", mode, "\n")
    cat("  K selection criterion:", criterion, "\n")
    if (criterion == "ICL") {
      cat("  ICL temperature:", icl_temperature, "\n")
    }
    cat("  Methods:", paste(methods, collapse = ", "), "\n")
    cat("  Strategies:", paste(strategies, collapse = ", "), "\n")
    cat("  K candidates:", paste(K_candidates, collapse = ", "), "\n")
    cat("  r candidates:", paste(r_candidates, collapse = ", "), "\n")
    cat("  CPU cores:", cores, "\n")
    cat("  Parallel strategy:", parallel_strategy, "\n\n")
  }
  
  # Define function to run a single parameter combination
  run_single_config <- function(i, param_grid, methods, strategies,
                                K_candidates, r_candidates, selection_mode,
                                criterion, icl_temperature,
                                max_iter, nIterFA, nIterPCA, n_init, tol,
                                use_kmeans_init, n_threads,
                                save_selection_summaries, output_dir) {
    p <- param_grid[i, ]
    results_list <- list()
    
    # Get true K and r from param_grid (handle both naming conventions)
    K_true <- if (!is.null(p$K_true)) p$K_true else p$K
    r_true <- if (!is.null(p$r_true)) p$r_true else p$r
    
    tryCatch({
      # Generate v2 synergy data
      data <- generate_synergy_data(
        N = p$N,
        K = K_true,
        r = r_true,
        M = p$M,
        T_each = p$T_each,
        cluster_sep_spatial = p$cluster_sep_spatial,
        cluster_sep_temporal = p$cluster_sep_temporal,
        cluster_sep_stability = p$cluster_sep_stability,
        seed = p$seed
      )
      
      list_of_data <- data$list_of_data
      true_cluster <- data$true_cluster
      true_Lambda_list <- data$Lambda_list
      
      # Determine effective K and r candidates based on selection mode
      effective_K_candidates <- switch(selection_mode,
        "oracle_Kr" = K_true,
        K_candidates
      )
      effective_r_candidates <- switch(selection_mode,
        "oracle_r" = r_true,
        "oracle_Kr" = r_true,
        r_candidates
      )
      
      # Run model selection for each method x strategy combination
      for (method in methods) {
        for (strategy in strategies) {
          result_row <- tryCatch({
            # Time the model selection
            start_time <- Sys.time()
            
            # Run model selection
            select_result <- select_K_r(
              list_of_data = list_of_data,
              method = method,
              strategy = strategy,
              Kvec = effective_K_candidates,
              rvec = effective_r_candidates,
              criterion = criterion,
              icl_temperature = icl_temperature,
              max_iter = max_iter,
              nIterFA = nIterFA,
              nIterPCA = nIterPCA,
              tol = tol,
              n_init = n_init,
              use_kmeans_init = use_kmeans_init,
              mc_cores_grid = 1,  # No nested parallelization
              mc_cores = 1,
              n_threads = n_threads,
              verbose = FALSE
            )
            
            runtime_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
            
            # Evaluate results
            eval_result <- evaluate_model_selection(
              list_of_data = list_of_data,
              true_cluster = true_cluster,
              K_true = K_true,
              r_true = r_true,
              select_result = select_result,
              runtime_sec = runtime_sec
            )
            
            # Compute subspace similarity
            subspace_sim <- NA_real_
            best_model <- select_result$best_model
            # MFA uses Lambda, MPCA uses W
            loading_matrix <- if (!is.null(best_model$Lambda)) {
              best_model$Lambda
            } else if (!is.null(best_model$W)) {
              best_model$W
            } else {
              NULL
            }
            if (!is.null(best_model) && !is.null(loading_matrix)) {
              # Get estimated cluster labels
              pred_cluster <- extract_cluster_labels(best_model)
              if (!any(is.na(pred_cluster))) {
                # Align clusters to truth
                alignment <- align_clusters_to_truth(true_cluster, pred_cluster)
                if (!is.na(alignment$accuracy)) {
                  # Compute cluster-wise subspace score
                  subspace_result <- compute_clusterwise_subspace_score(
                    true_Lambda_list,
                    loading_matrix,
                    alignment$mapping
                  )
                  subspace_sim <- subspace_result$mean_score
                }
              }
            }
            
            # Save selection summary if requested
            if (save_selection_summaries) {
              summary_file <- file.path(
                output_dir, "selection_summaries",
                sprintf("summary_run%d_seed%d_%s_%s.rds",
                        p$run_id, p$seed, method, strategy)
              )
              saveRDS(list(
                params = as.list(p),
                method = method,
                strategy = strategy,
                mode = selection_mode,
                summary_df = select_result$summary_df
              ), summary_file)
            }
            
            # Build result row
            data.frame(
              run_id = p$run_id,
              seed = p$seed,
              N = p$N,
              K_true = K_true,
              r_true = r_true,
              M = p$M,
              T_each = p$T_each,
              cluster_sep_spatial = p$cluster_sep_spatial,
              cluster_sep_temporal = p$cluster_sep_temporal,
              cluster_sep_stability = p$cluster_sep_stability,
              method = method,
              selection_strategy = strategy,
              mode = selection_mode,
              r_mode = "global",
              K_hat = eval_result$K_hat,
              r_hat = eval_result$r_hat,
              K_correct = eval_result$K_correct,
              r_correct = eval_result$r_correct,
              Kr_correct = eval_result$Kr_correct,
              ARI = eval_result$ARI,
              VAF = eval_result$VAF,
              logLik = eval_result$logLik,
              BIC = eval_result$BIC,
              subspace_similarity = subspace_sim,
              runtime_sec = eval_result$runtime_sec,
              error_msg = "",
              stringsAsFactors = FALSE
            )
          }, error = function(e) {
            # Return error row
            data.frame(
              run_id = p$run_id,
              seed = p$seed,
              N = p$N,
              K_true = K_true,
              r_true = r_true,
              M = p$M,
              T_each = p$T_each,
              cluster_sep_spatial = p$cluster_sep_spatial,
              cluster_sep_temporal = p$cluster_sep_temporal,
              cluster_sep_stability = p$cluster_sep_stability,
              method = method,
              selection_strategy = strategy,
              mode = selection_mode,
              r_mode = "global",
              K_hat = NA_integer_,
              r_hat = NA_integer_,
              K_correct = NA,
              r_correct = NA,
              Kr_correct = NA,
              ARI = NA_real_,
              VAF = NA_real_,
              logLik = NA_real_,
              BIC = NA_real_,
              subspace_similarity = NA_real_,
              runtime_sec = NA_real_,
              error_msg = as.character(e$message),
              stringsAsFactors = FALSE
            )
          })
          
          results_list[[length(results_list) + 1]] <- result_row
        }
      }
      
      dplyr::bind_rows(results_list)
      
    }, error = function(e) {
      # Return error rows for all method x strategy combinations
      error_rows <- list()
      for (method in methods) {
        for (strategy in strategies) {
          error_rows[[length(error_rows) + 1]] <- data.frame(
            run_id = p$run_id,
            seed = p$seed,
            N = p$N,
            K_true = K_true,
            r_true = r_true,
            M = p$M,
            T_each = p$T_each,
            cluster_sep_spatial = p$cluster_sep_spatial,
            cluster_sep_temporal = p$cluster_sep_temporal,
            cluster_sep_stability = p$cluster_sep_stability,
            method = method,
            selection_strategy = strategy,
            mode = selection_mode,
            r_mode = "global",
            K_hat = NA_integer_,
            r_hat = NA_integer_,
            K_correct = NA,
            r_correct = NA,
            Kr_correct = NA,
            ARI = NA_real_,
            VAF = NA_real_,
            logLik = NA_real_,
            BIC = NA_real_,
            subspace_similarity = NA_real_,
            runtime_sec = NA_real_,
            error_msg = as.character(e$message),
            stringsAsFactors = FALSE
          )
        }
      }
      dplyr::bind_rows(error_rows)
    })
  }
  
  # Run parameter sweep
  if (verbose) cat("Starting parameter sweep...\n")
  start_time <- Sys.time()
  
  if (parallel_strategy == "grid" && cores > 1) {
    if (verbose) cat("Using parallel processing with", cores, "cores\n\n")
    
    # Run in parallel using mclapply
    results_list <- parallel::mclapply(
      seq_len(nrow(param_grid)),
      function(i) {
        run_single_config(
          i, param_grid, methods, strategies,
          K_candidates, r_candidates, mode,
          criterion, icl_temperature,
          max_iter, nIterFA, nIterPCA, n_init, tol,
          use_kmeans_init, n_threads,
          save_selection_summaries, output_dir
        )
      },
      mc.cores = cores
    )
  } else {
    if (verbose) cat("Using sequential processing\n\n")
    
    # Run sequentially with progress
    results_list <- vector("list", nrow(param_grid))
    for (i in seq_len(nrow(param_grid))) {
      results_list[[i]] <- run_single_config(
        i, param_grid, methods, strategies,
        K_candidates, r_candidates, mode,
        criterion, icl_temperature,
        max_iter, nIterFA, nIterPCA, n_init, tol,
        use_kmeans_init, n_threads,
        save_selection_summaries, output_dir
      )
      
      if (verbose) {
        print_progress(i, nrow(param_grid), start_time)
      }
    }
  }
  
  if (verbose) {
    cat("\nParameter sweep completed!\n")
    elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    cat("Total time:", round(elapsed_time, 1), "minutes\n\n")
  }
  
  # Combine results
  if (verbose) cat("Combining results...\n")
  all_results <- dplyr::bind_rows(results_list)
  
  if (verbose) {
    cat("Total rows:", nrow(all_results), "\n")
    n_errors <- sum(all_results$error_msg != "", na.rm = TRUE)
    if (n_errors > 0) {
      cat("Runs with errors:", n_errors, "\n")
    }
    cat("\n")
  }
  
  # Save raw results
  if (verbose) cat("Saving raw results...\n")
  timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
  raw_file <- file.path(output_dir, paste0("sweep_modelsel_raw_", timestamp_str, ".rds"))
  saveRDS(all_results, raw_file)
  
  raw_csv <- file.path(output_dir, paste0("sweep_modelsel_raw_", timestamp_str, ".csv"))
  write.csv(all_results, raw_csv, row.names = FALSE)
  
  if (verbose) cat("Raw results saved to:", raw_file, "\n\n")
  
  # Aggregate results
  if (verbose) cat("Aggregating results...\n")
  aggregated_results <- aggregate_modelsel_results(all_results)
  
  agg_file <- file.path(output_dir, paste0("sweep_modelsel_aggregated_", timestamp_str, ".rds"))
  saveRDS(aggregated_results, agg_file)
  
  agg_csv <- file.path(output_dir, "sweep_modelsel_aggregated.csv")
  write.csv(aggregated_results, agg_csv, row.names = FALSE)
  
  if (verbose) cat("Aggregated results saved to:", agg_csv, "\n\n")
  
  # Compute win rates
  if (verbose) cat("Computing win rates...\n")
  win_rates <- compute_modelsel_win_rates(aggregated_results)
  
  win_rates_file <- file.path(output_dir, "win_rates.csv")
  write.csv(win_rates, win_rates_file, row.names = FALSE)
  
  if (verbose) {
    cat("Win rates saved to:", win_rates_file, "\n\n")
    cat("=== Win Rates Summary ===\n")
    print(win_rates)
    cat("\n")
  }
  
  if (verbose) {
    cat("=== PARAMETER SWEEP COMPLETE ===\n")
    cat("Results saved to:", output_dir, "\n")
    cat("Next steps:\n")
    cat("  1. Run plot_param_sweep_modelsel.R to generate figures\n")
    cat("  2. Review results in", output_dir, "\n\n")
  }
  
  # Return results invisibly
  invisible(list(
    raw_results = all_results,
    aggregated_results = aggregated_results,
    win_rates = win_rates,
    output_dir = output_dir
  ))
}


#' Create Parameter Grid for Model Selection Sweep (Full)
#'
#' @return A data frame with parameter combinations for full sweep.
#' @keywords internal
make_param_grid_modelsel_full <- function() {
  param_grid <- expand.grid(
    N = c(20, 30, 50),
    K = c(2, 3),
    r = c(2, 3),
    M = 8,
    T_each = 200,
    cluster_sep_spatial = c(0.5, 1.0, 1.5),
    cluster_sep_temporal = c(0.5, 1.0, 1.5),
    cluster_sep_stability = c(0.5, 1.0, 1.5),
    seed = 1:5,
    stringsAsFactors = FALSE
  )
  param_grid$run_id <- seq_len(nrow(param_grid))
  return(param_grid)
}


#' Create Parameter Grid for Model Selection Sweep (Test)
#'
#' @return A data frame with parameter combinations for test sweep.
#' @keywords internal
make_param_grid_modelsel_test <- function() {
  param_grid <- expand.grid(
    N = c(20),
    K = c(2),
    r = c(2),
    M = 8,
    T_each = 100,
    cluster_sep_spatial = c(0.5, 1.0),
    cluster_sep_temporal = c(1.0),
    cluster_sep_stability = c(1.0),
    seed = 1:2,
    stringsAsFactors = FALSE
  )
  param_grid$run_id <- seq_len(nrow(param_grid))
  return(param_grid)
}


#' Aggregate Model Selection Results
#'
#' @param results_df Raw results data frame from sweep.
#' @return Aggregated data frame.
#' @keywords internal
aggregate_modelsel_results <- function(results_df) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Please install it.")
  }
  
  # Group by all parameters except seed and run_id
  grouping_vars <- c("N", "K_true", "r_true", "M", "T_each",
                     "cluster_sep_spatial", "cluster_sep_temporal",
                     "cluster_sep_stability", "method", "selection_strategy",
                     "mode", "r_mode")
  
  grouping_vars <- grouping_vars[grouping_vars %in% names(results_df)]
  
  aggregated <- results_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) |>
    dplyr::summarize(
      ARI_mean = mean(ARI, na.rm = TRUE),
      ARI_sd = sd(ARI, na.rm = TRUE),
      VAF_mean = mean(VAF, na.rm = TRUE),
      VAF_sd = sd(VAF, na.rm = TRUE),
      BIC_mean = mean(BIC, na.rm = TRUE),
      BIC_sd = sd(BIC, na.rm = TRUE),
      subspace_similarity_mean = mean(subspace_similarity, na.rm = TRUE),
      subspace_similarity_sd = sd(subspace_similarity, na.rm = TRUE),
      runtime_mean = mean(runtime_sec, na.rm = TRUE),
      runtime_sd = sd(runtime_sec, na.rm = TRUE),
      K_correct_rate = mean(K_correct, na.rm = TRUE),
      r_correct_rate = mean(r_correct, na.rm = TRUE),
      Kr_correct_rate = mean(Kr_correct, na.rm = TRUE),
      n_reps = dplyr::n(),
      n_errors = sum(error_msg != "", na.rm = TRUE),
      .groups = "drop"
    )
  
  return(aggregated)
}


#' Compute Win Rates for Model Selection
#'
#' @param aggregated_df Aggregated results data frame.
#' @return Data frame with win rates.
#' @keywords internal
compute_modelsel_win_rates <- function(aggregated_df) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Please install it.")
  }
  
  # Define metrics and whether higher is better
  metrics <- list(
    ARI = list(col = "ARI_mean", maximize = TRUE),
    VAF = list(col = "VAF_mean", maximize = TRUE),
    subspace_similarity = list(col = "subspace_similarity_mean", maximize = TRUE),
    Kr_accuracy = list(col = "Kr_correct_rate", maximize = TRUE),
    BIC = list(col = "BIC_mean", maximize = FALSE),
    runtime = list(col = "runtime_mean", maximize = FALSE)
  )
  
  # Grouping variables for conditions
  condition_vars <- c("N", "K_true", "r_true", "M", "T_each",
                      "cluster_sep_spatial", "cluster_sep_temporal",
                      "cluster_sep_stability", "mode", "r_mode")
  condition_vars <- condition_vars[condition_vars %in% names(aggregated_df)]
  
  win_rates_list <- list()
  
  for (metric_name in names(metrics)) {
    metric_info <- metrics[[metric_name]]
    col_name <- metric_info$col
    maximize <- metric_info$maximize
    
    if (!col_name %in% names(aggregated_df)) next
    
    # Find winner for each condition
    if (maximize) {
      winners <- aggregated_df |>
        dplyr::group_by(dplyr::across(dplyr::all_of(condition_vars))) |>
        dplyr::filter(!!rlang::sym(col_name) == max(!!rlang::sym(col_name), na.rm = TRUE)) |>
        dplyr::ungroup()
    } else {
      winners <- aggregated_df |>
        dplyr::group_by(dplyr::across(dplyr::all_of(condition_vars))) |>
        dplyr::filter(!!rlang::sym(col_name) == min(!!rlang::sym(col_name), na.rm = TRUE)) |>
        dplyr::ungroup()
    }
    
    # Count wins
    win_counts <- winners |>
      dplyr::group_by(method, selection_strategy) |>
      dplyr::summarize(n_wins = dplyr::n(), .groups = "drop")

    # Total conditions
    total_conditions <- aggregated_df |>
      dplyr::select(dplyr::all_of(condition_vars)) |>
      dplyr::distinct() |>
      nrow()
    
    win_counts$win_rate <- win_counts$n_wins / total_conditions * 100
    win_counts$metric <- metric_name
    win_counts$total_conditions <- total_conditions
    
    win_rates_list[[metric_name]] <- win_counts
  }
  
  win_rates <- dplyr::bind_rows(win_rates_list)
  win_rates <- win_rates |>
    dplyr::select(metric, method, selection_strategy, n_wins, total_conditions, win_rate) |>
    dplyr::arrange(metric, dplyr::desc(win_rate))
  
  return(win_rates)
}
