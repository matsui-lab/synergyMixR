# PBS Parallelization Utilities
#
# Functions for submitting parameter sweeps to PBS/OpenPBS job schedulers.
# Supports array jobs for efficient multi-node parallelization on HPC clusters.

#' Run a Single Task from Parameter Grid
#'
#' Executes a single parameter combination from the parameter grid.
#' This function is called by PBS array job workers.
#'
#' @param task_id Integer. The task ID corresponding to run_id in the parameter grid.
#' @param param_grid_file Character. Path to the RDS file containing the parameter grid.
#' @param output_dir Character. Directory to save results.
#' @param sweep_type Character. Type of sweep: "modelsel" or "baselines".
#' @param methods Character vector. Methods to evaluate (e.g., c("MFA", "MPCA")).
#' @param strategies Character vector. Selection strategies (e.g., c("grid", "twostage")).
#' @param K_candidates Integer vector. Candidate K values for model selection.
#' @param r_candidates Integer vector. Candidate r values for model selection.
#' @param mode Character. Selection mode: "modelsel", "oracle_r", or "oracle_Kr".
#' @param criterion Character. Model selection criterion: "ICL" or "BIC".
#' @param icl_temperature Numeric. Temperature parameter for ICL entropy scaling.
#' @param max_iter Integer. Maximum EM iterations.
#' @param tol Numeric. Convergence tolerance.
#' @param n_init Integer. Number of random initializations.
#' @param use_kmeans_init Logical. Whether to use k-means initialization.
#' @param align_basis Logical. For baselines sweep: whether to align basis vectors.
#' @param refine_assignments Logical. For baselines sweep: whether to refine cluster assignments.
#' @param verbose Logical. Whether to print progress.
#'
#' @return Invisibly returns TRUE on success.
#'
#' @details
#' This function reads the parameter grid from a file, extracts the row
#' corresponding to the given task_id, runs the computation, and saves
#' the result to a separate file. This enables PBS array jobs where each
#' task processes one parameter combination independently.
#'
#' @examples
#' \dontrun{
#' # Called from PBS script:
#' # Rscript -e "synergyMixR::run_single_task(
#' #   task_id = as.integer(Sys.getenv('PBS_ARRAY_INDEX')),
#' #   param_grid_file = 'results/param_grid.rds',
#' #   output_dir = 'results'
#' # )"
#' }
#'
#' @export
run_single_task <- function(
    task_id,
    param_grid_file,
    output_dir,
    sweep_type = c("modelsel", "baselines"),
    methods = c("MFA", "MPCA"),
    strategies = c("grid", "twostage"),
    K_candidates = 1:5,
    r_candidates = 1:5,
    mode = c("modelsel", "oracle_r", "oracle_Kr"),
    criterion = c("ICL", "BIC"),
    icl_temperature = 2.0,
    max_iter = 50,
    tol = 1e-3,
    n_init = 1,
    use_kmeans_init = TRUE,
    align_basis = TRUE,
    refine_assignments = TRUE,
    verbose = TRUE
) {
  sweep_type <- match.arg(sweep_type)
  mode <- match.arg(mode)
  criterion <- match.arg(criterion)

  # Validate inputs
  stopifnot(
    "task_id must be a positive integer" = is.numeric(task_id) && task_id > 0,
    "param_grid_file must exist" = file.exists(param_grid_file)
  )

  # Read parameter grid
  param_grid <- readRDS(param_grid_file)

  # Find the row with matching run_id
 p <- param_grid[param_grid$run_id == task_id, ]

  if (nrow(p) == 0) {
    stop("task_id ", task_id, " not found in param_grid")
  }
  if (nrow(p) > 1) {
    warning("Multiple rows found for task_id ", task_id, ", using first row")
    p <- p[1, ]
  }

  if (verbose) {
    message("=== Task ", task_id, " ===")
    message("Parameters: N=", p$N, ", K=", p$K, ", r=", p$r,
            ", seed=", p$seed)
  }

  # Create partial results directory
  partial_dir <- file.path(output_dir, "partial_results")
  dir.create(partial_dir, recursive = TRUE, showWarnings = FALSE)

  # Get true K and r (handle both naming conventions)
  K_true <- if (!is.null(p$K_true)) p$K_true else p$K
  r_true <- if (!is.null(p$r_true)) p$r_true else p$r

  # Run the appropriate sweep type
  result <- tryCatch({
    if (sweep_type == "modelsel") {
      run_single_task_modelsel(
        p = p,
        K_true = K_true,
        r_true = r_true,
        methods = methods,
        strategies = strategies,
        K_candidates = K_candidates,
        r_candidates = r_candidates,
        mode = mode,
        criterion = criterion,
        icl_temperature = icl_temperature,
        max_iter = max_iter,
        tol = tol,
        n_init = n_init,
        use_kmeans_init = use_kmeans_init,
        verbose = verbose
      )
    } else {
      run_baselines_task_internal(
        p = p,
        K_true = K_true,
        r_true = r_true,
        align_basis = align_basis,
        refine_assignments = refine_assignments,
        verbose = verbose
      )
    }
  }, error = function(e) {
    warning("Error in task ", task_id, ": ", e$message)
    # Return error result
    data.frame(
      run_id = task_id,
      error_msg = as.character(e$message),
      stringsAsFactors = FALSE
    )
  })

  # Save result
  result_file <- file.path(partial_dir, sprintf("result_%06d.rds", task_id))
  saveRDS(result, result_file)

  if (verbose) {
    message("Result saved: ", result_file)
  }

  invisible(TRUE)
}


#' Run Single Task for Model Selection Sweep
#'
#' Internal function to run model selection for a single parameter combination.
#'
#' @param p Data frame row with parameter values.
#' @param K_true Integer. True number of clusters.
#' @param r_true Integer. True number of synergies.
#' @param methods Character vector. Methods to evaluate.
#' @param strategies Character vector. Selection strategies.
#' @param K_candidates Integer vector. Candidate K values.
#' @param r_candidates Integer vector. Candidate r values.
#' @param mode Character. Selection mode.
#' @param max_iter Integer. Maximum EM iterations.
#' @param tol Numeric. Convergence tolerance.
#' @param n_init Integer. Number of initializations.
#' @param use_kmeans_init Logical. Use k-means initialization.
#' @param verbose Logical. Print progress.
#'
#' @return Data frame with results for all method x strategy combinations.
#' @keywords internal
run_single_task_modelsel <- function(
    p,
    K_true,
    r_true,
    methods,
    strategies,
    K_candidates,
    r_candidates,
    mode,
    criterion = "ICL",
    icl_temperature = 2.0,
    max_iter,
    tol,
    n_init,
    use_kmeans_init,
    verbose
) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.")
  }

  results_list <- list()

  # Generate v2 synergy data
  data <- generate_synergy_data(
    N = p$N,
    K = K_true,
    r = r_true,
    M = p$M,
    T_each = p$T_each,
    cluster_sep_spatial = if (!is.null(p$cluster_sep_spatial)) p$cluster_sep_spatial else p$sep_spatial,
    cluster_sep_temporal = if (!is.null(p$cluster_sep_temporal)) p$cluster_sep_temporal else p$sep_temporal,
    cluster_sep_stability = if (!is.null(p$cluster_sep_stability)) p$cluster_sep_stability else p$sep_stability,
    seed = p$seed
  )

  list_of_data <- data$list_of_data
  true_cluster <- data$true_cluster
  true_Lambda_list <- data$Lambda_list

  # Determine effective candidates based on mode
  effective_K_candidates <- switch(mode,
    "oracle_Kr" = K_true,
    K_candidates
  )
  effective_r_candidates <- switch(mode,
    "oracle_r" = r_true,
    "oracle_Kr" = r_true,
    r_candidates
  )

  # Run model selection for each method x strategy combination
  for (method in methods) {
    for (strategy in strategies) {
      if (verbose) {
        message("  Running ", method, " with ", strategy, " strategy...")
      }

      result_row <- tryCatch({
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
          tol = tol,
          n_init = n_init,
          use_kmeans_init = use_kmeans_init,
          mc_cores_grid = 1,
          mc_cores = 1,
          n_threads = 1,
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
          pred_cluster <- extract_cluster_labels(best_model)
          if (!any(is.na(pred_cluster))) {
            alignment <- align_clusters_to_truth(true_cluster, pred_cluster)
            if (!is.na(alignment$accuracy)) {
              subspace_result <- compute_clusterwise_subspace_score(
                true_Lambda_list,
                loading_matrix,
                alignment$mapping
              )
              subspace_sim <- subspace_result$mean_score
            }
          }
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
          cluster_sep_spatial = if (!is.null(p$cluster_sep_spatial)) p$cluster_sep_spatial else p$sep_spatial,
          cluster_sep_temporal = if (!is.null(p$cluster_sep_temporal)) p$cluster_sep_temporal else p$sep_temporal,
          cluster_sep_stability = if (!is.null(p$cluster_sep_stability)) p$cluster_sep_stability else p$sep_stability,
          method = method,
          selection_strategy = strategy,
          mode = mode,
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
          runtime_sec = runtime_sec,
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
          cluster_sep_spatial = if (!is.null(p$cluster_sep_spatial)) p$cluster_sep_spatial else p$sep_spatial,
          cluster_sep_temporal = if (!is.null(p$cluster_sep_temporal)) p$cluster_sep_temporal else p$sep_temporal,
          cluster_sep_stability = if (!is.null(p$cluster_sep_stability)) p$cluster_sep_stability else p$sep_stability,
          method = method,
          selection_strategy = strategy,
          mode = mode,
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
}


#' Run Single Task for Baselines Sweep
#'
#' Internal function to run baseline methods for a single parameter combination.
#'
#' @param p Data frame row with parameter values.
#' @param K_true Integer. True number of clusters.
#' @param r_true Integer. True number of synergies.
#' @param align_basis Logical. Whether to align basis vectors.
#' @param refine_assignments Logical. Whether to refine cluster assignments.
#' @param verbose Logical. Print progress.
#'
#' @return Data frame with results.
#' @keywords internal
run_baselines_task_internal <- function(p, K_true, r_true, align_basis = TRUE,
                                        refine_assignments = TRUE, verbose = TRUE) {
  # Run comparison using compare_baselines (consistent with local mode)
  result <- compare_baselines(
    N = p$N,
    K = K_true,
    r = r_true,
    M = p$M,
    T_each = p$T_each,
    cluster_sep_spatial = if (!is.null(p$cluster_sep_spatial)) p$cluster_sep_spatial else p$sep_spatial,
    cluster_sep_temporal = if (!is.null(p$cluster_sep_temporal)) p$cluster_sep_temporal else p$sep_temporal,
    cluster_sep_stability = if (!is.null(p$cluster_sep_stability)) p$cluster_sep_stability else p$sep_stability,
    seed = p$seed,
    n_init = 1,
    use_kmeans_init = FALSE,
    mc_cores = 1,
    align_basis = align_basis,
    refine_assignments = refine_assignments
  )

  result$run_id <- p$run_id
  result
}


#' Generate PBS Job Script
#'
#' Generates a PBS array job script for submitting parameter sweep tasks.
#'
#' @param n_tasks Integer. Total number of tasks (rows in parameter grid).
#' @param output_dir Character. Directory for output files.
#' @param pbs_config List. PBS configuration options:
#'   \itemize{
#'     \item \code{job_name}: Job name (default: "param_sweep")
#'     \item \code{queue}: Queue name (default: NULL, use system default)
#'     \item \code{walltime}: Wall time limit (default: "04:00:00")
#'     \item \code{ncpus}: CPUs per task (default: 1)
#'     \item \code{mem}: Memory per task (default: "4gb")
#'     \item \code{r_module}: R module to load (default: NULL)
#'     \item \code{extra_modules}: Additional modules to load (default: NULL)
#'     \item \code{r_path}: Path to Rscript (default: "Rscript")
#'     \item \code{project}: Project/account code (default: NULL)
#'     \item \code{extra_pbs_options}: Additional PBS directives (default: NULL)
#'   }
#' @param sweep_type Character. Type of sweep: "modelsel" or "baselines".
#' @param sweep_config List. Sweep-specific configuration passed to run_single_task.
#'
#' @return Character. Path to the generated PBS script.
#'
#' @examples
#' \dontrun{
#' script_path <- generate_pbs_script(
#'   n_tasks = 100,
#'   output_dir = "results",
#'   pbs_config = list(
#'     queue = "batch",
#'     walltime = "02:00:00",
#'     r_module = "R/4.3.0"
#'   )
#' )
#' }
#'
#' @export
generate_pbs_script <- function(
    n_tasks,
    output_dir,
    pbs_config = list(),
    sweep_type = c("modelsel", "baselines"),
    sweep_config = list()
) {
  sweep_type <- match.arg(sweep_type)

  # Default PBS configuration
  defaults <- list(
    job_name = "param_sweep",
    queue = NULL,
    walltime = "04:00:00",
    ncpus = 1,
    mem = "4gb",
    r_module = NULL,
    extra_modules = NULL,
    r_path = "Rscript",
    project = NULL,
    extra_pbs_options = NULL
  )
  config <- modifyList(defaults, pbs_config)

  # Default sweep configuration based on sweep_type
  if (sweep_type == "baselines") {
    sweep_defaults <- list(
      align_basis = TRUE,
      refine_assignments = TRUE
    )
  } else {
    sweep_defaults <- list(
      methods = c("MFA", "MPCA"),
      strategies = c("grid", "twostage"),
      K_candidates = 1:5,
      r_candidates = 1:5,
      mode = "modelsel"
    )
  }
  sweep_cfg <- modifyList(sweep_defaults, sweep_config)

  # Build PBS directives
  pbs_directives <- c(
    sprintf("#PBS -N %s", config$job_name),
    sprintf("#PBS -J 1-%d", n_tasks),
    sprintf("#PBS -l select=1:ncpus=%d:mem=%s", config$ncpus, config$mem),
    sprintf("#PBS -l walltime=%s", config$walltime),
    "#PBS -j oe",  # Combine stdout and stderr
    sprintf("#PBS -o %s/logs/", output_dir)
  )

  if (!is.null(config$queue)) {
    pbs_directives <- c(pbs_directives, sprintf("#PBS -q %s", config$queue))
  }

  if (!is.null(config$project)) {
    pbs_directives <- c(pbs_directives, sprintf("#PBS -A %s", config$project))
  }

  if (!is.null(config$extra_pbs_options)) {
    pbs_directives <- c(pbs_directives, config$extra_pbs_options)
  }

  # Build module load commands
  module_commands <- ""
  if (!is.null(config$r_module)) {
    module_commands <- paste0(module_commands, "module load ", config$r_module, "\n")
  }
  if (!is.null(config$extra_modules)) {
    for (mod in config$extra_modules) {
      module_commands <- paste0(module_commands, "module load ", mod, "\n")
    }
  }

  # Build R command arguments based on sweep type
  if (sweep_type == "baselines") {
    r_args <- sprintf(
      'task_id = as.integer(Sys.getenv("PBS_ARRAY_INDEX")),
    param_grid_file = "%s/param_grid.rds",
    output_dir = "%s",
    sweep_type = "baselines",
    align_basis = %s,
    refine_assignments = %s',
      output_dir,
      output_dir,
      if (isTRUE(sweep_cfg$align_basis)) "TRUE" else "FALSE",
      if (isTRUE(sweep_cfg$refine_assignments)) "TRUE" else "FALSE"
    )
  } else {
    # Model selection sweep
    # Format K_candidates and r_candidates as valid R vectors
    k_cand_str <- paste0("c(", paste0(sweep_cfg$K_candidates, collapse = ", "), ")")
    r_cand_str <- paste0("c(", paste0(sweep_cfg$r_candidates, collapse = ", "), ")")

    # Get criterion and temperature with defaults
    criterion <- if (!is.null(sweep_cfg$criterion)) sweep_cfg$criterion else "ICL"
    icl_temperature <- if (!is.null(sweep_cfg$icl_temperature)) sweep_cfg$icl_temperature else 2.0

    r_args <- sprintf(
      'task_id = as.integer(Sys.getenv("PBS_ARRAY_INDEX")),
    param_grid_file = "%s/param_grid.rds",
    output_dir = "%s",
    sweep_type = "%s",
    methods = c(%s),
    strategies = c(%s),
    K_candidates = %s,
    r_candidates = %s,
    mode = "%s",
    criterion = "%s",
    icl_temperature = %s',
      output_dir,
      output_dir,
      sweep_type,
      paste0('"', sweep_cfg$methods, '"', collapse = ", "),
      paste0('"', sweep_cfg$strategies, '"', collapse = ", "),
      k_cand_str,
      r_cand_str,
      sweep_cfg$mode,
      criterion,
      icl_temperature
    )
  }

  # Build full script
  script <- paste0(
    "#!/bin/bash\n",
    paste(pbs_directives, collapse = "\n"), "\n\n",
    "# Change to submission directory\n",
    "cd $PBS_O_WORKDIR\n\n",
    "# Load modules\n",
    module_commands, "\n",
    "# Set environment variables\n",
    "export OMP_NUM_THREADS=1\n",
    "export OPENBLAS_NUM_THREADS=1\n",
    "export MKL_NUM_THREADS=1\n\n",
    "# Run single task\n",
    config$r_path, " --vanilla -e '\n",
    "library(synergyMixR)\n",
    "run_single_task(\n",
    r_args, "\n",
    ")\n",
    "'\n"
  )

  # Save script
  dir.create(file.path(output_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
  script_file <- file.path(output_dir, "submit_array.pbs")
  writeLines(script, script_file)

  message("PBS script generated: ", script_file)
  return(script_file)
}


#' Submit PBS Parameter Sweep
#'
#' Prepares and optionally submits a PBS array job for parameter sweep.
#'
#' @param param_grid Data frame. Parameter grid to sweep.
#' @param output_dir Character. Directory for output files.
#' @param pbs_config List. PBS configuration options (see \code{\link{generate_pbs_script}}).
#' @param sweep_type Character. Type of sweep: "modelsel" or "baselines".
#' @param sweep_config List. Sweep-specific configuration.
#' @param submit Logical. If TRUE, submit the job. If FALSE, only generate scripts.
#' @param verbose Logical. Print progress messages.
#'
#' @return List containing:
#'   \itemize{
#'     \item \code{job_id}: PBS job ID (if submitted)
#'     \item \code{script_file}: Path to the PBS script
#'     \item \code{param_grid_file}: Path to the saved parameter grid
#'     \item \code{output_dir}: Output directory
#'     \item \code{n_tasks}: Number of tasks
#'   }
#'
#' @examples
#' \dontrun{
#' # Generate scripts only (dry run)
#' result <- submit_pbs_sweep(
#'   param_grid = make_param_grid_test(),
#'   output_dir = "results",
#'   submit = FALSE
#' )
#'
#' # Submit job
#' result <- submit_pbs_sweep(
#'   param_grid = make_param_grid_paper_core(),
#'   output_dir = "results",
#'   pbs_config = list(queue = "batch"),
#'   submit = TRUE
#' )
#' }
#'
#' @export
submit_pbs_sweep <- function(
    param_grid,
    output_dir,
    pbs_config = list(),
    sweep_type = c("modelsel", "baselines"),
    sweep_config = list(),
    submit = FALSE,
    verbose = TRUE
) {
  sweep_type <- match.arg(sweep_type)

  # Create output directories
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "partial_results"), showWarnings = FALSE)
  dir.create(file.path(output_dir, "logs"), showWarnings = FALSE)

  n_tasks <- nrow(param_grid)

  if (verbose) {
    message("=== PBS Parameter Sweep Setup ===")
    message("Output directory: ", output_dir)
    message("Total tasks: ", n_tasks)
    message("Sweep type: ", sweep_type)
  }

  # Save parameter grid
  param_grid_file <- file.path(output_dir, "param_grid.rds")
  saveRDS(param_grid, param_grid_file)
  if (verbose) {
    message("Parameter grid saved: ", param_grid_file)
  }

  # Also save as CSV for inspection
  param_grid_csv <- file.path(output_dir, "param_grid.csv")
  write.csv(param_grid, param_grid_csv, row.names = FALSE)

  # Generate PBS script
  script_file <- generate_pbs_script(
    n_tasks = n_tasks,
    output_dir = output_dir,
    pbs_config = pbs_config,
    sweep_type = sweep_type,
    sweep_config = sweep_config
  )

  result <- list(
    job_id = NULL,
    script_file = script_file,
    param_grid_file = param_grid_file,
    output_dir = output_dir,
    n_tasks = n_tasks
  )

  if (!submit) {
    if (verbose) {
      message("\n=== Dry Run Mode ===")
      message("To submit manually:")
      message("  cd ", normalizePath(output_dir))
      message("  qsub submit_array.pbs")
      message("\nTo aggregate results after completion:")
      message("  aggregate_partial_results('", output_dir, "')")
    }
    return(invisible(result))
  }

  # Check if qsub is available
  qsub_check <- suppressWarnings(system2("which", "qsub", stdout = TRUE, stderr = TRUE))
  if (length(qsub_check) == 0 || !file.exists(qsub_check[1])) {
    warning("qsub command not found. Please submit manually:\n  qsub ", script_file)
    return(invisible(result))
  }

  # Submit job
  if (verbose) {
    message("\nSubmitting PBS job...")
  }

  submit_result <- system2(
    "qsub",
    script_file,
    stdout = TRUE,
    stderr = TRUE
  )

  if (length(submit_result) > 0) {
    # Extract job ID (format varies: "12345[]" or "12345.server[]")
    job_id <- gsub("\\[\\].*", "", submit_result[1])
    result$job_id <- job_id

    if (verbose) {
      message("Job submitted: ", job_id)
      message("\nMonitor with:")
      message("  qstat -t ", job_id)
      message("\nAfter completion, aggregate results with:")
      message("  aggregate_partial_results('", output_dir, "')")
    }
  } else {
    warning("Job submission may have failed. Check PBS logs.")
  }

  invisible(result)
}


#' Aggregate Partial Results from PBS Jobs
#'
#' Collects and combines results from individual PBS array job tasks.
#'
#' @param output_dir Character. Directory containing partial results.
#' @param pattern Character. Regex pattern for result files.
#' @param remove_partial Logical. If TRUE, remove partial result files after aggregation.
#' @param verbose Logical. Print progress messages.
#'
#' @return Data frame with combined results, or NULL if no results found.
#'
#' @examples
#' \dontrun{
#' # After PBS jobs complete
#' results <- aggregate_partial_results("results_sweep_modelsel")
#'
#' # View results
#' head(results)
#' }
#'
#' @export
aggregate_partial_results <- function(
    output_dir,
    pattern = "result_.*\\.rds$",
    remove_partial = FALSE,
    verbose = TRUE
) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.")
  }

  partial_dir <- file.path(output_dir, "partial_results")

  if (!dir.exists(partial_dir)) {
    warning("Partial results directory not found: ", partial_dir)
    return(NULL)
  }

  # Find result files
  result_files <- list.files(partial_dir, pattern = pattern, full.names = TRUE)

  if (length(result_files) == 0) {
    warning("No partial result files found in: ", partial_dir)
    return(NULL)
  }

  if (verbose) {
    message("Found ", length(result_files), " partial result files")
  }

  # Read all results
  results_list <- lapply(result_files, function(f) {
    tryCatch({
      readRDS(f)
    }, error = function(e) {
      warning("Failed to read: ", f, " - ", e$message)
      NULL
    })
  })

  # Remove NULLs
  results_list <- results_list[!sapply(results_list, is.null)]

  if (length(results_list) == 0) {
    warning("No valid results could be read")
    return(NULL)
  }

  if (verbose) {
    message("Successfully read ", length(results_list), " result files")
  }

  # Combine results
  all_results <- dplyr::bind_rows(results_list)

  if (verbose) {
    message("Total rows: ", nrow(all_results))
    n_errors <- sum(all_results$error_msg != "", na.rm = TRUE)
    if (n_errors > 0) {
      message("Rows with errors: ", n_errors)
    }
  }

  # Save combined results
  timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")

  raw_rds <- file.path(output_dir, paste0("sweep_raw_", timestamp_str, ".rds"))
  saveRDS(all_results, raw_rds)
  if (verbose) message("Raw results saved: ", raw_rds)

  raw_csv <- file.path(output_dir, paste0("sweep_raw_", timestamp_str, ".csv"))
  write.csv(all_results, raw_csv, row.names = FALSE)
  if (verbose) message("CSV saved: ", raw_csv)

  # Check if this is model selection results (has specific columns)
  if ("selection_strategy" %in% names(all_results)) {
    # Aggregate model selection results
    aggregated <- aggregate_modelsel_results(all_results)
  } else if ("Method" %in% names(all_results)) {
    # Aggregate baseline results
    aggregated <- aggregate_results(all_results)
  } else {
    aggregated <- all_results
  }

  agg_csv <- file.path(output_dir, "sweep_aggregated.csv")
  write.csv(aggregated, agg_csv, row.names = FALSE)
  if (verbose) message("Aggregated results saved: ", agg_csv)

  # Optionally remove partial files
 if (remove_partial) {
    if (verbose) message("Removing partial result files...")
    file.remove(result_files)
  }

  if (verbose) {
    message("\n=== Aggregation Complete ===")
    message("Results saved to: ", output_dir)
  }

  return(all_results)
}


#' Check PBS Job Status
#'
#' Checks the status of a PBS job.
#'
#' @param job_id Character. PBS job ID.
#'
#' @return List with job status information.
#'
#' @examples
#' \dontrun{
#' status <- check_pbs_job_status("12345")
#' }
#'
#' @export
check_pbs_job_status <- function(job_id) {
  if (is.null(job_id) || job_id == "") {
    return(list(running = FALSE, completed = FALSE, error = "No job ID provided"))
  }

  # Run qstat - try different formats for compatibility
  # Some PBS versions need -t for array jobs, others don't
  result <- suppressWarnings(
    system2("qstat", c("-t", job_id), stdout = TRUE, stderr = TRUE)
  )

  # If -t failed, try without it

  if (length(result) == 0 || any(grepl("illegal", result, ignore.case = TRUE))) {
    result <- suppressWarnings(
      system2("qstat", job_id, stdout = TRUE, stderr = TRUE)
    )
  }

  # Parse output
  if (length(result) == 0 || any(grepl("Unknown Job Id", result, ignore.case = TRUE))) {
    # Job not found - might be completed
    return(list(
      running = FALSE,
      completed = TRUE,
      job_id = job_id,
      message = "Job not in queue (possibly completed)"
    ))
  }

  # Count running/queued tasks
  # PBS Pro / OpenPBS format: job state in column (R=running, Q=queued, H=held, E=exiting, F=finished)
  # Match state codes at word boundaries or end of fields
  n_running <- sum(grepl("\\bR\\b|\\sR\\s|\\sR$", result))
  n_queued <- sum(grepl("\\bQ\\b|\\sQ\\s|\\sQ$", result))
  n_held <- sum(grepl("\\bH\\b|\\sH\\s|\\sH$", result))
  n_exiting <- sum(grepl("\\bE\\b|\\sE\\s|\\sE$", result))
  n_finished <- sum(grepl("\\bF\\b|\\sF\\s|\\sF$", result))

  # Job is still running if there are running, queued, or exiting tasks
  still_active <- (n_running + n_queued + n_exiting) > 0

  list(
    running = still_active,
    completed = !still_active,
    job_id = job_id,
    n_running = n_running,
    n_queued = n_queued,
    n_held = n_held,
    n_exiting = n_exiting,
    n_finished = n_finished,
    raw_output = result
  )
}


#' Wait for PBS Job Completion
#'
#' Blocks until a PBS job completes.
#'
#' @param job_id Character. PBS job ID.
#' @param check_interval Integer. Seconds between status checks.
#' @param timeout Integer. Maximum seconds to wait (NULL for no timeout).
#' @param initial_wait Integer. Seconds to wait before first status check.
#'   This gives time for the job to be registered in the queue.
#' @param output_dir Character. If provided, also monitors partial_results
#'   directory to track progress.
#' @param n_tasks Integer. Total number of expected tasks (used with output_dir).
#' @param verbose Logical. Print progress messages.
#'
#' @return Logical. TRUE if job completed, FALSE if timeout.
#'
#' @examples
#' \dontrun{
#' wait_for_pbs_job("12345", timeout = 3600)  # Wait up to 1 hour
#' }
#'
#' @export
wait_for_pbs_job <- function(
    job_id,
    check_interval = 60,
    timeout = NULL,
    initial_wait = 5,
    output_dir = NULL,
    n_tasks = NULL,
    verbose = TRUE
) {
  start_time <- Sys.time()

  # Initial wait for job to be registered in the queue
  if (initial_wait > 0) {
    if (verbose) {
      message(sprintf("[%s] Waiting %d seconds for job to be registered...",
                      format(Sys.time(), "%H:%M:%S"), initial_wait))
    }
    Sys.sleep(initial_wait)
  }

  while (TRUE) {
    status <- check_pbs_job_status(job_id)

    # Check partial results if output_dir is provided
    tasks_completed <- NA
    if (!is.null(output_dir)) {
      partial_dir <- file.path(output_dir, "partial_results")
      if (dir.exists(partial_dir)) {
        result_files <- list.files(partial_dir, pattern = "result_.*\\.rds$")
        tasks_completed <- length(result_files)
      }
    }

    if (!status$running) {
      if (verbose) message("Job completed: ", job_id)
      if (!is.na(tasks_completed)) {
        if (verbose) message("Tasks with results: ", tasks_completed)
      }
      return(TRUE)
    }

    if (verbose) {
      progress_msg <- sprintf(
        "[%s] Job %s: %d running, %d queued",
        format(Sys.time(), "%H:%M:%S"),
        job_id,
        status$n_running,
        status$n_queued
      )
      if (!is.null(n_tasks) && !is.na(tasks_completed)) {
        progress_msg <- paste0(progress_msg,
                               sprintf(" | Completed: %d/%d", tasks_completed, n_tasks))
      }
      message(progress_msg)
    }

    # Check timeout
    if (!is.null(timeout)) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      if (elapsed > timeout) {
        warning("Timeout waiting for job: ", job_id)
        return(FALSE)
      }
    }

    Sys.sleep(check_interval)
  }
}


#' Count Completed PBS Tasks
#'
#' Counts how many tasks have completed by checking partial result files.
#'
#' @param output_dir Character. Output directory.
#' @param n_total Integer. Total number of expected tasks (optional).
#'
#' @return List with completion statistics.
#'
#' @export
count_completed_tasks <- function(output_dir, n_total = NULL) {
  partial_dir <- file.path(output_dir, "partial_results")

  if (!dir.exists(partial_dir)) {
    return(list(n_completed = 0, n_total = n_total, pct_complete = 0))
  }

  result_files <- list.files(partial_dir, pattern = "result_.*\\.rds$")
  n_completed <- length(result_files)

  if (is.null(n_total)) {
    # Try to read from param_grid
    param_grid_file <- file.path(output_dir, "param_grid.rds")
    if (file.exists(param_grid_file)) {
      param_grid <- readRDS(param_grid_file)
      n_total <- nrow(param_grid)
    }
  }

  pct_complete <- if (!is.null(n_total) && n_total > 0) {
    round(100 * n_completed / n_total, 1)
  } else {
    NA
  }

  list(
    n_completed = n_completed,
    n_total = n_total,
    pct_complete = pct_complete
  )
}


# =============================================================================
# Future-based Parallelization Support
# =============================================================================

#' Setup Future Backend for HPC
#'
#' Configures the future backend for parallel execution. Supports local
#' multicore, PBS/TORQUE, SLURM, and other HPC schedulers via future.batchtools.
#'
#' @param backend Character. Backend type:
#'   \itemize{
#'     \item "multisession": Local parallel sessions (default)
#'     \item "multicore": Local forked processes (Unix only)
#'     \item "cluster": Manual cluster setup
#'     \item "batchtools_torque": PBS/TORQUE scheduler
#'     \item "batchtools_slurm": SLURM scheduler
#'     \item "batchtools_sge": SGE scheduler
#'     \item "batchtools_lsf": LSF scheduler
#'   }
#' @param workers Integer. Number of workers for local backends.
#' @param resources List. HPC resource specifications for batchtools backends:
#'   \itemize{
#'     \item \code{walltime}: Wall time limit (e.g., "04:00:00")
#'     \item \code{memory}: Memory per job (e.g., "4gb" or "4096")
#'     \item \code{ncpus}: CPUs per job
#'     \item \code{queue}: Queue/partition name
#'   }
#' @param template Character. Path to batchtools template file (optional).
#' @param verbose Logical. Print setup information.
#'
#' @return Invisibly returns the previous plan.
#'
#' @details
#' This function is a wrapper around \code{future::plan()} that simplifies
#' HPC backend configuration. For PBS/TORQUE, it uses
#' \code{future.batchtools::batchtools_torque}.
#'
#' @examples
#' \dontrun{
#' # Local multicore
#' setup_future_backend("multisession", workers = 4)
#'
#' # PBS/TORQUE with resources
#' setup_future_backend(
#'   "batchtools_torque",
#'   resources = list(walltime = "02:00:00", memory = "8gb")
#' )
#'
#' # Reset to sequential
#' setup_future_backend("sequential")
#' }
#'
#' @export
setup_future_backend <- function(
    backend = c("multisession", "multicore", "cluster", "sequential",
                "batchtools_torque", "batchtools_slurm",
                "batchtools_sge", "batchtools_lsf"),
    workers = NULL,
    resources = list(),
    template = NULL,
    verbose = TRUE
) {
  backend <- match.arg(backend)

  # Check for required packages
  if (!requireNamespace("future", quietly = TRUE)) {
    stop("Package 'future' is required. Install with: install.packages('future')")
  }

  is_batchtools <- grepl("^batchtools_", backend)

  if (is_batchtools) {
    if (!requireNamespace("future.batchtools", quietly = TRUE)) {
      stop("Package 'future.batchtools' is required for HPC backends. ",
           "Install with: install.packages('future.batchtools')")
    }
  }

  # Build plan arguments
  plan_args <- list()

  # Set workers for local backends
  if (!is.null(workers) && backend %in% c("multisession", "multicore", "cluster")) {
    plan_args$workers <- workers
  }

  # Set resources for batchtools backends
  if (is_batchtools && length(resources) > 0) {
    plan_args$resources <- resources
  }

  # Set template for batchtools backends
  if (is_batchtools && !is.null(template)) {
    plan_args$template <- template
  }

  # Get the plan function
  plan_func <- switch(backend,
    "sequential" = future::sequential,
    "multisession" = future::multisession,
    "multicore" = future::multicore,
    "cluster" = future::cluster,
    "batchtools_torque" = future.batchtools::batchtools_torque,
    "batchtools_slurm" = future.batchtools::batchtools_slurm,
    "batchtools_sge" = future.batchtools::batchtools_sge,
    "batchtools_lsf" = future.batchtools::batchtools_lsf
  )

  # Set the plan
  if (verbose) {
    message("Setting up future backend: ", backend)
    if (length(plan_args) > 0) {
      message("  Options: ", paste(names(plan_args), "=",
                                   sapply(plan_args, function(x) {
                                     if (is.list(x)) paste(names(x), collapse = ",")
                                     else as.character(x)
                                   }), collapse = ", "))
    }
  }

  old_plan <- do.call(future::plan, c(list(plan_func), plan_args))

  if (verbose) {
    message("Future backend ready.")
  }

  invisible(old_plan)
}


#' Run Parameter Sweep Using Future
#'
#' Executes a parameter sweep using the future framework for parallelization.
#' This works with any future backend (local, PBS, SLURM, etc.).
#'
#' @param param_grid Data frame. Parameter grid to sweep.
#' @param run_func Function. Function to run for each parameter combination.
#'   Should accept a single row of param_grid and return a data frame.
#' @param ... Additional arguments passed to run_func.
#' @param future_backend Character. Backend to use (see \code{\link{setup_future_backend}}).
#' @param workers Integer. Number of workers for local backends.
#' @param resources List. HPC resources for batchtools backends.
#' @param chunk_size Integer. Number of tasks per future (for load balancing).
#' @param verbose Logical. Print progress information.
#'
#' @return Data frame with combined results from all parameter combinations.
#'
#' @examples
#' \dontrun{
#' # Define a simple run function
#' my_run_func <- function(p) {
#'   Sys.sleep(1)  # Simulate work
#'   data.frame(run_id = p$run_id, result = p$x^2)
#' }
#'
#' # Create parameter grid
#' param_grid <- data.frame(run_id = 1:10, x = 1:10)
#'
#' # Run with local parallelization
#' results <- run_sweep_future(
#'   param_grid, my_run_func,
#'   future_backend = "multisession",
#'   workers = 2
#' )
#'
#' # Run with PBS
#' results <- run_sweep_future(
#'   param_grid, my_run_func,
#'   future_backend = "batchtools_torque",
#'   resources = list(walltime = "01:00:00")
#' )
#' }
#'
#' @export
run_sweep_future <- function(
    param_grid,
    run_func,
    ...,
    future_backend = "multisession",
    workers = NULL,
    resources = list(),
    chunk_size = NULL,
    verbose = TRUE
) {
  # Check for required packages
  if (!requireNamespace("future", quietly = TRUE)) {
    stop("Package 'future' is required.")
  }
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Package 'future.apply' is required. ",
         "Install with: install.packages('future.apply')")
  }

  n_tasks <- nrow(param_grid)

  if (verbose) {
    message("=== Future-based Parameter Sweep ===")
    message("Total tasks: ", n_tasks)
    message("Backend: ", future_backend)
  }

  # Setup backend
  old_plan <- setup_future_backend(
    backend = future_backend,
    workers = workers,
    resources = resources,
    verbose = verbose
  )

  # Ensure we restore the old plan on exit
  on.exit(future::plan(old_plan), add = TRUE)

  # Prepare run function wrapper
  run_wrapper <- function(i) {
    p <- param_grid[i, , drop = FALSE]
    tryCatch({
      run_func(p, ...)
    }, error = function(e) {
      data.frame(
        run_id = if ("run_id" %in% names(p)) p$run_id else i,
        error_msg = as.character(e$message),
        stringsAsFactors = FALSE
      )
    })
  }

  # Run with future_lapply
  if (verbose) message("Running tasks...")
  start_time <- Sys.time()

  # Set chunk size for better load balancing with HPC
  future_options <- list()
  if (!is.null(chunk_size)) {
    future_options$future.chunk.size <- chunk_size
  }

  results_list <- do.call(
    future.apply::future_lapply,
    c(
      list(
        X = seq_len(n_tasks),
        FUN = run_wrapper,
        future.seed = TRUE  # Ensure reproducibility
      ),
      future_options
    )
  )

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

  if (verbose) {
    message("Completed in ", round(elapsed, 1), " minutes")
    message("Combining results...")
  }

  # Combine results
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    all_results <- do.call(rbind, results_list)
  } else {
    all_results <- dplyr::bind_rows(results_list)
  }

  if (verbose) {
    message("Total rows: ", nrow(all_results))
    if ("error_msg" %in% names(all_results)) {
      n_errors <- sum(all_results$error_msg != "", na.rm = TRUE)
      if (n_errors > 0) {
        message("Tasks with errors: ", n_errors)
      }
    }
  }

  all_results
}


#' Check if Future Packages are Available
#'
#' Checks whether the future ecosystem packages are installed.
#'
#' @param batchtools Logical. Also check for future.batchtools.
#'
#' @return Logical. TRUE if packages are available.
#'
#' @export
check_future_available <- function(batchtools = FALSE) {

  has_future <- requireNamespace("future", quietly = TRUE)
  has_apply <- requireNamespace("future.apply", quietly = TRUE)

  if (!batchtools) {
    return(has_future && has_apply)
  }

  has_batchtools <- requireNamespace("future.batchtools", quietly = TRUE)
  has_future && has_apply && has_batchtools
}


# =============================================================================
# Convenience Wrapper for Simple HPC Usage
# =============================================================================

#' Submit Parameter Sweep to HPC (Simple Interface)
#'
#' A simplified wrapper for submitting parameter sweeps to HPC clusters.
#' Automatically selects the appropriate parallelization strategy based on
#' the specified cluster type.
#'
#' @param cluster_type Character. Type of HPC cluster:
#'   \itemize{
#'     \item "local": Local multicore execution (no HPC)
#'     \item "pbs": PBS/OpenPBS/TORQUE cluster
#'     \item "slurm": SLURM cluster
#'     \item "sge": SGE cluster
#'   }
#' @param queue Character. Queue/partition name for the cluster. Required for
#'   HPC clusters.
#' @param walltime Character. Maximum wall time per job (default: "04:00:00").
#' @param memory Character. Memory per job (default: "4gb" for PBS, "4096" for SLURM).
#' @param test Logical. If TRUE, use test parameter grid. Default is FALSE.
#' @param cores Integer. Number of cores for local execution. Default is 4.
#' @param output_dir Character. Output directory for results.
#' @param submit Logical. If TRUE, submit jobs. If FALSE, just generate scripts.
#'   Default is TRUE.
#' @param verbose Logical. Print progress messages. Default is TRUE.
#'
#' @return Invisibly returns the result from \code{run_param_sweep_modelsel}.
#'
#' @details
#' This function provides a simplified interface for the most common HPC use cases.
#' It automatically:
#' \itemize{
#'   \item Selects the appropriate future backend for the cluster type
#'   \item Configures resource specifications
#'   \item Handles PBS vs SLURM differences
#' }
#'
#' For more control over the submission, use \code{\link{run_param_sweep_modelsel}}
#' directly with \code{parallel_strategy = "future"} or \code{"pbs"}.
#'
#' @examples
#' \dontrun{
#' # Local execution (8 cores)
#' submit_hpc_sweep("local", cores = 8, test = TRUE)
#'
#' # PBS cluster
#' submit_hpc_sweep("pbs", queue = "batch", walltime = "02:00:00")
#'
#' # SLURM cluster
#' submit_hpc_sweep("slurm", queue = "compute", memory = "8192")
#'
#' # Dry run (generate scripts without submitting)
#' submit_hpc_sweep("pbs", queue = "batch", submit = FALSE)
#' }
#'
#' @seealso \code{\link{run_param_sweep_modelsel}}, \code{\link{setup_future_backend}}
#'
#' @export
submit_hpc_sweep <- function(
    cluster_type = c("local", "pbs", "slurm", "sge"),
    queue = NULL,
    walltime = "04:00:00",
    memory = NULL,
    test = FALSE,
    cores = 4,
    output_dir = NULL,
    submit = TRUE,
    verbose = TRUE
) {
  cluster_type <- match.arg(cluster_type)

  # Set default memory based on cluster type
  if (is.null(memory)) {
    memory <- switch(cluster_type,
      "slurm" = "4096",  # SLURM uses MB
      "4gb"  # PBS/SGE use gb notation
    )
  }

  # Validate queue for HPC clusters
  if (cluster_type != "local" && is.null(queue)) {
    warning("No queue specified for HPC cluster. Using system default.")
  }

  if (verbose) {
    message("=== HPC Parameter Sweep Submission ===")
    message("Cluster type: ", cluster_type)
    if (cluster_type != "local") {
      message("Queue: ", if (is.null(queue)) "(default)" else queue)
      message("Walltime: ", walltime)
      message("Memory: ", memory)
    } else {
      message("Cores: ", cores)
    }
    message("Test mode: ", test)
    message("")
  }

  # Route to appropriate execution method
  if (cluster_type == "local") {
    # Local execution
    result <- run_param_sweep_modelsel(
      test = test,
      cores = cores,
      output_dir = output_dir,
      parallel_strategy = "grid",
      verbose = verbose
    )
  } else {
    # HPC execution via future.batchtools
    backend <- switch(cluster_type,
      "pbs" = "batchtools_torque",
      "slurm" = "batchtools_slurm",
      "sge" = "batchtools_sge"
    )

    # Build resources list
    resources <- list(
      walltime = walltime,
      memory = memory
    )
    if (!is.null(queue)) {
      # Different field names for different schedulers
      if (cluster_type == "slurm") {
        resources$partition <- queue
      } else {
        resources$queue <- queue
      }
    }

    # Check for required packages
    if (!check_future_available(batchtools = TRUE)) {
      stop("Required packages not installed. Run:\n",
           "  install.packages(c('future', 'future.apply', 'future.batchtools'))")
    }

    result <- run_param_sweep_modelsel(
      test = test,
      cores = cores,
      output_dir = output_dir,
      parallel_strategy = "future",
      future_backend = backend,
      future_resources = resources,
      verbose = verbose
    )
  }

  invisible(result)
}


#' Print Summary of Available Parallelization Options
#'
#' Displays a summary of available parallelization strategies and
#' their current availability status.
#'
#' @return Invisibly returns a list with availability status.
#'
#' @examples
#' show_parallel_options()
#'
#' @export
show_parallel_options <- function() {
  cat("\n")
  cat("=== Available Parallelization Options ===\n")
  cat("\n")

  # Check package availability
  has_parallel <- requireNamespace("parallel", quietly = TRUE)
  has_future <- requireNamespace("future", quietly = TRUE)
  has_apply <- requireNamespace("future.apply", quietly = TRUE)
  has_batchtools <- requireNamespace("future.batchtools", quietly = TRUE)

  # Local options
  cat("Local Execution:\n")
  cat(sprintf("  %-20s %s\n", "grid (mclapply)",
              if (has_parallel) "[OK]" else "[Not available]"))
  cat(sprintf("  %-20s %s\n", "sequential",
              "[OK] (always available)"))
  cat(sprintf("  %-20s %s\n", "future:multisession",
              if (has_future && has_apply) "[OK]" else "[Need: future, future.apply]"))
  cat(sprintf("  %-20s %s\n", "future:multicore",
              if (has_future && has_apply) "[OK] (Unix only)" else "[Need: future, future.apply]"))
  cat("\n")

  # HPC options
  cat("HPC Execution (requires future.batchtools):\n")
  hpc_status <- if (has_batchtools) "[OK]" else "[Need: future.batchtools]"
  cat(sprintf("  %-20s %s\n", "PBS/TORQUE", hpc_status))
  cat(sprintf("  %-20s %s\n", "SLURM", hpc_status))
  cat(sprintf("  %-20s %s\n", "SGE", hpc_status))
  cat(sprintf("  %-20s %s\n", "LSF", hpc_status))
  cat("\n")

  # Direct PBS submission
  qsub_available <- suppressWarnings(
    length(system2("which", "qsub", stdout = TRUE, stderr = TRUE)) > 0
  )
  cat("Direct PBS Submission (qsub):\n")
  cat(sprintf("  %-20s %s\n", "pbs",
              if (qsub_available) "[OK] (qsub found)" else "[qsub not in PATH]"))
  cat("\n")

  # Installation hints
  if (!has_future || !has_apply) {
    cat("To enable future-based parallelization:\n")
    cat("  install.packages(c('future', 'future.apply'))\n\n")
  }
  if (!has_batchtools) {
    cat("To enable HPC cluster support:\n")
    cat("  install.packages('future.batchtools')\n\n")
  }

  invisible(list(
    parallel = has_parallel,
    future = has_future,
    future.apply = has_apply,
    future.batchtools = has_batchtools,
    qsub = qsub_available
  ))
}
