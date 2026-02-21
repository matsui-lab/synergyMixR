#' Run K Selection Method Comparison (v3)
#'
#' Compares multiple K selection methods (BIC, ICL with various temperatures,
#' Ensemble, Adaptive, Silhouette) across a parameter grid.
#'
#' @param test Logical. If TRUE, use smaller test grid. Default FALSE.
#' @param cores Integer. Number of CPU cores. Default 4.
#' @param output_dir Character. Output directory. Default auto-detected.
#' @param parallel_strategy Character. "local", "pbs", or "future".
#' @param pbs_config List. PBS configuration (queue, walltime, mem, r_module).
#' @param pbs_submit Logical. If TRUE, submit PBS job. Default FALSE (dry run).
#' @param pbs_wait Logical. If TRUE, wait for PBS completion. Default FALSE.
#' @param Kvec Integer vector. Candidate K values. Default 2:6.
#' @param r_fixed Integer. Fixed r value for two-stage. Default 2.
#' @param max_iter Integer. Maximum EM iterations. Default 50.
#' @param verbose Logical. Print progress. Default TRUE.
#'
#' @return Data frame with results, or PBS job info if using PBS.
#'
#' @details
#' This function compares the following K selection methods:
#' \itemize{
#'   \item BIC: Bayesian Information Criterion
#'   \item ICL_T1: ICL with temperature = 1.0
#'   \item ICL_T2: ICL with temperature = 2.0 (recommended)
#'   \item ICL_T3: ICL with temperature = 3.0
#'   \item Ensemble: Majority vote of BIC, ICL_T1, ICL_T2
#'   \item Adaptive: Adaptive temperature selection
#'   \item Silhouette: Silhouette-based selection
#' }
#'
#' @examples
#' \dontrun{
#' # Local test run
#' results <- run_k_selection_comparison(test = TRUE, cores = 2)
#'
#' # Full simulation
#' results <- run_k_selection_comparison(test = FALSE, cores = 8)
#'
#' # PBS submission
#' run_k_selection_comparison(
#'   test = FALSE,
#'   parallel_strategy = "pbs",
#'   pbs_config = list(queue = "batch"),
#'   pbs_submit = TRUE
#' )
#' }
#'
#' @export
run_k_selection_comparison <- function(
    test = FALSE,
    cores = 4,
    output_dir = NULL,
    parallel_strategy = c("local", "pbs", "future"),
    pbs_config = list(),
    pbs_submit = FALSE,
    pbs_wait = FALSE,
    Kvec = 2:6,
    r_fixed = 2,
    max_iter = 50,
    verbose = TRUE
) {
  parallel_strategy <- match.arg(parallel_strategy)

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.")
  }

  # Set output directory
 if (is.null(output_dir)) {
    if (dir.exists("inst/paper")) {
      output_dir <- file.path("inst", "paper", "v3", "results_k_selection")
    } else {
      output_dir <- file.path(getwd(), "results_k_selection")
    }
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (verbose) {
    cat("\n")
    cat("=================================================\n")
    cat("  K Selection Method Comparison (v3)\n")
    cat("=================================================\n\n")
    cat("Output directory:", output_dir, "\n")
    cat("Mode:", if (test) "TEST" else "FULL", "\n")
    cat("Parallel strategy:", parallel_strategy, "\n\n")
  }

  # Define parameter grid
  if (test) {
    param_grid <- expand.grid(
      N = 40,
      K_true = c(2, 3),
      r_true = 2,
      M = 8,
      T_each = 100,
      cluster_sep_spatial = c(0.6, 1.0),
      cluster_sep_temporal = 0.6,
      cluster_sep_stability = 0.6,
      seed = 1:2,
      stringsAsFactors = FALSE
    )
  } else {
    param_grid <- expand.grid(
      N = c(40, 60),
      K_true = c(2, 3, 4, 5),
      r_true = 2,
      M = 8,
      T_each = 100,
      cluster_sep_spatial = c(0.3, 0.6, 1.0),
      cluster_sep_temporal = 0.6,
      cluster_sep_stability = 0.6,
      seed = 1:5,
      stringsAsFactors = FALSE
    )
  }
  param_grid$run_id <- seq_len(nrow(param_grid))

  if (verbose) {
    cat("Parameter combinations:", nrow(param_grid), "\n\n")
  }

  # K selection methods
  k_methods <- list(
    BIC = list(criterion = "BIC", icl_temperature = 1.0),
    ICL_T1 = list(criterion = "ICL", icl_temperature = 1.0),
    ICL_T2 = list(criterion = "ICL", icl_temperature = 2.0),
    ICL_T3 = list(criterion = "ICL", icl_temperature = 3.0)
  )

  # =========================================================================
  # PBS MODE
  # =========================================================================
  if (parallel_strategy == "pbs") {
    if (verbose) cat("=== PBS Mode ===\n")

    # Create directories
    dir.create(file.path(output_dir, "partial_results"), showWarnings = FALSE)
    dir.create(file.path(output_dir, "logs"), showWarnings = FALSE)

    # Save parameter grid
    param_grid_file <- file.path(output_dir, "param_grid.rds")
    saveRDS(param_grid, param_grid_file)
    write.csv(param_grid, file.path(output_dir, "param_grid.csv"), row.names = FALSE)

    if (verbose) cat("Parameter grid saved:", param_grid_file, "\n")

    # PBS defaults
    pbs_defaults <- list(
      queue = NULL,
      walltime = "04:00:00",
      mem = "4gb",
      r_module = NULL
    )
    pbs_cfg <- modifyList(pbs_defaults, pbs_config)

    # Generate PBS script
    n_tasks <- nrow(param_grid)
    pbs_directives <- c(
      "#PBS -N k_selection_v3",
      sprintf("#PBS -J 1-%d", n_tasks),
      sprintf("#PBS -l select=1:ncpus=1:mem=%s", pbs_cfg$mem),
      sprintf("#PBS -l walltime=%s", pbs_cfg$walltime),
      "#PBS -j oe",
      sprintf("#PBS -o %s/logs/", output_dir)
    )

    if (!is.null(pbs_cfg$queue)) {
      pbs_directives <- c(pbs_directives, sprintf("#PBS -q %s", pbs_cfg$queue))
    }

    module_cmd <- if (!is.null(pbs_cfg$r_module)) {
      paste0("module load ", pbs_cfg$r_module, "\n")
    } else ""

    script_path <- "inst/paper/v3/run_k_selection_comparison.R"
    pbs_script <- paste0(
      "#!/bin/bash\n",
      paste(pbs_directives, collapse = "\n"), "\n\n",
      "cd $PBS_O_WORKDIR\n\n",
      module_cmd, "\n",
      "export OMP_NUM_THREADS=1\n",
      "export OPENBLAS_NUM_THREADS=1\n",
      "export MKL_NUM_THREADS=1\n\n",
      "Rscript ", script_path, " --task-id $PBS_ARRAY_INDEX --output ", output_dir, "\n"
    )

    script_file <- file.path(output_dir, "submit_array.pbs")
    writeLines(pbs_script, script_file)
    if (verbose) cat("PBS script generated:", script_file, "\n\n")

    if (!pbs_submit) {
      if (verbose) {
        cat("=== Dry Run Mode ===\n")
        cat("To submit manually:\n")
        cat("  cd", normalizePath(output_dir), "\n")
        cat("  qsub submit_array.pbs\n\n")
        cat("After completion:\n")
        cat("  aggregate_partial_results('", output_dir, "')\n")
      }
      return(invisible(list(
        output_dir = output_dir,
        script_file = script_file,
        n_tasks = n_tasks
      )))
    }

    # Submit
    if (verbose) cat("Submitting PBS job...\n")
    submit_result <- system2("qsub", script_file, stdout = TRUE, stderr = TRUE)

    if (length(submit_result) > 0) {
      job_id <- gsub("\\[\\].*", "", submit_result[1])
      if (verbose) cat("Job submitted:", job_id, "\n")

      if (pbs_wait) {
        wait_for_pbs_job(job_id, output_dir = output_dir, n_tasks = n_tasks)
        return(aggregate_partial_results(output_dir))
      }

      return(invisible(list(
        job_id = job_id,
        output_dir = output_dir,
        n_tasks = n_tasks
      )))
    }

    return(invisible(NULL))
  }

  # =========================================================================
  # LOCAL MODE
  # =========================================================================
  run_single <- function(p) {
    results_list <- list()

    set.seed(p$seed)
    sim <- tryCatch({
      simulate_dynamic_synergy_data(
        N = p$N, K = p$K_true, r = p$r_true, M = p$M, T_each = p$T_each,
        cluster_sep_spatial = p$cluster_sep_spatial,
        cluster_sep_temporal = p$cluster_sep_temporal,
        cluster_sep_stability = p$cluster_sep_stability
      )
    }, error = function(e) NULL)

    if (is.null(sim)) return(NULL)

    list_of_data <- sim$X_list

    for (method_name in names(k_methods)) {
      cfg <- k_methods[[method_name]]
      start_time <- Sys.time()

      K_hat <- tryCatch({
        suppressMessages({
          res <- select_optimal_K_r_mfa_twostage(
            list_of_data = list_of_data,
            Kvec = Kvec, rvec = 1:3, r_fixed = r_fixed,
            criterion = cfg$criterion,
            icl_max_iter = 10,
            icl_temperature = cfg$icl_temperature,
            max_iter = max_iter,
            n_init = 1, use_kmeans_init = TRUE
          )
          res$best_K
        })
      }, error = function(e) NA)

      runtime_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

      results_list[[length(results_list) + 1]] <- data.frame(
        run_id = p$run_id, seed = p$seed, N = p$N,
        K_true = p$K_true, r_true = p$r_true,
        cluster_sep_spatial = p$cluster_sep_spatial,
        method = method_name,
        K_hat = K_hat,
        K_correct = (!is.na(K_hat) && K_hat == p$K_true),
        runtime_sec = runtime_sec,
        stringsAsFactors = FALSE
      )
    }

    # Advanced methods
    if (exists("select_K_ensemble", mode = "function")) {
      start_time <- Sys.time()
      K_hat <- tryCatch({
        suppressMessages({
          res <- select_K_ensemble(list_of_data, Kvec, r_fixed,
                                   c("BIC", "ICL_T1", "ICL_T2"), max_iter, FALSE)
          res$best_K
        })
      }, error = function(e) NA)
      results_list[[length(results_list) + 1]] <- data.frame(
        run_id = p$run_id, seed = p$seed, N = p$N,
        K_true = p$K_true, r_true = p$r_true,
        cluster_sep_spatial = p$cluster_sep_spatial,
        method = "Ensemble", K_hat = K_hat,
        K_correct = (!is.na(K_hat) && K_hat == p$K_true),
        runtime_sec = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
        stringsAsFactors = FALSE
      )
    }

    if (exists("select_K_adaptive_temperature", mode = "function")) {
      start_time <- Sys.time()
      K_hat <- tryCatch({
        suppressMessages({
          res <- select_K_adaptive_temperature(list_of_data, Kvec, r_fixed, max_iter, FALSE)
          res$best_K
        })
      }, error = function(e) NA)
      results_list[[length(results_list) + 1]] <- data.frame(
        run_id = p$run_id, seed = p$seed, N = p$N,
        K_true = p$K_true, r_true = p$r_true,
        cluster_sep_spatial = p$cluster_sep_spatial,
        method = "Adaptive", K_hat = K_hat,
        K_correct = (!is.na(K_hat) && K_hat == p$K_true),
        runtime_sec = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
        stringsAsFactors = FALSE
      )
    }

    if (exists("select_K_silhouette", mode = "function")) {
      start_time <- Sys.time()
      K_hat <- tryCatch({
        suppressMessages({
          res <- select_K_silhouette(list_of_data, Kvec, r_fixed, max_iter, FALSE)
          res$best_K
        })
      }, error = function(e) NA)
      results_list[[length(results_list) + 1]] <- data.frame(
        run_id = p$run_id, seed = p$seed, N = p$N,
        K_true = p$K_true, r_true = p$r_true,
        cluster_sep_spatial = p$cluster_sep_spatial,
        method = "Silhouette", K_hat = K_hat,
        K_correct = (!is.na(K_hat) && K_hat == p$K_true),
        runtime_sec = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
        stringsAsFactors = FALSE
      )
    }

    dplyr::bind_rows(results_list)
  }

  # Run sweep
  if (verbose) cat("Starting parameter sweep...\n")
  start_time <- Sys.time()

  if (cores > 1 && .Platform$OS.type != "windows") {
    if (verbose) cat("Using", cores, "cores\n\n")
    results_list <- parallel::mclapply(
      seq_len(nrow(param_grid)),
      function(i) run_single(param_grid[i, ]),
      mc.cores = cores
    )
  } else {
    if (verbose) cat("Sequential processing\n\n")
    results_list <- lapply(seq_len(nrow(param_grid)), function(i) {
      if (verbose && i %% 10 == 0) {
        cat(sprintf("  %d/%d\n", i, nrow(param_grid)))
      }
      run_single(param_grid[i, ])
    })
  }

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  if (verbose) cat(sprintf("\nCompleted in %.1f minutes\n\n", elapsed))

  all_results <- dplyr::bind_rows(results_list)

  # Save results
  timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
  raw_file <- file.path(output_dir, paste0("k_selection_raw_", timestamp_str, ".csv"))
  write.csv(all_results, raw_file, row.names = FALSE)

  # Summary
  if (verbose) {
    cat("=== Summary ===\n\n")
    overall_acc <- all_results |>
      dplyr::group_by(method) |>
      dplyr::summarise(
        n = dplyr::n(),
        accuracy = mean(K_correct, na.rm = TRUE) * 100,
        .groups = "drop"
      ) |>
      dplyr::arrange(dplyr::desc(accuracy))
    print(as.data.frame(overall_acc), row.names = FALSE)
    cat("\nResults saved to:", raw_file, "\n")
  }

  invisible(all_results)
}
