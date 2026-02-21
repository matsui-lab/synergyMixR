# Tests for inst/paper directory structure validation
# Ensures both legacy compatibility wrappers and v2 scripts exist and are parseable
#
# These tests work in two contexts:
# 1. Installed package: system.file("paper/...", package="synergyMixR") returns valid paths
# 2. Development mode: Falls back to file.path(repo_root, "inst", ...) using DESCRIPTION detection
#
# Tests skip gracefully if neither context provides valid paths.

local({

  # ---------------------------------------------------------------------------
  # Helper: Find repository root by searching for DESCRIPTION file
  # ---------------------------------------------------------------------------

  find_repo_root <- function(start = getwd(), max_depth = 10) {
    cur <- normalizePath(start, winslash = "/", mustWork = FALSE)
    for (i in seq_len(max_depth)) {
      if (file.exists(file.path(cur, "DESCRIPTION"))) {
        return(cur)
      }
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
    NULL
  }

  # Cache repo root once for all tests
  repo_root <- find_repo_root()

  # ---------------------------------------------------------------------------
  # Helper: Resolve paths in both installed and development contexts
  # ---------------------------------------------------------------------------
  resolve_paper_path <- function(rel_path) {
    # Try installed package path first
    pkg_path <- system.file(rel_path, package = "synergyMixR")
    if (nzchar(pkg_path) && file.exists(pkg_path)) {
      return(pkg_path)
    }
    # Fall back to development path (repo root based, not WD dependent)
    if (!is.null(repo_root)) {
      dev_path <- file.path(repo_root, "inst", rel_path)
      if (file.exists(dev_path)) {
        return(dev_path)
      }
    }
    NULL
  }

  # ---------------------------------------------------------------------------
  # Helper: Assert files exist (with skip if none found)
  # ---------------------------------------------------------------------------
  assert_files_exist <- function(paths, label) {
    resolved <- lapply(paths, resolve_paper_path)
    existing <- Filter(Negate(is.null), resolved)

    skip_if(length(existing) == 0,
            paste("No", label, "found in installed or development context"))

    for (i in seq_along(paths)) {
      path <- resolved[[i]]
      if (!is.null(path)) {
        expect_true(file.exists(path),
                    info = paste(label, "missing:", paths[i]))
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Helper: Assert R files are parseable (with skip if none found)
  # Automatically filters to .R files only
  # ---------------------------------------------------------------------------
  assert_files_parseable <- function(paths, label) {
    # Filter to .R files only
    r_paths <- paths[grepl("\\.R$", paths, ignore.case = TRUE)]
    if (length(r_paths) == 0) {
      skip(paste("No .R files in", label, "to parse"))
    }

    resolved <- lapply(r_paths, resolve_paper_path)
    existing <- Filter(Negate(is.null), resolved)

    skip_if(length(existing) == 0,
            paste("No", label, "found in installed or development context"))

    for (i in seq_along(r_paths)) {
      path <- resolved[[i]]
      if (!is.null(path)) {
        # Use tryCatch to provide informative error messages
        result <- tryCatch(
          { parse(file = path); TRUE },
          error = function(e) e$message
        )
        expect_true(
          isTRUE(result),
          info = paste("Failed to parse:", r_paths[i],
                       if (!isTRUE(result)) paste("-", result) else "")
        )
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Path definitions (centralized for easy maintenance)
  # ---------------------------------------------------------------------------
  legacy_wrappers <- c(
    "paper/config.R",
    "paper/01_simulation_design.R",
    "paper/02_run_grid_search.R",
    "paper/03_aggregate_results.R",
    "paper/03b_compare_baselines.R",
    "paper/04_make_figures.R",
    "paper/bic_test.R"
  )

  legacy_scripts <- c(
    "paper/legacy/config.R",
    "paper/legacy/01_simulation_design.R",
    "paper/legacy/02_run_grid_search.R",
    "paper/legacy/03_aggregate_results.R",
    "paper/legacy/03b_compare_baselines.R",
    "paper/legacy/04_make_figures.R",
    "paper/legacy/bic_test.R",
    "paper/legacy/README_legacy.md"
  )

  v2_scripts <- c(
    "paper/v2/config_param_grid.R",
    "paper/v2/helper_functions.R",
    "paper/v2/run_param_sweep.R",
    "paper/v2/plot_param_sweep_results.R",
    "paper/v2/run_param_sweep_modelsel.R",
    "paper/v2/plot_param_sweep_modelsel.R",
    "paper/v2/methods_comparison_figure.R",
    "paper/v2/method_comparison_figure.R",
    "paper/v2/JBHI_figure_plan.md"
  )

  v2_wrappers <- c(
    "paper/config_param_grid.R",
    "paper/helper_functions.R",
    "paper/run_param_sweep.R",
    "paper/plot_param_sweep_results.R",
    "paper/run_param_sweep_modelsel.R",
    "paper/plot_param_sweep_modelsel.R",
    "paper/methods_comparison_figure.R",
    "paper/method_comparison_figure.R"
  )

  # ---------------------------------------------------------------------------
  # Tests: Legacy compatibility wrappers
  # ---------------------------------------------------------------------------
  test_that("legacy compatibility wrappers exist at original paths", {
    assert_files_exist(legacy_wrappers, "legacy wrappers")
  })

  test_that("legacy compatibility wrappers are parseable", {
    assert_files_parseable(legacy_wrappers, "legacy wrappers")
  })

  # ---------------------------------------------------------------------------
  # Tests: Legacy scripts in legacy/ directory
  # ---------------------------------------------------------------------------
  test_that("legacy scripts exist in legacy/ directory", {
    assert_files_exist(legacy_scripts, "legacy scripts")
  })

  test_that("legacy scripts in legacy/ directory are parseable", {
    assert_files_parseable(legacy_scripts, "legacy scripts")
  })

  # ---------------------------------------------------------------------------
  # Tests: v2 scripts in v2/ directory
  # ---------------------------------------------------------------------------
  test_that("v2 scripts exist in v2/ directory", {
    assert_files_exist(v2_scripts, "v2 scripts")
  })

  test_that("v2 R scripts are parseable", {
    assert_files_parseable(v2_scripts, "v2 scripts")
  })

  # ---------------------------------------------------------------------------
  # Tests: v2 compatibility wrappers
  # ---------------------------------------------------------------------------
  test_that("v2 compatibility wrappers exist at original paths", {
    assert_files_exist(v2_wrappers, "v2 wrappers")
  })

  test_that("v2 compatibility wrappers are parseable", {
    assert_files_parseable(v2_wrappers, "v2 wrappers")
  })

  # ---------------------------------------------------------------------------
  # Tests: README files
  # ---------------------------------------------------------------------------
  test_that("README.md exists and mentions v2 workflow", {
    readme_path <- resolve_paper_path("paper/README.md")

    skip_if(is.null(readme_path),
            "README.md not found in installed or development context")

    expect_true(file.exists(readme_path), info = "README.md missing")

    readme_content <- readLines(readme_path)
    readme_text <- paste(readme_content, collapse = "\n")

    expect_true(
      grepl("v2", readme_text, ignore.case = TRUE),
      info = "README should mention v2 workflow"
    )

    expect_true(
      grepl("legacy", readme_text, ignore.case = TRUE),
      info = "README should mention legacy workflow"
    )
  })

  test_that("legacy/README_legacy.md exists", {
    assert_files_exist("paper/legacy/README_legacy.md", "legacy README")
  })

})
