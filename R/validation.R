# Input Validation Functions for synergyMixR
#
# This file provides standardized input validation functions used across
# the package to ensure consistent error messages and behavior.

#' Validate list_of_data Input
#'
#' Validates that the input is a proper list of matrices with consistent
#' dimensions suitable for muscle synergy analysis.
#'
#' @param list_of_data A list of matrices to validate
#' @param arg_name Character string for the argument name in error messages.
#'   Default: "list_of_data"
#' @param min_subjects Minimum number of subjects required. Default: 1
#' @param require_numeric Logical. If TRUE, require all elements to be numeric.
#'   Default: TRUE
#'
#' @return Invisibly returns TRUE if validation passes
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' data <- list(matrix(rnorm(100), 10, 10), matrix(rnorm(100), 10, 10))
#' validate_list_of_data(data)  # Returns TRUE
#'
#' validate_list_of_data(list())  # Error: must not be empty
#' }
validate_list_of_data <- function(list_of_data,
                                   arg_name = "list_of_data",
                                   min_subjects = 1,
                                   require_numeric = TRUE) {
  # Check if list

if (!is.list(list_of_data)) {
    stop(sprintf("'%s' must be a list, not %s",
                 arg_name, class(list_of_data)[1]),
         call. = FALSE)
  }

  # Check if empty
  if (length(list_of_data) == 0) {
    stop(sprintf("'%s' must not be empty", arg_name),
         call. = FALSE)
  }

  # Check minimum subjects
  if (length(list_of_data) < min_subjects) {
    stop(sprintf("'%s' must have at least %d subject(s), got %d",
                 arg_name, min_subjects, length(list_of_data)),
         call. = FALSE)
  }

  # Check all elements are matrices
  is_matrix <- vapply(list_of_data, is.matrix, logical(1))
  if (!all(is_matrix)) {
    bad_idx <- which(!is_matrix)[1]
    stop(sprintf("All elements of '%s' must be matrices. Element %d is %s",
                 arg_name, bad_idx, class(list_of_data[[bad_idx]])[1]),
         call. = FALSE)
  }

  # Check all matrices have same number of columns (muscles)
  n_cols <- vapply(list_of_data, ncol, integer(1))
  if (length(unique(n_cols)) != 1) {
    stop(sprintf(
      "All matrices in '%s' must have the same number of columns (muscles). Found: %s",
      arg_name, paste(unique(n_cols), collapse = ", ")),
      call. = FALSE)
  }

  # Check numeric if required
  if (require_numeric) {
    is_numeric <- vapply(list_of_data, function(x) is.numeric(x), logical(1))
    if (!all(is_numeric)) {
      bad_idx <- which(!is_numeric)[1]
      stop(sprintf("All matrices in '%s' must be numeric. Element %d is not numeric",
                   arg_name, bad_idx),
           call. = FALSE)
    }
  }

  # Check for NA/NaN/Inf values (warning only)
  has_na <- vapply(list_of_data, function(x) any(is.na(x)), logical(1))
  has_inf <- vapply(list_of_data, function(x) any(is.infinite(x)), logical(1))

  if (any(has_na)) {
    warning(sprintf("'%s' contains NA values in element(s): %s",
                    arg_name, paste(which(has_na), collapse = ", ")),
            call. = FALSE)
  }

  if (any(has_inf)) {
    warning(sprintf("'%s' contains Inf values in element(s): %s",
                    arg_name, paste(which(has_inf), collapse = ", ")),
            call. = FALSE)
  }

  invisible(TRUE)
}


#' Validate Model Parameters
#'
#' Validates parameters K (number of clusters), r (number of synergies/factors),
#' and max_iter (maximum iterations).
#'
#' @param K Number of clusters. Must be a positive integer.
#' @param r Number of synergies/factors. Must be a positive integer.
#' @param max_iter Maximum number of iterations. Must be a positive number.
#' @param n_init Number of initializations. Must be a positive integer if provided.
#' @param tol Convergence tolerance. Must be a positive number if provided.
#'
#' @return Invisibly returns TRUE if validation passes
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' validate_model_params(K = 2, r = 3, max_iter = 100)  # OK
#' validate_model_params(K = 0, r = 3, max_iter = 100)  # Error
#' validate_model_params(K = 2.5, r = 3, max_iter = 100)  # Error
#' }
validate_model_params <- function(K = NULL,
                                   r = NULL,
                                   max_iter = NULL,
                                   n_init = NULL,
                                   tol = NULL) {
  # Validate K
  if (!is.null(K)) {
    if (!is.numeric(K) || length(K) != 1) {
      stop("'K' must be a single numeric value", call. = FALSE)
    }
    if (K < 1 || K != floor(K)) {
      stop("'K' must be a positive integer (>= 1)", call. = FALSE)
    }
  }

  # Validate r
  if (!is.null(r)) {
    if (!is.numeric(r) || length(r) != 1) {
      stop("'r' must be a single numeric value", call. = FALSE)
    }
    if (r < 1 || r != floor(r)) {
      stop("'r' must be a positive integer (>= 1)", call. = FALSE)
    }
  }

  # Validate max_iter
  if (!is.null(max_iter)) {
    if (!is.numeric(max_iter) || length(max_iter) != 1) {
      stop("'max_iter' must be a single numeric value", call. = FALSE)
    }
    if (max_iter <= 0) {
      stop("'max_iter' must be positive", call. = FALSE)
    }
  }

  # Validate n_init
  if (!is.null(n_init)) {
    if (!is.numeric(n_init) || length(n_init) != 1) {
      stop("'n_init' must be a single numeric value", call. = FALSE)
    }
    if (n_init < 1 || n_init != floor(n_init)) {
      stop("'n_init' must be a positive integer (>= 1)", call. = FALSE)
    }
  }

  # Validate tol
  if (!is.null(tol)) {
    if (!is.numeric(tol) || length(tol) != 1) {
      stop("'tol' must be a single numeric value", call. = FALSE)
    }
    if (tol <= 0) {
      stop("'tol' must be positive", call. = FALSE)
    }
  }

  invisible(TRUE)
}


#' Validate K Against Number of Subjects
#'
#' Ensures that the number of clusters K does not exceed the number of subjects.
#'
#' @param K Number of clusters
#' @param N Number of subjects
#'
#' @return Invisibly returns TRUE if validation passes
#' @keywords internal
validate_K_vs_N <- function(K, N) {
  if (K > N) {
    stop(sprintf(
      "Number of clusters K (%d) cannot exceed number of subjects N (%d)",
      K, N),
      call. = FALSE)
  }
  invisible(TRUE)
}


#' Validate r Against Number of Muscles
#'
#' Ensures that the number of synergies r does not exceed the number of muscles.
#'
#' @param r Number of synergies/factors
#' @param M Number of muscles/variables
#'
#' @return Invisibly returns TRUE if validation passes
#' @keywords internal
validate_r_vs_M <- function(r, M) {
  if (r > M) {
    stop(sprintf(
      "Number of synergies r (%d) cannot exceed number of muscles M (%d)",
      r, M),
      call. = FALSE)
  }
  invisible(TRUE)
}


#' Validate Model Type
#'
#' Validates that model_type is one of the supported types.
#'
#' @param model_type Character string specifying the model type
#' @param allowed Character vector of allowed model types.
#'   Default: c("MFA", "MPCA")
#'
#' @return The validated (and possibly standardized) model_type
#' @keywords internal
validate_model_type <- function(model_type,
                                 allowed = c("MFA", "MPCA")) {
  if (is.null(model_type) || length(model_type) != 1) {
    stop(sprintf("'model_type' must be one of: %s",
                 paste(allowed, collapse = ", ")),
         call. = FALSE)
  }

  # Case-insensitive matching
  model_type_upper <- toupper(model_type)

  if (!model_type_upper %in% allowed) {
    stop(sprintf("'model_type' must be one of: %s. Got: '%s'",
                 paste(allowed, collapse = ", "), model_type),
         call. = FALSE)
  }

  model_type_upper
}


#' Validate Fitted Model Object
#'
#' Validates that a fitted model object has the expected structure.
#'
#' @param fit Fitted model object
#' @param model_type Expected model type ("MFA" or "MPCA"). If NULL, auto-detect.
#' @param required_fields Character vector of fields that must be present
#'
#' @return The detected or validated model_type
#' @keywords internal
validate_fit_object <- function(fit,
                                 model_type = NULL,
                                 required_fields = NULL) {
  if (!is.list(fit)) {
    stop("'fit' must be a list object", call. = FALSE)
  }

  # Auto-detect model type if not specified
  detected_type <- NULL
  if ("Lambda" %in% names(fit)) {
    detected_type <- "MFA"
  } else if ("W" %in% names(fit)) {
    detected_type <- "MPCA"
  }

  if (is.null(model_type)) {
    if (is.null(detected_type)) {
      stop("Cannot determine model type from fit object. Missing 'Lambda' (MFA) or 'W' (MPCA)",
           call. = FALSE)
    }
    model_type <- detected_type
  } else {
    model_type <- validate_model_type(model_type)
    if (!is.null(detected_type) && detected_type != model_type) {
      warning(sprintf(
        "Specified model_type '%s' does not match detected type '%s'",
        model_type, detected_type),
        call. = FALSE)
    }
  }

  # Check required fields
  if (!is.null(required_fields)) {
    missing_fields <- setdiff(required_fields, names(fit))
    if (length(missing_fields) > 0) {
      stop(sprintf("'fit' is missing required fields: %s",
                   paste(missing_fields, collapse = ", ")),
           call. = FALSE)
    }
  }

  model_type
}


#' Handle Deprecated Arguments
#'
#' Checks for deprecated argument names in `...` and returns the value
#' with a deprecation warning if found.
#'
#' @param dots Named list from `list(...)`
#' @param old_name Character. The deprecated argument name
#' @param new_name Character. The new argument name
#' @param new_value Current value of the new argument
#' @param pkg_name Package name for the warning message
#'
#' @return The value to use (from new_name if set, otherwise from old_name)
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' my_func <- function(new_arg = 5, ...) {
#'   dots <- list(...)
#'   new_arg <- handle_deprecated_arg(dots, "old_arg", "new_arg", new_arg)
#'   # ... use new_arg
#' }
#' }
handle_deprecated_arg <- function(dots,
                                   old_name,
                                   new_name,
                                   new_value,
                                   pkg_name = "synergyMixR") {
  if (old_name %in% names(dots)) {
    warning(sprintf(
      "Argument '%s' is deprecated. Use '%s' instead.",
      old_name, new_name),
      call. = FALSE)

    # If new_value wasn't explicitly set, use the deprecated value
    if (missing(new_value) || is.null(new_value)) {
      return(dots[[old_name]])
    }
    # If both are set, warn and prefer new value
    warning(sprintf(
      "Both '%s' and '%s' provided. Using '%s'.",
      old_name, new_name, new_name),
      call. = FALSE)
  }
  new_value
}


#' Standardized Argument Names Mapping
#'
#' Returns a named list mapping old argument names to new names.
#' Used internally for consistent deprecation handling.
#'
#' @return Named list where names are old arguments and values are new arguments
#' @keywords internal
get_deprecated_arg_map <- function() {
  list(
    nIterFA = "iter_fa",
    nIterPCA = "iter_pca",
    max_iterations = "max_iter",
    subject_rdim_for_kmeans = "kmeans_rdim",
    nc_cores = "n_cores"
  )
}
