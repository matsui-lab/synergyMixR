#' Global variables to avoid R CMD check NOTEs
#' 
#' This file declares global variables used in ggplot2 and dplyr expressions
#' to avoid "no visible binding for global variable" warnings in R CMD check.
#' 
#' @noRd
utils::globalVariables(c(
  # Variables used in ggplot2 aesthetics and data manipulation
  "Time", "Score", "Subject", "Cluster", "Factor", "Channel", "Amplitude",
  "Series", "Loading", "Muscle", "Variable", "K", "r", "BIC", "Method", "ARI",
  "sep", "noise", "K_true", "r_true", "N", "M", "BIC_MFA", "ARI_PCA", "BIC_PCA",
  "ARI_MFA", "Factor_A", "Factor_B", "Correlation", "label_both",
  # Additional variables for plot functions
  "Dim1", "Dim2", "cluster_true", "cluster_est", "cluster_id", "value",
  "time_norm", "synergy", "muscle", "subject_id", "sigma_white", "sigma_low",
  "filter_length", "activation", "Estimated", "True", "Value",
  "Est_Component", "True_Component", "MethodType", ".data",
  "SSE", "SSE_plot", "subject", "vaf", "quantile", "component",
  "eigenvalue_obs", "eigenvalue_perm_95", "model_type", "win_rate",
  "selection_strategy", "metric", "n_wins", "method",
  "K_correct", "r_correct", "Kr_correct", "error_msg", "runtime_sec",
  "subspace_similarity", "n", "total", "y", "ymin", "ymax", "Metric",
  "Mean", "SD", "Lower", "Upper", "Type", "Component", "aic", "bic",
  "VAF_mean", "ARI_mean", "BIC_mean", "SSE_mean",
  # Additional NSE variables for ggplot2 and dplyr contexts
  "VAF", "accuracy", "mean_value", "time", "cv_mean", "cv_se",
  "sd_value", "min_value", "max_value",
  # Functions checked with exists() before calling - not actual globals but 
  # referenced conditionally for optional functionality
  "apply_temperature", "select_K_ensemble", "select_K_adaptive_temperature",
  "select_K_silhouette", "make_param_grid_paper_core", "make_param_grid_paper_stress"
))

#' Strip All Attributes Except Essential Ones
#'
#' Removes all attributes from an object except for names, row.names, and class.
#' This is useful for testing scenarios where you want to simulate data frames
#' without custom attributes (e.g., fitted model objects stored as attributes).
#'
#' @param x An R object (typically a data frame)
#'
#' @return The same object with all non-essential attributes removed
#'
#' @examples
#' \dontrun{
#' # Create a data frame with custom attributes
#' df <- data.frame(x = 1:3, y = 4:6)
#' attr(df, "custom_attr") <- "some_value"
#' 
#' # Strip custom attributes
#' df_clean <- strip_attributes(df)
#' # df_clean will have only names, row.names, and class attributes
#' }
#'
#' @export
strip_attributes <- function(x) {
  # Get essential attributes to preserve
  essential_attrs <- list(
    names = attr(x, "names"),
    row.names = attr(x, "row.names"),
    class = attr(x, "class")
  )
  
  # Remove all attributes
  attributes(x) <- NULL
  
  # Restore essential attributes
  if (!is.null(essential_attrs$names)) {
    attr(x, "names") <- essential_attrs$names
  }
  if (!is.null(essential_attrs$row.names)) {
    attr(x, "row.names") <- essential_attrs$row.names
  }
  if (!is.null(essential_attrs$class)) {
    attr(x, "class") <- essential_attrs$class
  }
  
  return(x)
}
