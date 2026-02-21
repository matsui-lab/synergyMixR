#' Utility Functions for Visualization
#'
#' This file contains utility functions used across visualization functions
#' to support consistent formatting and labeling.

#' Apply Muscle Names to Y-axis Labels (Internal)
#'
#' Helper function to apply muscle names to ggplot y-axis labels if provided.
#' **Internal function** - not intended for direct user access.
#'
#' @param p A ggplot object
#' @param muscle_names Optional character vector of muscle names for y-axis labels
#' @param channel_col_name Name of the column containing channel information
#'
#' @return Modified ggplot object with updated y-axis labels if muscle_names provided
#' @noRd
apply_muscle_names_to_plot <- function(p, muscle_names = NULL, channel_col_name = "Channel") {
  if (!is.null(muscle_names)) {
    p <- p + ggplot2::scale_y_discrete(labels = muscle_names)
  }
  return(p)
}

#' Get Default Channel Names (Internal)
#'
#' Generate default channel names for visualization when muscle names are not provided.
#' **Internal function** - not intended for direct user access.
#'
#' @param M Number of channels
#'
#' @return Character vector of default channel names
#' @noRd
get_default_channel_names <- function(M) {
  paste0("Channel", seq_len(M))
}
