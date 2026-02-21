# Plotting Functions for Parameter Sweep Results
#
# This file contains exported plotting functions for visualizing parameter sweep
# results interactively in R sessions, notebooks, or Shiny apps.

# Internal helper: Validate required columns
require_columns <- function(df, cols, fn_name) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(
      "Function '", fn_name, "' requires the following columns: ",
      paste(cols, collapse = ", "), "\n",
      "Missing columns: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Theme for Parameter Sweep Plots
#'
#' A consistent theme for parameter sweep visualization that ensures readable
#' text, clear facets, and prevents layout collapse.
#'
#' @param base_size Base font size. Default is 12.
#' @param base_family Base font family. Default is "".
#'
#' @return A ggplot2 theme object
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point() +
#'   theme_sweep()
#' }
#'
#' @export
theme_sweep <- function(base_size = 12, base_family = "") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.")
  }
  
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(size = base_size),
      axis.title = ggplot2::element_text(size = base_size + 2, face = "bold"),
      axis.text = ggplot2::element_text(size = base_size - 1),
      legend.title = ggplot2::element_text(size = base_size, face = "bold"),
      legend.text = ggplot2::element_text(size = base_size - 1),
      strip.text = ggplot2::element_text(size = base_size, face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "gray80", fill = NA, linewidth = 0.5)
    )
}

#' Plot Metric vs Axis for Parameter Sweep
#'
#' Creates a line plot showing how a performance metric varies with a parameter
#' axis (e.g., cluster separation). Useful for exploring parameter sensitivity.
#'
#' @param aggregated_df Data frame with aggregated results from \code{run_param_sweep()}.
#'   Must contain columns: Method, model_type, the specified metric columns
#'   (metric_mean, metric_sd), and the specified axis column.
#' @param metric Character. Metric to plot. One of "ARI", "VAF", "BIC", "SSE".
#'   Default is "ARI".
#' @param axis Character. Parameter axis to plot on x-axis. One of
#'   "cluster_sep_spatial", "cluster_sep_temporal", "cluster_sep_stability",
#'   "N", "K", "r". Default is "cluster_sep_spatial".
#' @param facet_by Character. Variable to facet by. One of "model_type" or "Method".
#'   Default is "model_type".
#' @param color_by Character. Variable to color by. One of "Method" or "model_type".
#'   Default is "Method".
#' @param show_error_bars Logical. Whether to show error bars (mean +/- sd).
#'   Default is TRUE.
#'
#' @return A ggplot2 object
#'
#' @examples
#' \dontrun{
#' library(synergyMixR)
#' results <- run_param_sweep(test = TRUE, cores = 2)
#' 
#' # Plot ARI vs spatial separation
#' p1 <- plot_metric_vs_axis(
#'   results$aggregated_results,
#'   metric = "ARI",
#'   axis = "cluster_sep_spatial"
#' )
#' print(p1)
#' 
#' # Plot VAF vs temporal separation
#' p2 <- plot_metric_vs_axis(
#'   results$aggregated_results,
#'   metric = "VAF",
#'   axis = "cluster_sep_temporal"
#' )
#' print(p2)
#' }
#'
#' @export
plot_metric_vs_axis <- function(aggregated_df,
                                metric = "ARI",
                                axis = "cluster_sep_spatial",
                                facet_by = "model_type",
                                color_by = "Method",
                                show_error_bars = TRUE) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Please install it.")
  }
  
  # Validate inputs
  metric <- match.arg(metric, c("ARI", "VAF", "BIC", "SSE"))
  axis <- match.arg(axis, c("cluster_sep_spatial", "cluster_sep_temporal",
                            "cluster_sep_stability", "N", "K", "r"))
  facet_by <- match.arg(facet_by, c("model_type", "Method"))
  color_by <- match.arg(color_by, c("Method", "model_type"))
  
  # Check required columns
  metric_mean_col <- paste0(metric, "_mean")
  metric_sd_col <- paste0(metric, "_sd")
  required_cols <- c("Method", "model_type", axis, metric_mean_col)
  if (show_error_bars) {
    required_cols <- c(required_cols, metric_sd_col)
  }
  require_columns(aggregated_df, required_cols, "plot_metric_vs_axis")
  
  # Create base plot
  p <- ggplot2::ggplot(
    aggregated_df,
    ggplot2::aes(
      x = .data[[axis]],
      y = .data[[metric_mean_col]],
      color = .data[[color_by]],
      group = .data[[color_by]]
    )
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2)
  
  # Add error bars if requested
  if (show_error_bars) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = .data[[metric_mean_col]] - .data[[metric_sd_col]],
        ymax = .data[[metric_mean_col]] + .data[[metric_sd_col]]
      ),
      width = 0.05,
      alpha = 0.5
    )
  }
  
  # Add faceting
  p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by)), ncol = 2)
  
  # Add color scale
  if (color_by == "Method") {
    p <- p + scale_color_mixsynergy()
  }
  
  # Add labels
  axis_label <- switch(
    axis,
    cluster_sep_spatial = "Spatial Separation",
    cluster_sep_temporal = "Temporal Separation",
    cluster_sep_stability = "Stability Separation",
    N = "Number of Subjects (N)",
    K = "Number of Clusters (K)",
    r = "Number of Synergies (r)"
  )
  
  metric_label <- switch(
    metric,
    ARI = "Adjusted Rand Index (ARI)",
    VAF = "Variance Accounted For (VAF)",
    BIC = "Bayesian Information Criterion (BIC)",
    SSE = "Sum of Squared Errors (SSE)"
  )
  
  p <- p + ggplot2::labs(
    title = paste(metric_label, "vs", axis_label),
    x = axis_label,
    y = metric_label,
    color = color_by
  ) +
    theme_sweep()
  
  return(p)
}

#' Plot Summary Comparison for Parameter Sweep
#'
#' Creates a bar plot comparing method performance across all parameter
#' combinations, showing mean performance with error bars.
#'
#' @param aggregated_df Data frame with aggregated results from \code{run_param_sweep()}.
#'   Must contain columns: Method, model_type, and the specified metric columns.
#' @param metric Character. Metric to plot. One of "ARI", "VAF", "BIC", "SSE".
#'   Default is "ARI".
#' @param facet_by Character. Variable to facet by. One of "model_type" or "Method".
#'   Default is "model_type".
#' @param summary Character. Type of summary to show. One of "mean_sd" (mean +/- sd)
#'   or "min_max" (min to max range). Default is "mean_sd".
#'
#' @return A ggplot2 object
#'
#' @examples
#' \dontrun{
#' library(synergyMixR)
#' results <- run_param_sweep(test = TRUE, cores = 2)
#' 
#' # Plot overall ARI comparison
#' p <- plot_sweep_summary(results$aggregated_results, metric = "ARI")
#' print(p)
#' }
#'
#' @export
plot_sweep_summary <- function(aggregated_df,
                               metric = "ARI",
                               facet_by = "model_type",
                               summary = "mean_sd") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Please install it.")
  }
  
  # Validate inputs
  metric <- match.arg(metric, c("ARI", "VAF", "BIC", "SSE"))
  facet_by <- match.arg(facet_by, c("model_type", "Method"))
  summary <- match.arg(summary, c("mean_sd", "min_max"))
  
  # Check required columns
  metric_mean_col <- paste0(metric, "_mean")
  metric_sd_col <- paste0(metric, "_sd")
  required_cols <- c("Method", "model_type", metric_mean_col, metric_sd_col)
  require_columns(aggregated_df, required_cols, "plot_sweep_summary")
  
  # Compute summary across all parameter combinations
  summary_df <- aggregated_df %>%
    dplyr::group_by(.data$Method, .data$model_type) %>%
    dplyr::summarize(
      mean_value = mean(.data[[metric_mean_col]], na.rm = TRUE),
      sd_value = stats::sd(.data[[metric_mean_col]], na.rm = TRUE),
      min_value = min(.data[[metric_mean_col]], na.rm = TRUE),
      max_value = max(.data[[metric_mean_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create plot based on summary type
  if (summary == "mean_sd") {
    p <- ggplot2::ggplot(
      summary_df,
      ggplot2::aes(x = .data$Method, y = .data$mean_value, fill = .data$model_type)
    ) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7, alpha = 0.8) +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = .data$mean_value - .data$sd_value,
          ymax = .data$mean_value + .data$sd_value
        ),
        position = ggplot2::position_dodge(width = 0.8),
        width = 0.25,
        alpha = 0.6
      )
  } else {  # min_max
    p <- ggplot2::ggplot(
      summary_df,
      ggplot2::aes(x = .data$Method, y = .data$mean_value, fill = .data$model_type)
    ) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7, alpha = 0.8) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$min_value, ymax = .data$max_value),
        position = ggplot2::position_dodge(width = 0.8),
        width = 0.25,
        alpha = 0.6
      )
  }
  
  # Add faceting if requested
  if (facet_by == "model_type") {
    p <- p + ggplot2::facet_wrap(~ model_type, ncol = 2, scales = "free_x")
  }
  
  # Add fill scale (discrete for model_type)
  p <- p + scale_fill_mixsynergy_discrete()
  
  # Add labels
  metric_label <- switch(
    metric,
    ARI = "Adjusted Rand Index (ARI)",
    VAF = "Variance Accounted For (VAF)",
    BIC = "Bayesian Information Criterion (BIC)",
    SSE = "Sum of Squared Errors (SSE)"
  )
  
  error_label <- if (summary == "mean_sd") "Mean +/- SD" else "Min to Max"
  
  p <- p + ggplot2::labs(
    title = paste("Performance Comparison -", metric),
    subtitle = paste("Error bars show", error_label),
    x = "Method",
    y = metric_label,
    fill = "Model Type"
  ) +
    theme_sweep() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  return(p)
}

#' Plot BIC Comparison for Parameter Sweep
#'
#' Creates a boxplot or violin plot comparing BIC values across methods and
#' model types. Lower BIC values indicate better model fit.
#'
#' @param aggregated_df Data frame with aggregated results from \code{run_param_sweep()}.
#'   Must contain columns: Method, model_type, BIC_mean.
#' @param plot_type Character. Type of plot. One of "boxplot" or "violin".
#'   Default is "boxplot".
#'
#' @return A ggplot2 object
#'
#' @examples
#' \dontrun{
#' library(synergyMixR)
#' results <- run_param_sweep(test = TRUE, cores = 2)
#' 
#' # Plot BIC comparison
#' p <- plot_bic_comparison(results$aggregated_results)
#' print(p)
#' }
#'
#' @export
plot_bic_comparison <- function(aggregated_df, plot_type = "boxplot") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.")
  }
  
  # Validate inputs
  plot_type <- match.arg(plot_type, c("boxplot", "violin"))
  
  # Check required columns
  require_columns(aggregated_df, c("Method", "model_type", "BIC_mean"), "plot_bic_comparison")
  
  # Create base plot
  p <- ggplot2::ggplot(
    aggregated_df,
    ggplot2::aes(x = .data$Method, y = .data$BIC_mean, fill = .data$model_type)
  )
  
  # Add geom based on plot type
  if (plot_type == "boxplot") {
    p <- p + ggplot2::geom_boxplot(alpha = 0.7)
  } else {
    p <- p + ggplot2::geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75))
  }
  
  # Add fill scale (discrete for model_type)
  p <- p + scale_fill_mixsynergy_discrete()
  
  # Add labels
  p <- p + ggplot2::labs(
    title = "BIC Comparison Across Methods",
    subtitle = "Lower BIC indicates better model fit",
    x = "Method",
    y = "Bayesian Information Criterion (BIC)",
    fill = "Model Type"
  ) +
    theme_sweep() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  return(p)
}

#' Plot Win Rate Heatmap for Parameter Sweep
#'
#' Creates a heatmap showing the percentage of parameter combinations where
#' each method performed best for each metric.
#'
#' @param win_rates_df Data frame with win rates from \code{run_param_sweep()}.
#'   Must contain columns: Method, model_type, metric, win_rate.
#'   This is typically the data from the win_rates.csv file.
#' @param facet_by Character. Variable to facet by. One of "model_type" or "metric".
#'   Default is "metric".
#'
#' @return A ggplot2 object
#'
#' @examples
#' \dontrun{
#' library(synergyMixR)
#' results <- run_param_sweep(test = TRUE, cores = 2)
#' 
#' # Read win rates
#' win_rates <- read.csv(results$win_rates_file)
#' 
#' # Plot win rate heatmap
#' p <- plot_winrate_heatmap(win_rates)
#' print(p)
#' }
#'
#' @export
plot_winrate_heatmap <- function(win_rates_df, facet_by = "metric") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.")
  }
  
  # Validate inputs
  facet_by <- match.arg(facet_by, c("model_type", "metric"))
  
  # Check required columns
  require_columns(win_rates_df, c("Method", "model_type", "metric", "win_rate"),
                 "plot_winrate_heatmap")
  
  # Create heatmap
  p <- ggplot2::ggplot(
    win_rates_df,
    ggplot2::aes(x = .data$Method, y = .data$metric, fill = .data$win_rate)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 1) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f%%", .data$win_rate)),
      color = "black",
      size = 3.5,
      fontface = "bold"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#D73027",
      mid = "#FEE08B",
      high = "#1A9850",
      midpoint = 50,
      limits = c(0, 100),
      name = "Win Rate (%)"
    )
  
  # Add faceting (use fixed scales to allow coord_fixed)
  p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by)), scales = "fixed")
  
  # Add labels and theme
  p <- p + ggplot2::labs(
    title = "Method Win Rates Across Metrics",
    subtitle = "Percentage of parameter combinations where method performed best",
    x = "Method",
    y = "Metric"
  ) +
    theme_sweep() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::coord_fixed(ratio = 0.8)
  
  return(p)
}
