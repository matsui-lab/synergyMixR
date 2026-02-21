#' Unified Theme and Style Definitions for synergyMixR Visualizations
#'
#' This file provides consistent theme and style definitions for all
#' visualization functions in the synergyMixR package, ensuring
#' publication-quality graphics with unified aesthetics.

#' Get Unified Theme for synergyMixR Plots
#'
#' Returns a consistent ggplot2 theme for all package visualizations.
#' Based on theme_minimal() with customizations for publication quality.
#'
#' @param base_size Base font size (default: 11)
#' @param base_family Base font family (default: "")
#'
#' @return A ggplot2 theme object
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(data, aes(x, y)) + geom_point() + theme_mixsynergy()
#' }
theme_mixsynergy <- function(base_size = 11, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = base_size * 1.2),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = base_size * 0.9),
      axis.title = ggplot2::element_text(face = "bold", size = base_size),
      axis.text = ggplot2::element_text(size = base_size * 0.9),
      legend.title = ggplot2::element_text(face = "bold", size = base_size),
      legend.text = ggplot2::element_text(size = base_size * 0.9),
      strip.text = ggplot2::element_text(face = "bold", size = base_size),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, color = "grey70", linewidth = 0.5)
    )
}

#' Get Unified Heatmap Theme for synergyMixR
#'
#' Returns a specialized theme for heatmap visualizations with minimal grid lines.
#'
#' @param base_size Base font size (default: 11)
#' @param base_family Base font family (default: "")
#'
#' @return A ggplot2 theme object optimized for heatmaps
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(data, aes(x, y, fill = value)) + 
#'   geom_tile() + 
#'   theme_mixsynergy_heatmap()
#' }
theme_mixsynergy_heatmap <- function(base_size = 11, base_family = "") {
  theme_mixsynergy(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(color = "grey70")
    )
}

#' Get Unified Viridis Color Scale for Continuous Values
#'
#' Returns a consistent viridis color scale for continuous fill aesthetics.
#' Uses the "viridis" option by default for colorblind-friendly visualization.
#'
#' @param option Viridis color option: "viridis" (default), "magma", "plasma", 
#'   "inferno", "cividis", "rocket", "mako", "turbo"
#' @param direction Direction of color scale: 1 (default) or -1 (reversed)
#' @param ... Additional arguments passed to scale_fill_viridis_c()
#'
#' @return A ggplot2 scale object
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(data, aes(x, y, fill = value)) + 
#'   geom_tile() + 
#'   scale_fill_mixsynergy()
#' }
scale_fill_mixsynergy <- function(option = "viridis", direction = 1, ...) {
  ggplot2::scale_fill_viridis_c(option = option, direction = direction, ...)
}

#' Get Unified Viridis Color Scale for Diverging Values
#'
#' Returns a consistent color scale for diverging data (e.g., correlations, loadings).
#' Uses a diverging palette centered at zero.
#'
#' @param option Viridis color option for diverging scale (default: "RdBu")
#' @param midpoint Midpoint value for diverging scale (default: 0)
#' @param limits Limits for the color scale (default: NULL for automatic)
#' @param ... Additional arguments passed to scale_fill_gradient2()
#'
#' @return A ggplot2 scale object
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(data, aes(x, y, fill = correlation)) + 
#'   geom_tile() + 
#'   scale_fill_mixsynergy_diverging()
#' }
scale_fill_mixsynergy_diverging <- function(option = "RdBu", 
                                            midpoint = 0, 
                                            limits = NULL, 
                                            ...) {
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",      # Blue
    mid = "#F7F7F7",      # Light gray
    high = "#B2182B",     # Red
    midpoint = midpoint,
    limits = limits,
    ...
  )
}

#' Get Unified Fill Scale for Discrete Categories
#'
#' Returns a consistent fill scale for discrete categories using viridis.
#' Use this for bar plots, boxplots, and other geoms with discrete fill aesthetics.
#'
#' @param option Viridis color option: "viridis" (default), "magma", "plasma", 
#'   "inferno", "cividis", "rocket", "mako", "turbo"
#' @param ... Additional arguments passed to scale_fill_viridis_d()
#'
#' @return A ggplot2 scale object
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(data, aes(x, y, fill = cluster)) + 
#'   geom_col() + 
#'   scale_fill_mixsynergy_discrete()
#' }
scale_fill_mixsynergy_discrete <- function(option = "viridis", ...) {
  ggplot2::scale_fill_viridis_d(option = option, ...)
}

#' Get Unified Color Scale for Discrete Categories
#'
#' Returns a consistent color scale for discrete categories using viridis.
#'
#' @param option Viridis color option: "viridis" (default), "magma", "plasma", 
#'   "inferno", "cividis", "rocket", "mako", "turbo"
#' @param ... Additional arguments passed to scale_color_viridis_d()
#'
#' @return A ggplot2 scale object
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(data, aes(x, y, color = cluster)) + 
#'   geom_point() + 
#'   scale_color_mixsynergy()
#' }
scale_color_mixsynergy <- function(option = "viridis", ...) {
  ggplot2::scale_color_viridis_d(option = option, ...)
}

#' Apply Consistent Faceting for Cluster and Factor Visualizations
#'
#' Helper function to apply consistent facet_grid() styling for cluster-factor plots.
#'
#' @param rows Variable(s) for facet rows (e.g., "Cluster")
#' @param cols Variable(s) for facet columns (e.g., "Factor")
#' @param scales Scale behavior: "fixed" (default), "free", "free_x", "free_y"
#' @param ... Additional arguments passed to facet_grid()
#'
#' @return A ggplot2 facet_grid object
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(data, aes(x, y)) + 
#'   geom_line() + 
#'   facet_mixsynergy(rows = "Cluster", cols = "Factor", scales = "free_y")
#' }
facet_mixsynergy <- function(rows = NULL, cols = NULL, scales = "fixed", ...) {
  if (is.null(rows) && is.null(cols)) {
    stop("At least one of 'rows' or 'cols' must be specified")
  }
  
  facet_formula <- NULL
  if (!is.null(rows) && !is.null(cols)) {
    facet_formula <- as.formula(paste(rows, "~", cols))
  } else if (!is.null(rows)) {
    facet_formula <- as.formula(paste(rows, "~ ."))
  } else {
    facet_formula <- as.formula(paste(". ~", cols))
  }
  
  ggplot2::facet_grid(facet_formula, scales = scales, ...)
}
