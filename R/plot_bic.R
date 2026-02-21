#' Enhanced BIC Plotting with Minimum Point Highlighting (Internal)
#'
#' Creates a BIC plot with enhanced visualization showing the minimum BIC point
#' with special highlighting and text labels. **Internal function**.
#'
#' @param df_summary A data frame containing BIC results with columns K, r, BIC
#' @param x_axis Which variable to use on x-axis: "K" or "r"
#'
#' @return A ggplot object showing BIC values with minimum point highlighted
#'
#' @examples
#' \dontrun{
#' # Assuming df_summary has columns K, r, BIC
#' p <- plot_bic_enhanced(df_summary, x_axis = "K")
#' print(p)
#' }
#'
#' @noRd
plot_bic_enhanced <- function(df_summary, x_axis = c("K","r")) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required for enhanced BIC plotting. Please install it.")
  }
  
  x_axis <- match.arg(x_axis)
  
  min_bic_idx <- which.min(df_summary$BIC)
  min_bic_row <- df_summary[min_bic_idx, ]
  
  if (x_axis == "K") {
    p <- ggplot2::ggplot(df_summary, ggplot2::aes(x = K, y = BIC, color = factor(r))) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::geom_point(data = min_bic_row, size = 4, color = "red") +
      ggrepel::geom_text_repel(
        data = min_bic_row,
        ggplot2::aes(label = paste0("Min BIC: ", round(BIC, 2))),
        color = "red",
        size = 3.5
      ) +
      ggplot2::labs(color = "r") +
      ggplot2::scale_x_continuous(breaks = sort(unique(df_summary$K)))
  } else {
    p <- ggplot2::ggplot(df_summary, ggplot2::aes(x = r, y = BIC, color = factor(K))) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::geom_point(data = min_bic_row, size = 4, color = "red") +
      ggrepel::geom_text_repel(
        data = min_bic_row,
        ggplot2::aes(label = paste0("Min BIC: ", round(BIC, 2))),
        color = "red",
        size = 3.5
      ) +
      ggplot2::labs(color = "K") +
      ggplot2::scale_x_continuous(breaks = sort(unique(df_summary$r)))
  }

  p <- p +
    ggplot2::labs(title = "BIC vs. (K, r) with Minimum Highlighted", x = x_axis, y = "BIC") +
    theme_mixsynergy()
  
  return(p)
}
