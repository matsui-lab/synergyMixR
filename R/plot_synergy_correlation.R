#' Plot Synergy Correlation Heatmap Between Clusters (Internal)
#'
#' Creates a correlation heatmap comparing synergy patterns (factor loadings)
#' between two clusters from a fitted model. **Internal function**.
#'
#' @param model A fitted MFA or MPCA model containing Lambda or W matrices
#' @param clusterA Index of the first cluster to compare
#' @param clusterB Index of the second cluster to compare  
#' @param method Correlation method: "pearson", "kendall", or "spearman"
#'
#' @return A ggplot object showing the correlation heatmap
#'
#' @examples
#' \dontrun{
#' # For MFA model
#' p <- plot_synergy_correlation(mfa_model, clusterA = 1, clusterB = 2)
#' print(p)
#' 
#' # For MPCA model  
#' p <- plot_synergy_correlation(mpca_model, clusterA = 1, clusterB = 3, method = "spearman")
#' print(p)
#' }
#'
#' @noRd
plot_synergy_correlation <- function(model, clusterA, clusterB, method = "pearson") {
  method <- match.arg(method, c("pearson", "kendall", "spearman"))
  
  if (!is.null(model$Lambda)) {
    loadings_A <- model$Lambda[[clusterA]]
    loadings_B <- model$Lambda[[clusterB]]
    model_type <- "MFA"
  } else if (!is.null(model$W)) {
    loadings_A <- model$W[[clusterA]]
    loadings_B <- model$W[[clusterB]]
    model_type <- "MPCA"
  } else {
    stop("Model must contain either Lambda (MFA) or W (MPCA) matrices")
  }
  
  if (is.null(loadings_A) || is.null(loadings_B)) {
    stop("Specified clusters not found in model")
  }
  
  if (ncol(loadings_A) != ncol(loadings_B)) {
    stop("Clusters must have the same number of factors for correlation analysis")
  }
  
  r_factors <- ncol(loadings_A)
  cor_matrix <- matrix(0, nrow = r_factors, ncol = r_factors)
  
  for (i in seq_len(r_factors)) {
    for (j in seq_len(r_factors)) {
      cor_matrix[i, j] <- cor(loadings_A[, i], loadings_B[, j], method = method)
    }
  }
  
  rownames(cor_matrix) <- paste0("Cluster", clusterA, "_F", seq_len(r_factors))
  colnames(cor_matrix) <- paste0("Cluster", clusterB, "_F", seq_len(r_factors))
  
  cor_df <- expand.grid(
    Factor_A = rownames(cor_matrix),
    Factor_B = colnames(cor_matrix),
    stringsAsFactors = FALSE
  )
  cor_df$Correlation <- as.vector(cor_matrix)
  cor_df$Factor_A <- factor(cor_df$Factor_A, levels = rownames(cor_matrix))
  cor_df$Factor_B <- factor(cor_df$Factor_B, levels = colnames(cor_matrix))
  
  p <- ggplot2::ggplot(cor_df, ggplot2::aes(x = Factor_B, y = Factor_A, fill = Correlation)) +
    ggplot2::geom_tile() +
    scale_fill_mixsynergy_diverging(
      midpoint = 0,
      limits = c(-1, 1),
      name = paste0(stringr::str_to_title(method), "\nCorrelation")
    ) +
    ggplot2::geom_text(ggplot2::aes(label = round(Correlation, 2)), size = 3) +
    ggplot2::labs(
      title = paste0("Synergy Correlation: Cluster ", clusterA, " vs Cluster ", clusterB),
      subtitle = paste0(model_type, " Model - ", stringr::str_to_title(method), " Correlation"),
      x = paste0("Cluster ", clusterB, " Factors"),
      y = paste0("Cluster ", clusterA, " Factors")
    ) +
    theme_mixsynergy_heatmap() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}
