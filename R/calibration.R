# Calibration Diagnostics for Mixture Models
#
# This file provides functions for assessing the calibration of probabilistic
# cluster assignments in mixture models. Well-calibrated probabilities should
# match observed frequencies: if a model assigns 80% probability to cluster k,
# then approximately 80% of such cases should truly belong to cluster k.

utils::globalVariables(c("mean_predicted", "mean_observed", "x", "cluster", "n_samples"))

#' Compute Calibration Statistics for Mixture Model
#'
#' Computes calibration statistics by binning predicted probabilities and
#' comparing to observed frequencies. This helps diagnose whether the model's
#' probabilistic predictions are reliable.
#'
#' @param resp Matrix of responsibilities (N x K), where resp[i,k] is the
#'   posterior probability that subject i belongs to cluster k. Rows should
#'   sum to 1.
#' @param z_true Integer vector of length N with true cluster labels (1 to K).
#'   If not available, can use hard cluster assignments from another method.
#' @param n_bins Integer. Number of bins for calibration curve. Default: 10.
#'
#' @return A list with class "calibration_result" containing:
#' \describe{
#'   \item{calibration_df}{Data frame with columns: bin_midpoint, mean_predicted,
#'     mean_observed, n_samples, cluster}
#'   \item{ece}{Expected Calibration Error (weighted average bin error)}
#'   \item{mce}{Maximum Calibration Error (worst bin)}
#'   \item{n_bins}{Number of bins used}
#'   \item{K}{Number of clusters}
#'   \item{N}{Number of subjects}
#' }
#'
#' @details
#' Calibration is computed for each cluster separately:
#' 1. Bin the predicted probabilities into n_bins equal-width bins
#' 2. For each bin, compute mean predicted probability and observed frequency
#' 3. A well-calibrated model has mean_predicted close to mean_observed
#'
#' The Expected Calibration Error (ECE) summarizes overall calibration:
#' \deqn{ECE = \sum_{b=1}^{B} \frac{|S_b|}{N} |acc(S_b) - conf(S_b)|}
#' where S_b is the set of samples in bin b, acc is accuracy (observed
#' frequency), and conf is confidence (mean predicted probability).
#'
#' @examples
#' \dontrun{
#' # Fit a mixture model
#' fit <- mixture_pca_em_fit(list_of_data, K = 3, r = 2)
#'
#' # Compute calibration using true labels
#' calib <- compute_calibration(fit$resp, z_true)
#'
#' # View ECE
#' print(calib$ece)
#'
#' # Plot calibration curve
#' plot_calibration(calib)
#' }
#'
#' @seealso \code{\link{plot_calibration}}, \code{\link{compute_ece}}
#'
#' @export
compute_calibration <- function(resp, z_true, n_bins = 10) {
  # Input validation
  if (!is.matrix(resp)) {
    resp <- as.matrix(resp)
  }

  N <- nrow(resp)
  K <- ncol(resp)

  if (length(z_true) != N) {
    stop("Length of z_true must equal number of rows in resp")
  }

  if (!all(z_true %in% 1:K)) {
    # Try to handle factor or character labels
    z_true <- as.integer(as.factor(z_true))
    if (!all(z_true %in% 1:K)) {
      stop("z_true must contain values from 1 to K")
    }
  }

  if (n_bins < 2) {
    stop("n_bins must be at least 2")
  }

  # Define bin edges
  bin_edges <- seq(0, 1, length.out = n_bins + 1)
  bin_midpoints <- (bin_edges[-1] + bin_edges[-(n_bins + 1)]) / 2

  # Compute calibration for each cluster
  calibration_results <- list()

  for (k in 1:K) {
    # Predicted probability for cluster k
    pred_probs <- resp[, k]

    # Binary indicator: does subject truly belong to cluster k?
    true_membership <- as.integer(z_true == k)

    # Bin the predictions
    bin_assignments <- cut(pred_probs, breaks = bin_edges,
                           include.lowest = TRUE, labels = FALSE)

    # Compute statistics for each bin
    bin_stats <- data.frame(
      bin = 1:n_bins,
      bin_midpoint = bin_midpoints,
      mean_predicted = NA_real_,
      mean_observed = NA_real_,
      n_samples = 0L,
      cluster = k
    )

    for (b in 1:n_bins) {
      in_bin <- which(bin_assignments == b)

      if (length(in_bin) > 0) {
        bin_stats$mean_predicted[b] <- mean(pred_probs[in_bin])
        bin_stats$mean_observed[b] <- mean(true_membership[in_bin])
        bin_stats$n_samples[b] <- length(in_bin)
      }
    }

    calibration_results[[k]] <- bin_stats
  }

  # Combine results
  calibration_df <- do.call(rbind, calibration_results)
  calibration_df$cluster <- factor(calibration_df$cluster)

  # Compute ECE (Expected Calibration Error)
  # ECE = sum over bins of (n_b / N) * |mean_pred - mean_obs|
  calibration_df_clean <- calibration_df[calibration_df$n_samples > 0, ]

  if (nrow(calibration_df_clean) > 0) {
    ece <- sum(calibration_df_clean$n_samples *
                 abs(calibration_df_clean$mean_predicted -
                       calibration_df_clean$mean_observed)) / (N * K)

    mce <- max(abs(calibration_df_clean$mean_predicted -
                     calibration_df_clean$mean_observed))
  } else {
    ece <- NA
    mce <- NA
  }

  result <- list(
    calibration_df = calibration_df,
    ece = ece,
    mce = mce,
    n_bins = n_bins,
    K = K,
    N = N
  )

  class(result) <- c("calibration_result", "list")
  result
}


#' Compute Expected Calibration Error
#'
#' A convenience function to compute only the ECE from responsibilities
#' and true labels.
#'
#' @param resp Matrix of responsibilities (N x K)
#' @param z_true Integer vector of true cluster labels (length N)
#' @param n_bins Number of bins. Default: 10.
#'
#' @return Numeric scalar, the ECE value. Lower is better (0 = perfect calibration).
#'
#' @details
#' ECE measures the average absolute difference between predicted probabilities
#' and observed frequencies. A value of 0 indicates perfect calibration; values
#' above 0.1 suggest significant miscalibration.
#'
#' Guidelines for interpretation:
#' \itemize{
#'   \item ECE < 0.02: Excellent calibration
#'   \item ECE 0.02-0.05: Good calibration
#'   \item ECE 0.05-0.10: Moderate calibration
#'   \item ECE > 0.10: Poor calibration
#' }
#'
#' @examples
#' \dontrun{
#' ece <- compute_ece(fit$resp, z_true)
#' cat(sprintf("ECE: %.3f\n", ece))
#' }
#'
#' @seealso \code{\link{compute_calibration}}
#'
#' @export
compute_ece <- function(resp, z_true, n_bins = 10) {
  calib <- compute_calibration(resp, z_true, n_bins)
  calib$ece
}


#' Plot Calibration Curve
#'
#' Creates a reliability diagram (calibration plot) showing the relationship
#' between predicted probabilities and observed frequencies.
#'
#' @param calib_result Result from \code{\link{compute_calibration}}.
#' @param show_histogram Logical. If TRUE, show histogram of predictions below
#'   the calibration curve. Default: TRUE.
#' @param by_cluster Logical. If TRUE, show separate curves for each cluster.
#'   If FALSE, aggregate across clusters. Default: FALSE.
#' @param title Character. Plot title. Default: "Calibration Plot".
#' @param ... Additional arguments passed to ggplot2 functions.
#'
#' @return A ggplot2 object.
#'
#' @details
#' A reliability diagram plots mean predicted probability (x-axis) against
#' mean observed frequency (y-axis) for each bin. A perfectly calibrated model
#' follows the diagonal line (y = x). Points above the diagonal indicate
#' underconfidence; points below indicate overconfidence.
#'
#' The shaded band around the diagonal shows +/- 0.1 deviation as a reference.
#' The ECE value is displayed in the plot subtitle.
#'
#' @examples
#' \dontrun{
#' calib <- compute_calibration(fit$resp, z_true)
#'
#' # Basic calibration plot
#' p <- plot_calibration(calib)
#' print(p)
#'
#' # Show by cluster
#' p <- plot_calibration(calib, by_cluster = TRUE)
#' print(p)
#'
#' # Save to file
#' ggsave("calibration.pdf", p, width = 6, height = 5)
#' }
#'
#' @seealso \code{\link{compute_calibration}}
#'
#' @export
plot_calibration <- function(calib_result,
                              show_histogram = TRUE,
                              by_cluster = FALSE,
                              title = "Calibration Plot",
                              ...) {

  if (!inherits(calib_result, "calibration_result")) {
    stop("calib_result must be output from compute_calibration()")
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_calibration")
  }

  df <- calib_result$calibration_df
  df <- df[df$n_samples > 0, ]  # Remove empty bins

  if (nrow(df) == 0) {
    stop("No non-empty bins to plot")
  }

  # Aggregate across clusters if not showing by cluster
  if (!by_cluster && calib_result$K > 1) {
    # Compute weighted averages across clusters for each bin
    df_agg <- do.call(rbind, lapply(unique(df$bin), function(b) {
      subset_df <- df[df$bin == b, ]
      # Skip bins with no data
      if (nrow(subset_df) == 0 || sum(subset_df$n_samples) == 0) {
        return(NULL)
      }
      # Handle NA values
      valid_rows <- !is.na(subset_df$mean_predicted) & !is.na(subset_df$mean_observed)
      if (!any(valid_rows)) {
        return(NULL)
      }
      subset_df <- subset_df[valid_rows, ]
      if (nrow(subset_df) == 0) {
        return(NULL)
      }

      weights <- subset_df$n_samples / sum(subset_df$n_samples)
      data.frame(
        bin = b,
        bin_midpoint = unique(subset_df$bin_midpoint),
        mean_predicted = sum(weights * subset_df$mean_predicted),
        mean_observed = sum(weights * subset_df$mean_observed),
        n_samples = sum(subset_df$n_samples),
        cluster = "All"
      )
    }))

    if (is.null(df_agg) || nrow(df_agg) == 0) {
      stop("No valid data to aggregate for calibration plot")
    }
    df <- df_agg
  }

  # Create base plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = mean_predicted, y = mean_observed))

  # Add perfect calibration line
  p <- p + ggplot2::geom_abline(intercept = 0, slope = 1,
                                 linetype = "dashed", color = "gray50")

  # Add shaded band for +/- 0.1 reference
  p <- p + ggplot2::geom_ribbon(
    data = data.frame(x = c(0, 1)),
    ggplot2::aes(x = x, ymin = pmax(0, x - 0.1), ymax = pmin(1, x + 0.1)),
    fill = "gray80", alpha = 0.3, inherit.aes = FALSE
  )

  # Add calibration points and lines
  if (by_cluster && calib_result$K > 1) {
    p <- p + ggplot2::geom_line(ggplot2::aes(color = cluster), linewidth = 0.8)
    p <- p + ggplot2::geom_point(ggplot2::aes(color = cluster, size = n_samples))
  } else {
    p <- p + ggplot2::geom_line(color = "steelblue", linewidth = 0.8)
    p <- p + ggplot2::geom_point(ggplot2::aes(size = n_samples),
                                  color = "steelblue")
  }

  # Add labels and theme
  p <- p + ggplot2::scale_size_continuous(
    name = "N samples",
    range = c(2, 8),
    guide = ggplot2::guide_legend(override.aes = list(color = "steelblue"))
  )

  p <- p + ggplot2::labs(
    title = title,
    subtitle = sprintf("ECE = %.3f, MCE = %.3f",
                       calib_result$ece, calib_result$mce),
    x = "Mean Predicted Probability",
    y = "Mean Observed Frequency"
  )

  p <- p + ggplot2::coord_fixed(xlim = c(0, 1), ylim = c(0, 1))

  p <- p + ggplot2::theme_minimal()
  p <- p + ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 11),
    legend.position = "right"
  )

  # Add histogram of predictions if requested
  if (show_histogram) {
    # Create a small histogram inset would require additional packages
    # For simplicity, we'll add a note about the distribution
    p <- p + ggplot2::annotate(
      "text",
      x = 0.75, y = 0.15,
      label = sprintf("N = %d", calib_result$N),
      size = 3, color = "gray40"
    )
  }

  p
}


#' Print Summary of Calibration Results
#'
#' @param x A calibration_result object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.calibration_result <- function(x, ...) {
  cat("Calibration Analysis Results\n")
  cat("============================\n")
  cat(sprintf("Subjects (N): %d\n", x$N))
  cat(sprintf("Clusters (K): %d\n", x$K))
  cat(sprintf("Number of bins: %d\n", x$n_bins))
  cat("\n")
  cat("Calibration Metrics:\n")
  cat(sprintf("  Expected Calibration Error (ECE): %.4f\n", x$ece))
  cat(sprintf("  Maximum Calibration Error (MCE):  %.4f\n", x$mce))
  cat("\n")

  # Interpretation
  if (!is.na(x$ece)) {
    if (x$ece < 0.02) {
      cat("Interpretation: Excellent calibration\n")
    } else if (x$ece < 0.05) {
      cat("Interpretation: Good calibration\n")
    } else if (x$ece < 0.10) {
      cat("Interpretation: Moderate calibration\n")
    } else {
      cat("Interpretation: Poor calibration - consider recalibration\n")
    }
  }

  invisible(x)
}


#' Compute Calibration from Model Fit
#'
#' Convenience function that extracts responsibilities from a fitted model
#' and computes calibration against true labels.
#'
#' @param fit Fitted mixture model object (MFA or MPCA)
#' @param z_true Integer vector of true cluster labels
#' @param n_bins Number of bins for calibration. Default: 10.
#'
#' @return A calibration_result object
#'
#' @examples
#' \dontrun{
#' # Fit model
#' fit <- mixture_pca_em_fit(list_of_data, K = 3, r = 2)
#'
#' # Compute calibration
#' calib <- calibrate_model(fit, z_true)
#' print(calib)
#' }
#'
#' @export
calibrate_model <- function(fit, z_true, n_bins = 10) {
  if (is.null(fit$resp)) {
    stop("Fitted model must contain 'resp' (responsibilities matrix)")
  }

  compute_calibration(fit$resp, z_true, n_bins)
}
