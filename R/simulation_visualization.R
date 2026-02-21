#' Simulation Visualization Functions
#'
#' Functions for visualizing simulation results, including ARI plots, BIC heatmaps,
#' and VAF/reconstruction error visualizations.
#'
#' @name simulation_visualization
NULL

#' Read All Results from .rds Files
#'
#' This function searches a directory for \code{.rds} files, each presumably containing
#' one or more rows of simulation results, then merges them into a single data frame.
#' Each row includes a \code{file} column indicating the source file path.
#'
#' @param result_dir A character string specifying the directory to look for \code{.rds} files.
#'   Defaults to \code{"results"}.
#'
#' @return A combined \code{data.frame} with rows from all loaded files, including a \code{file} column.
#'
#' @details
#' If no files are found, the function raises an error. Otherwise, each file is read with
#' \code{\link{readRDS}} and merged using \code{dplyr::bind_rows}, which safely handles
#' heterogeneous column structures by filling missing columns with NA.
#'
#' @export
#' @examples
#' \dontrun{
#' df <- read_all_results("results")
#' head(df)
#' }
read_all_results <- function(result_dir = "results") {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Please install it.")
  }
  
  file_vec <- list.files(result_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(file_vec) == 0) {
    stop("No .rds files found in ", result_dir)
  }
  
  message("Reading ", length(file_vec), " .rds files from ", result_dir)
  
  df_list <- lapply(file_vec, function(f) {
    res <- tryCatch(readRDS(f), error = function(e) {
      warning("Failed to read file: ", f, " - ", e$message)
      return(NULL)
    })
    
    if (is.null(res)) {
      return(NULL)
    }
    
    # Handle both old format (data.frame) and new format (list with $summary)
    if (is.data.frame(res)) {
      res$file <- f
      return(res)
    } else if (is.list(res) && !is.null(res$summary)) {
      df <- res$summary
      df$file <- f
      return(df)
    } else {
      warning("Unrecognized format in file: ", f)
      return(NULL)
    }
  })
  
  # Remove NULL entries
  df_list <- df_list[!sapply(df_list, is.null)]
  
  if (length(df_list) == 0) {
    stop("No valid data frames found in .rds files")
  }
  
  # Use dplyr::bind_rows to safely merge data frames with different columns
  df_merged <- dplyr::bind_rows(df_list)
  message("Total results loaded: ", nrow(df_merged), " rows")
  
  df_merged
}

#' Create a Static Table of Simulation Conditions (Table 1)
#'
#' This function returns a small \code{data.frame} that describes the parameter
#' ranges used in the simulation (e.g., \code{K}, \code{r}, \code{N}, \code{M},
#' \code{T_each}, etc.).
#'
#' @return A \code{data.frame} with columns \code{Parameter} and \code{Values}.
#' @export
#' @examples
#' \dontrun{
#' table1 <- make_table1()
#' print(table1)
#' }
make_table1 <- function() {
  df_table1 <- data.frame(
    Parameter = c(
      "True cluster number (K)",
      "True factor (PC) number (r)",
      "Number of subjects (N)",
      "Channels (M)",
      "Time length (T_each)",
      "Separation parameter (sep)",
      "Noise scale (noise)"
    ),
    Values = c(
      "3, 4",
      "3, 5",
      "50, 100",
      "8, 12",
      "100",
      "0.5, 1.0",
      "1.0, 2.0"
    ),
    stringsAsFactors = FALSE
  )
  df_table1
}

#' Summarize Results for a Subset (Table 2)
#'
#' This function filters the data for certain separation/noise conditions
#' (e.g., \code{sep=1.0} and \code{noise=1.0}) and computes group-wise means and
#' standard deviations for ARI and BIC. It returns a small summary table.
#'
#' @param df A \code{data.frame} with columns like \code{ARI_MFA, ARI_PCA, BIC_MFA, BIC_PCA}.
#'
#' @return A \code{data.frame} summarizing ARI/BIC means and sds, grouped by
#'   \code{(K_true, r_true, N, M)}.
#'
#' @export
#' @examples
#' \dontrun{
#' df <- read_all_results("results")
#' table2 <- make_table2(df)
#' print(table2)
#' }
make_table2 <- function(df) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Please install it.")
  }
  
  df_subset <- df %>%
    dplyr::filter(abs(sep - 1.0) < 1e-9, abs(noise - 1.0) < 1e-9)

  df_table2 <- df_subset %>%
    dplyr::group_by(K_true, r_true, N, M) %>%
    dplyr::summarise(
      ARI_MFA_mean = mean(ARI_MFA),
      ARI_MFA_sd   = sd(ARI_MFA),
      BIC_MFA_mean = mean(BIC_MFA),
      ARI_PCA_mean = mean(ARI_PCA),
      ARI_PCA_sd   = sd(ARI_PCA),
      BIC_PCA_mean = mean(BIC_PCA),
      .groups = "drop"
    ) %>%
    dplyr::arrange(K_true, r_true, N, M)

  df_table2
}

#' Plot ARI vs. Separation/Noise (Figure 1)
#'
#' Creates a faceted plot showing ARI (Adjusted Rand Index) values for MFA and MPCA
#' across different separation and noise conditions.
#'
#' @param df A \code{data.frame} with columns \code{K_true, r_true, N, M, sep, noise, ARI_MFA, ARI_PCA}.
#' @return A ggplot object.
#' @export
#' @examples
#' \dontrun{
#' df <- read_all_results("results")
#' p <- plot_fig1_ARI(df)
#' print(p)
#' ggsave("figure1_ARI.pdf", p, width = 8, height = 6)
#' }
plot_fig1_ARI <- function(df) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Please install it.")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required. Please install it.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.")
  }
  
  df_long <- df %>%
    dplyr::select(K_true, r_true, N, M, sep, noise, ARI_MFA, ARI_PCA) %>%
    tidyr::pivot_longer(cols = c("ARI_MFA", "ARI_PCA"),
                        names_to = "Method", values_to = "ARI")

  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = factor(sep), y = ARI, 
                                              color = factor(noise),
                                              group = interaction(Method, noise))) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.3)) +
    ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.3)) +
    ggplot2::facet_wrap(~ Method + K_true + r_true + N + M, labeller = ggplot2::label_both) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Separation (sep)", y = "ARI", color = "Noise",
                  title = "(Figure 1) ARI vs sep/noise (MFA, MPCA)")

  return(p)
}

#' Plot BIC Heatmap for a Subset (Figure 2)
#'
#' Creates a heatmap of BIC values for a specific subset of conditions.
#'
#' @param df A \code{data.frame} with columns including \code{N, M, sep, noise, K_true, r_true, BIC_MFA}.
#' @return A ggplot object.
#' @export
#' @examples
#' \dontrun{
#' df <- read_all_results("results")
#' p <- plot_fig2_BIC_heatmap(df)
#' print(p)
#' ggsave("figure2_BIC.pdf", p, width = 6, height = 5)
#' }
plot_fig2_BIC_heatmap <- function(df) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Please install it.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.")
  }
  
  sub <- df %>%
    dplyr::filter(N == 50, M == 8, abs(sep - 1.0) < 1e-9, abs(noise - 1.0) < 1e-9)

  p <- ggplot2::ggplot(sub, ggplot2::aes(x = factor(r_true), y = factor(K_true), fill = BIC_MFA)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "r", y = "K", fill = "BIC",
                  title = "(Figure 2) BIC heatmap (MFA, N=50, M=8, sep=1.0, noise=1.0)")
  p
}

#' Plot Boxplot of ARI for Multiple Initial Seeds (Figure 3)
#'
#' Creates a boxplot showing the distribution of ARI values across multiple
#' random initializations.
#'
#' @param df A \code{data.frame} with columns \code{ARI_MFA, ARI_PCA}.
#' @return A ggplot object.
#' @export
#' @examples
#' \dontrun{
#' df <- read_all_results("results")
#' p <- plot_fig3_init_boxplot(df)
#' print(p)
#' ggsave("figure3_boxplot.pdf", p, width = 6, height = 4)
#' }
plot_fig3_init_boxplot <- function(df) {
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required. Please install it.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.")
  }
  
  df_long <- df %>%
    tidyr::pivot_longer(cols = c("ARI_MFA", "ARI_PCA"),
                        names_to = "Method", values_to = "ARI")

  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Method, y = ARI, color = Method)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "(Figure 3) Boxplot for multiple initial seeds")
  p
}

#' Plot Reconstruction Error and VAF (Figures 4-6)
#'
#' Creates visualizations of reconstruction error (SSE) and Variance Accounted For (VAF)
#' for different methods and conditions.
#'
#' @param df A \code{data.frame} with columns including VAF and SSE metrics.
#' @return A ggplot object.
#' @export
#' @examples
#' \dontrun{
#' df <- read_all_results("results")
#' p <- plot_fig4to6_reconstruction_vaf(df)
#' print(p)
#' ggsave("figure4to6_VAF.pdf", p, width = 10, height = 8)
#' }
plot_fig4to6_reconstruction_vaf <- function(df) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.")
  }
  
  # Check if VAF columns exist
  if (!all(c("VAF_MFA", "VAF_PCA") %in% colnames(df))) {
    warning("VAF columns not found in data frame. Returning empty plot.")
    return(ggplot2::ggplot() + ggplot2::theme_void() + 
           ggplot2::labs(title = "VAF data not available"))
  }
  
  # Create a simple VAF comparison plot
  df_long <- df %>%
    tidyr::pivot_longer(cols = c("VAF_MFA", "VAF_PCA"),
                        names_to = "Method", values_to = "VAF")
  
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = factor(sep), y = VAF, 
                                              color = Method,
                                              group = interaction(Method, noise))) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.3)) +
    ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.3)) +
    ggplot2::facet_wrap(~ K_true + r_true + N + M, labeller = ggplot2::label_both) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Separation", y = "VAF", color = "Method",
                  title = "(Figures 4-6) Variance Accounted For")
  
  p
}

#' Compute VAF for a Fitted MFA Model
#'
#' This function computes the Variance Accounted For (VAF) for a fitted MFA model
#' across all subjects, reconstructing each subject's data via factor scores and
#' summing squared errors.
#'
#' @param list_of_data A list of subject data matrices, each \code{(T_i x M)}.
#' @param mfa_fit An MFA fit object containing \code{Lambda}, \code{Psi}, \code{mu}, etc.
#' @param z An integer vector of length \code{N} specifying each subject's cluster ID (1..K).
#'
#' @return A numeric scalar representing the total VAF over all subjects.
#' @export
#' @examples
#' \dontrun{
#' sim <- generate_synergy_data(N = 50, K = 3, r = 2, M = 6, T_each = 100, seed = 123)
#' fit <- mfa_em_fit(sim$list_of_data, K = 3, r = 2)
#' vaf <- compute_VAF_mfa(sim$list_of_data, fit, fit$z)
#' print(vaf)
#' }
compute_VAF_mfa <- function(list_of_data, mfa_fit, z) {
  N <- length(list_of_data)
  K <- length(mfa_fit$Lambda)
  
  if (length(z) != N) stop("length(z) != N")

  total_num <- 0
  total_den <- 0

  for (i in seq_len(N)) {
    Xi <- list_of_data[[i]]
    k_i <- z[i]
    Lambda_k <- mfa_fit$Lambda[[k_i]]
    Psi_k    <- mfa_fit$Psi[[k_i]]
    mu_k     <- mfa_fit$mu[[k_i]]

    invPsi <- diag(1 / diag(Psi_k))

    M_ <- nrow(Lambda_k)
    r_ <- ncol(Lambda_k)
    A <- solve(diag(r_) + t(Lambda_k) %*% invPsi %*% Lambda_k)
    B <- t(Lambda_k) %*% invPsi

    X_centered <- sweep(Xi, 2, mu_k, FUN = "-")
    f_mat <- (X_centered %*% t(B)) %*% t(A)

    Xhat <- matrix(NA, nrow = nrow(Xi), ncol = ncol(Xi))
    for (t in seq_len(nrow(Xi))) {
      Xhat[t, ] <- mu_k + Lambda_k %*% f_mat[t, ]
    }

    ss_res <- sum((Xi - Xhat)^2)
    ss_total <- sum(Xi^2)

    total_num <- total_num + ss_res
    total_den <- total_den + ss_total
  }

  VAF <- 1 - total_num / total_den
  VAF
}

#' Compute VAF for a Mixture PCA (PPCA) Model
#'
#' This function computes the Variance Accounted For (VAF) for a fitted Mixture PCA model.
#'
#' @param list_of_data A list of subject data matrices.
#' @param mpca_fit A fitted MPCA model with \code{P, D, Psi, mu}.
#' @param z An integer vector of cluster assignments.
#'
#' @return A numeric scalar for the overall VAF.
#' @export
#' @examples
#' \dontrun{
#' sim <- generate_synergy_data(N = 50, K = 3, r = 2, M = 6, T_each = 100, seed = 123)
#' fit <- mixture_pca_em_fit(sim$list_of_data, K = 3, r = 2)
#' vaf <- compute_VAF_mpca(sim$list_of_data, fit, fit$z)
#' print(vaf)
#' }
compute_VAF_mpca <- function(list_of_data, mpca_fit, z) {
  N <- length(list_of_data)
  if (length(z) != N) {
    stop("length(z) != N (subject assignments do not match number of data sets).")
  }

  total_num <- 0
  total_den <- 0

  for (i in seq_len(N)) {
    Xi <- list_of_data[[i]]
    k_i <- z[i]

    P_k <- mpca_fit$P[[k_i]]
    D_k <- mpca_fit$D[[k_i]]
    Psi_k <- mpca_fit$Psi[[k_i]]
    mu_k  <- mpca_fit$mu[[k_i]]

    W_k <- P_k %*% D_k

    diagPsi <- diag(Psi_k)
    if (!all(abs(diff(diagPsi)) < 1e-10)) {
      # Non-isotropic
      invPsi <- diag(1 / diagPsi)
      M_ <- nrow(W_k)
      r_ <- ncol(W_k)
      X_centered <- sweep(Xi, 2, mu_k, FUN = "-")
      for (t in seq_len(nrow(Xi))) {
        x_t <- X_centered[t, ]
        tmp <- t(W_k) %*% invPsi %*% W_k + diag(r_)
        A <- solve(tmp)
        z_t <- A %*% (t(W_k) %*% invPsi %*% x_t)
        x_hat <- mu_k + W_k %*% z_t
        total_num <- total_num + sum((Xi[t, ] - x_hat)^2)
        total_den <- total_den + sum(Xi[t, ]^2)
      }
    } else {
      # Isotropic
      sigma2k <- diagPsi[1]
      M_ <- nrow(W_k)
      r_ <- ncol(W_k)
      A <- solve(t(W_k) %*% W_k + sigma2k * diag(r_))
      B <- t(W_k)

      X_centered <- sweep(Xi, 2, mu_k, FUN = "-")
      z_mat <- (X_centered %*% t(B)) %*% A

      Xhat <- (z_mat %*% t(W_k))
      for (t in seq_len(nrow(Xhat))) {
        Xhat[t, ] <- Xhat[t, ] + mu_k
      }

      ss_res <- sum((Xi - Xhat)^2)
      ss_total <- sum(Xi^2)

      total_num <- total_num + ss_res
      total_den <- total_den + ss_total
    }
  }

  VAF <- 1 - total_num / total_den
  VAF
}

#' Compile ARI and VAF into a Single Table
#'
#' This function takes a data frame with a \code{file} column (typically from
#' \code{\link{read_all_results}}), reads each result file to extract fitted models,
#' computes VAF for MFA and MPCA, and joins these back to the input data frame.
#'
#' @param df A data frame with a \code{file} column containing paths to \code{.rds} result files.
#'
#' @return The input data frame enriched with \code{VAF_MFA} and \code{VAF_PCA} columns.
#'   If models are not available in a file, VAF values will be \code{NA}.
#' @export
#' @examples
#' \dontrun{
#' df <- read_all_results("results")
#' df <- compute_ari_vaf_table(df)
#' head(df)
#' }
compute_ari_vaf_table <- function(df) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Please install it.")
  }
  
  # Check if file column exists
  if (is.null(df$file) || !"file" %in% colnames(df)) {
    message("Note: df$file column is missing. VAF computation skipped safely.")
    df$VAF_MFA <- NA
    df$VAF_PCA <- NA
    return(df)
  }
  
  message("Computing VAF for ", length(unique(df$file)), " result files...")
  
  # Compute VAF for each unique file
  vaf_list <- lapply(unique(df$file), function(f) {
    res_list <- tryCatch(readRDS(f), error = function(e) {
      warning("Failed to read file for VAF computation: ", f, " - ", e$message)
      return(NULL)
    })
    
    if (is.null(res_list)) {
      return(NULL)
    }
    
    # Initialize VAF values
    vaf_mfa <- NA
    vaf_pca <- NA
    
    # Compute VAF for MFA
    if (!is.null(res_list$best_model_mfa)) {
      vaf_mfa <- tryCatch({
        compute_VAF_mfa(
          list_of_data = res_list$best_model_mfa$list_of_data,
          mfa_fit      = res_list$best_model_mfa,
          z            = res_list$best_model_mfa$z
        )
      }, error = function(e) {
        warning("Failed to compute VAF_MFA for file: ", f, " - ", e$message)
        return(NA)
      })
    }
    
    # Compute VAF for MPCA
    if (!is.null(res_list$best_model_pca)) {
      vaf_pca <- tryCatch({
        compute_VAF_mpca(
          list_of_data = res_list$best_model_pca$list_of_data,
          mpca_fit     = res_list$best_model_pca,
          z            = res_list$best_model_pca$z
        )
      }, error = function(e) {
        warning("Failed to compute VAF_PCA for file: ", f, " - ", e$message)
        return(NA)
      })
    }
    
    data.frame(file = f, VAF_MFA = vaf_mfa, VAF_PCA = vaf_pca, stringsAsFactors = FALSE)
  })
  
  # Remove NULL entries
  vaf_list <- vaf_list[!sapply(vaf_list, is.null)]
  
  if (length(vaf_list) == 0) {
    message("Note: No models contained fit data. VAF computation skipped safely.")
    df$VAF_MFA <- NA
    df$VAF_PCA <- NA
    return(df)
  }
  
  # Combine VAF results
  vaf_df <- dplyr::bind_rows(vaf_list)
  
  # Join VAF back to original data frame by file
  df_out <- dplyr::left_join(df, vaf_df, by = "file")
  
  message("VAF computation complete.")
  df_out
}

#' Read All Results and Produce Tables/Figures
#'
#' This function reads all \code{.rds} result files in a specified directory, merges them
#' into a single data frame, prints summary tables, and generates several figures.
#'
#' @param result_dir A character string specifying the directory containing \code{.rds} result files.
#'   Defaults to \code{"results"}.
#'
#' @return Invisibly returns the merged data frame of all results.
#' @export
#' @examples
#' \dontrun{
#' df_all <- read_and_visualize_all("results")
#' head(df_all)
#' }
read_and_visualize_all <- function(result_dir = "results") {
  df <- read_all_results(result_dir = result_dir)

  df_table1 <- make_table1()
  message("=== Table 1: Simulation Conditions ===")
  message(paste(utils::capture.output(print(df_table1)), collapse = "\n"))

  df_table2 <- make_table2(df)
  message("=== Table 2: Summary of ARI etc. ===")
  message(paste(utils::capture.output(print(df_table2)), collapse = "\n"))

  p1 <- plot_fig1_ARI(df)
  ggplot2::ggsave("figure1_ARI_sep_noise.pdf", p1, width = 8, height = 6)

  p2 <- plot_fig2_BIC_heatmap(df)
  ggplot2::ggsave("figure2_BIC_heatmap.pdf", p2, width = 6, height = 5)

  if (nrow(df) > 1) {
    p3 <- plot_fig3_init_boxplot(df)
    ggplot2::ggsave("figure3_init_boxplot.pdf", p3, width = 6, height = 4)
  }

  invisible(df)
}

#' Plot Baseline Comparison (ARI and VAF)
#'
#' Creates a bar plot comparing ARI and VAF values across different methods,
#' including Two-step baselines. This visualization highlights the performance
#' differences between single-model approaches, two-step baselines, and mixture models.
#'
#' @param df A \code{data.frame} with columns \code{Method, ARI, VAF}.
#'   Typically the output from \code{simulate_and_compare_5methods}.
#' @param metric Character string specifying which metric to plot: "ARI" or "VAF".
#'   Defaults to "ARI".
#'
#' @return A ggplot object showing a bar plot of the specified metric by method.
#'
#' @details
#' The plot uses color coding to distinguish between method types:
#' \itemize{
#'   \item Single methods (SingleFA, SinglePCA) - shown in one color group
#'   \item Two-step methods (TwoStep_FA, TwoStep_PCA) - shown in another color group
#'   \item Mixture methods (MixtureFA, MixturePCA) - shown in a third color group
#' }
#'
#' This visualization is designed to demonstrate the superiority of mixture-based
#' approaches over simpler baseline methods.
#'
#' @examples
#' \dontrun{
#' df_comp <- simulate_and_compare_5methods(
#'   N = 100, K_true = 3, r_true = 4, M = 6, T_each = 100
#' )
#' p_ari <- plot_fig1_baseline_comparison(df_comp, metric = "ARI")
#' p_vaf <- plot_fig1_baseline_comparison(df_comp, metric = "VAF")
#' print(p_ari)
#' ggsave("baseline_comparison_ARI.pdf", p_ari, width = 8, height = 6)
#' }
#'
#' @export
plot_fig1_baseline_comparison <- function(df, metric = c("ARI", "VAF")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.")
  }
  
  metric <- match.arg(metric)
  
  # Check if the metric column exists
  if (!metric %in% colnames(df)) {
    stop(sprintf("Column '%s' not found in data frame", metric))
  }
  
  # Add method type for color grouping
  df$MethodType <- ifelse(grepl("^Single", df$Method), "Single",
                          ifelse(grepl("^TwoStep", df$Method), "TwoStep",
                                 ifelse(grepl("^Mixture", df$Method), "Mixture", "Other")))
  
  # Reorder methods for better visualization
  method_order <- c("SingleFA", "SinglePCA", "TwoStep_FA", "TwoStep_PCA", 
                    "MixtureFA", "MixturePCA")
  df$Method <- factor(df$Method, levels = method_order)
  
  # Create the plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Method, y = .data[[metric]], 
                                         fill = MethodType)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::scale_fill_manual(
      values = c("Single" = "#E69F00", "TwoStep" = "#56B4E9", "Mixture" = "#009E73"),
      name = "Method Type"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "top"
    ) +
    ggplot2::labs(
      x = "Method",
      y = metric,
      title = sprintf("Baseline Comparison: %s across Methods", metric),
      subtitle = "Comparing Single, Two-step, and Mixture approaches"
    )
  
  # Add a horizontal line at y=0 for reference
  if (metric == "ARI") {
    p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")
  }
  
  p
}
