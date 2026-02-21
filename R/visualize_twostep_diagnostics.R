#' Visualization Functions for Two-Step Algorithm Diagnostics
#'
#' This file contains visualization functions to diagnose failures in the
#' Two-step algorithm (fit_two_step_baseline). These functions help identify
#' issues with basis alignment, clustering, and reconstruction.
#'
#' @name visualize_twostep_diagnostics
NULL


#' Plot Aligned Basis Heatmaps
#'
#' Creates heatmaps of basis vectors (W matrices) for each subject and component,
#' showing how well the Two-step alignment procedure worked. Compares aligned
#' basis vectors with true synergies (if available) to diagnose alignment failures.
#'
#' @param W_list List of basis matrices, one per subject. Each matrix should be
#'   M x r (M channels/muscles, r components/synergies). This is typically the
#'   output from \code{fit_two_step_baseline()$W_list}.
#' @param W_true Optional list of true basis matrices for comparison (from simulation
#'   data, e.g., \code{sim$Lambda_list}). If provided, true synergies will be shown
#'   alongside estimated ones. Default is NULL.
#' @param z_true Optional vector of true cluster assignments (length N). Used to
#'   organize subjects by cluster. If NULL, subjects are shown in original order.
#' @param z_est Optional vector of estimated cluster assignments (length N). Used to
#'   show clustering results alongside basis vectors. Default is NULL.
#' @param max_subjects Maximum number of subjects to display. If N > max_subjects,
#'   a random sample will be shown. Default is 30.
#' @param theme_base_size Base font size for the plot theme. Default is 10.
#'
#' @return A ggplot object showing heatmaps of basis vectors organized by subject
#'   and component, with optional comparison to true synergies.
#'
#' @details
#' This function visualizes the output of the Two-step alignment procedure
#' (\code{align_basis_across_subjects()}) to diagnose alignment failures. Each
#' row represents one subject's basis vector for a specific component. Columns
#' represent muscles/channels. Colors indicate the weight/loading values.
#'
#' When \code{W_true} is provided, the function creates a side-by-side comparison
#' showing true synergies (left panel) and estimated aligned basis vectors (right
#' panel). This helps identify whether alignment failures are due to:
#' \itemize{
#'   \item Poor individual decomposition (Step 1)
#'   \item Alignment algorithm issues (Step 1.5)
#'   \item Fundamental differences between true and estimated structures
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate data with known ground truth
#' sim <- simulate_dynamic_synergy_data(N=30, K=3, r=3, M=8, T_each=200)
#' 
#' # Fit Two-step baseline
#' fit <- fit_two_step_baseline(
#'   sim$list_of_data, r=3, K=3, method="PCA",
#'   z_true=sim$z_true, align_basis=TRUE
#' )
#' 
#' # Visualize aligned basis vectors
#' plot_aligned_basis_heatmaps(fit$W_list, z_est=fit$z_est)
#' 
#' # Compare with true synergies
#' plot_aligned_basis_heatmaps(
#'   fit$W_list, 
#'   W_true=sim$Lambda_list,
#'   z_true=sim$z_true,
#'   z_est=fit$z_est
#' )
#' }
#'
#' @export
plot_aligned_basis_heatmaps <- function(W_list, W_true = NULL, z_true = NULL, 
                                        z_est = NULL, max_subjects = 30,
                                        theme_base_size = 10) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function.")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required for this function.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for this function.")
  }
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("Package 'viridis' is required for this function.")
  }
  
  N <- length(W_list)
  M <- nrow(W_list[[1]])
  r <- ncol(W_list[[1]])
  
  # Sample subjects if N is too large
  if (N > max_subjects) {
    set.seed(123)
    subject_indices <- sort(sample(seq_len(N), max_subjects))
    W_list <- W_list[subject_indices]
    if (!is.null(W_true)) W_true <- W_true[subject_indices]
    if (!is.null(z_true)) z_true <- z_true[subject_indices]
    if (!is.null(z_est)) z_est <- z_est[subject_indices]
    N <- max_subjects
  } else {
    subject_indices <- seq_len(N)
  }
  
  # Create data frame for estimated basis vectors
  W_df <- do.call(rbind, lapply(seq_len(N), function(i) {
    W_i <- W_list[[i]]
    data.frame(
      subject_id = i,
      subject_idx = subject_indices[i],
      component = rep(seq_len(r), each = M),
      muscle = rep(seq_len(M), r),
      value = as.vector(W_i),
      type = "Estimated"
    )
  }))
  
  # Add cluster information if available
  if (!is.null(z_est)) {
    W_df$cluster_est <- z_est[W_df$subject_id]
  }
  if (!is.null(z_true)) {
    W_df$cluster_true <- z_true[W_df$subject_id]
  }
  
  # If true basis is provided, add it to the data frame
  if (!is.null(W_true)) {
    W_true_df <- do.call(rbind, lapply(seq_len(N), function(i) {
      W_i <- W_true[[i]]
      data.frame(
        subject_id = i,
        subject_idx = subject_indices[i],
        component = rep(seq_len(r), each = M),
        muscle = rep(seq_len(M), r),
        value = as.vector(W_i),
        type = "True"
      )
    }))
    
    if (!is.null(z_est)) {
      W_true_df$cluster_est <- z_est[W_true_df$subject_id]
    }
    if (!is.null(z_true)) {
      W_true_df$cluster_true <- z_true[W_true_df$subject_id]
    }
    
    W_df <- rbind(W_true_df, W_df)
  }
  
  # Order subjects by true cluster if available, otherwise by estimated cluster
  if (!is.null(z_true)) {
    subject_order <- order(z_true[subject_indices])
    W_df$subject_id <- factor(W_df$subject_id, levels = subject_order)
  } else if (!is.null(z_est)) {
    subject_order <- order(z_est[subject_indices])
    W_df$subject_id <- factor(W_df$subject_id, levels = subject_order)
  } else {
    W_df$subject_id <- factor(W_df$subject_id, levels = seq_len(N))
  }
  
  W_df$component <- factor(W_df$component, levels = seq_len(r))
  W_df$muscle <- factor(W_df$muscle, levels = seq_len(M))
  W_df$type <- factor(W_df$type, levels = c("True", "Estimated"))
  
  # Create the plot
  p <- ggplot2::ggplot(W_df, ggplot2::aes(x = muscle, y = subject_id, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::facet_grid(component ~ type, 
                        scales = "free_y",
                        labeller = ggplot2::labeller(
                          component = function(x) paste0("Component ", x),
                          type = function(x) x
                        )) +
    viridis::scale_fill_viridis(option = "viridis", name = "Weight") +
    ggplot2::labs(
      title = "Aligned Basis Vectors Heatmap",
      subtitle = sprintf("N=%d subjects, r=%d components, M=%d muscles", N, r, M),
      x = "Muscle",
      y = "Subject"
    ) +
    ggplot2::theme_minimal(base_size = theme_base_size) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.6)),
      axis.text.x = ggplot2::element_text(size = ggplot2::rel(0.8)),
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = ggplot2::rel(0.9))
    )
  
  # Add cluster annotations if available
  if (!is.null(z_true) || !is.null(z_est)) {
    # Create annotation data
    if (!is.null(z_true)) {
      cluster_labels <- paste0("True: ", z_true[subject_indices][subject_order])
      if (!is.null(z_est)) {
        cluster_labels <- paste0(cluster_labels, " | Est: ", z_est[subject_indices][subject_order])
      }
    } else {
      cluster_labels <- paste0("Est: ", z_est[subject_indices][subject_order])
    }
    
    p <- p + ggplot2::labs(
      caption = "Subjects ordered by cluster assignment"
    )
  }
  
  return(p)
}


#' Plot Basis UMAP Embedding
#'
#' Creates a UMAP or PCA embedding of basis vectors at the component level
#' (not subject level) to visualize clustering quality. Each point represents
#' one basis vector (subject x component), colored by true cluster and shaped
#' by estimated cluster. This helps diagnose clustering failures in the Two-step
#' algorithm.
#'
#' @param W_list List of basis matrices, one per subject. Each matrix should be
#'   M x r (M channels/muscles, r components/synergies). This is typically the
#'   output from \code{fit_two_step_baseline()$W_list}.
#' @param z_true Vector of true cluster assignments (length N). Required for
#'   coloring points by true cluster.
#' @param z_est Optional vector of estimated cluster assignments (length N). If
#'   provided, points will be shaped by estimated cluster to show mismatches.
#'   Default is NULL.
#' @param method Character string specifying the embedding method: "umap" or
#'   "pca". Default is "umap".
#' @param theme_base_size Base font size for the plot theme. Default is 14.
#'
#' @return A ggplot object showing low-dimensional embedding of basis vectors
#'   (Nxr points), colored by true cluster and optionally shaped by estimated
#'   cluster.
#'
#' @details
#' This function extends \code{plot_cluster_embedding()} to work at the basis
#' vector level rather than the subject level. Instead of showing N points
#' (one per subject), it shows Nxr points (one per basis vector).
#'
#' This visualization is particularly useful for diagnosing Two-step clustering
#' failures because it reveals:
#' \itemize{
#'   \item Whether basis vectors from the same true cluster are close together
#'   \item Whether k-means clustering in Step 2 correctly separated clusters
#'   \item Which specific basis vectors were misclassified
#'   \item Whether alignment in Step 1.5 successfully matched corresponding components
#' }
#'
#' When \code{z_est} is provided, misclassified basis vectors are easy to spot
#' as points with different colors (true cluster) and shapes (estimated cluster).
#'
#' @examples
#' \dontrun{
#' # Simulate data with known ground truth
#' sim <- simulate_dynamic_synergy_data(N=30, K=3, r=3, M=8, T_each=200)
#' 
#' # Fit Two-step baseline
#' fit <- fit_two_step_baseline(
#'   sim$list_of_data, r=3, K=3, method="PCA",
#'   z_true=sim$z_true, align_basis=TRUE
#' )
#' 
#' # Visualize basis vector embedding
#' plot_basis_umap(fit$W_list, z_true=sim$z_true, z_est=fit$z_est)
#' 
#' # Use PCA instead of UMAP
#' plot_basis_umap(fit$W_list, z_true=sim$z_true, z_est=fit$z_est, method="pca")
#' }
#'
#' @export
plot_basis_umap <- function(W_list, z_true, z_est = NULL, method = "umap",
                            theme_base_size = 14) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function.")
  }
  
  N <- length(W_list)
  M <- nrow(W_list[[1]])
  r <- ncol(W_list[[1]])
  K <- length(unique(z_true))
  
  # Stack all basis vectors into a matrix (N*r rows, M columns)
  # Each row is one basis vector from one subject
  basis_mat <- do.call(rbind, lapply(seq_len(N), function(i) {
    W_list[[i]]  # M x r matrix, transpose to get r x M, then stack
  }))
  
  # Actually, we need to reshape properly: each basis vector is a column
  # So we want N*r rows, each with M features
  basis_mat <- do.call(rbind, lapply(seq_len(N), function(i) {
    W_i <- W_list[[i]]  # M x r
    # Each column is a basis vector, so transpose to get r rows of M features
    t(W_i)  # r x M
  }))
  
  # Create metadata for each basis vector
  basis_df <- data.frame(
    subject_id = rep(seq_len(N), each = r),
    component = rep(seq_len(r), N),
    cluster_true = rep(z_true, each = r)
  )
  
  if (!is.null(z_est)) {
    basis_df$cluster_est <- rep(z_est, each = r)
  }
  
  # Apply dimensionality reduction
  if (method == "umap") {
    if (!requireNamespace("uwot", quietly = TRUE)) {
      stop("Package 'uwot' is required for UMAP embedding. Install it with install.packages('uwot')")
    }
    
    # Set seed for reproducibility
    set.seed(123)
    
    # Use uwot for UMAP embedding
    n_neighbors <- min(15, nrow(basis_mat) - 1)
    embedding <- uwot::umap(basis_mat,
                           n_neighbors = n_neighbors,
                           min_dist = 0.1,
                           metric = "euclidean")
    embed_coords <- embedding
    colnames(embed_coords) <- c("Dim1", "Dim2")
    
    x_label <- "UMAP 1"
    y_label <- "UMAP 2"
    
  } else if (method == "pca") {
    pca_result <- stats::prcomp(basis_mat, center = TRUE, scale. = TRUE)
    embed_coords <- pca_result$x[, 1:2]
    colnames(embed_coords) <- c("Dim1", "Dim2")
    
    var_explained <- summary(pca_result)$importance[2, 1:2] * 100
    x_label <- sprintf("PC1 (%.1f%%)", var_explained[1])
    y_label <- sprintf("PC2 (%.1f%%)", var_explained[2])
    
  } else {
    stop("method must be either 'umap' or 'pca'")
  }
  
  # Add embedding coordinates to data frame
  basis_df$Dim1 <- embed_coords[, 1]
  basis_df$Dim2 <- embed_coords[, 2]
  
  # Convert to factors for plotting
  basis_df$cluster_true <- factor(basis_df$cluster_true, levels = seq_len(K))
  if (!is.null(z_est)) {
    basis_df$cluster_est <- factor(basis_df$cluster_est, levels = seq_len(K))
  }
  
  # Create the plot
  if (!is.null(z_est)) {
    # Show both true and estimated clusters
    p <- ggplot2::ggplot(basis_df, ggplot2::aes(x = Dim1, y = Dim2, 
                                                  color = cluster_true,
                                                  shape = cluster_est)) +
      ggplot2::geom_point(size = 3, alpha = 0.7) +
      viridis::scale_color_viridis(discrete = TRUE, option = "turbo", 
                                    name = "True\nCluster") +
      ggplot2::scale_shape_manual(
        values = c(16, 17, 15, 18, 3, 4, 8)[seq_len(K)],
        name = "Est.\nCluster"
      ) +
      ggplot2::labs(
        title = sprintf("Basis Vector Embedding (%s)", toupper(method)),
        subtitle = sprintf("N=%d subjects x r=%d components = %d basis vectors", N, r, N*r),
        x = x_label,
        y = y_label
      )
  } else {
    # Show only true clusters
    p <- ggplot2::ggplot(basis_df, ggplot2::aes(x = Dim1, y = Dim2, 
                                                  color = cluster_true)) +
      ggplot2::geom_point(size = 3, alpha = 0.7) +
      viridis::scale_color_viridis(discrete = TRUE, option = "turbo", 
                                    name = "True\nCluster") +
      ggplot2::labs(
        title = sprintf("Basis Vector Embedding (%s)", toupper(method)),
        subtitle = sprintf("N=%d subjects x r=%d components = %d basis vectors", N, r, N*r),
        x = x_label,
        y = y_label
      )
  }
  
  p <- p +
    ggplot2::theme_minimal(base_size = theme_base_size) +
    ggplot2::theme(
      legend.position = "right"
    )
  
  return(p)
}


#' Plot SSE Heatmap Across Methods and Subjects
#'
#' Creates a heatmap showing Sum of Squared Errors (SSE) for each subject
#' across different methods. This helps identify which methods perform poorly
#' for specific subjects and reveals the instability of the Two-step approach.
#'
#' @param results_list A named list of method results. Each element should be
#'   a fitted model object containing SSE information. Names should be method
#'   names (e.g., "TwoStep_PCA", "MixtureFA", "SinglePCA").
#' @param list_of_data List of data matrices, one per subject (T_i x M).
#'   Required to compute per-subject SSE if not already in results.
#' @param log_scale Logical. If TRUE, displays log10(SSE) instead of raw SSE.
#'   Recommended for better visualization when SSE values span multiple orders
#'   of magnitude. Default is TRUE.
#' @param theme_base_size Base font size for the plot theme. Default is 12.
#'
#' @return A ggplot object showing a heatmap of SSE values (subjects x methods).
#'
#' @details
#' This function computes per-subject reconstruction errors for multiple methods
#' and displays them as a heatmap. Darker colors indicate higher errors (worse
#' reconstruction). This visualization is useful for:
#' \itemize{
#'   \item Comparing reconstruction quality across methods
#'   \item Identifying subjects that are poorly reconstructed by Two-step methods
#'   \item Revealing instability in Two-step assignments (high SSE for misclassified subjects)
#'   \item Showing that mixture models provide more consistent reconstruction
#' }
#'
#' The function expects results from \code{fit_all_methods_on_data()} or similar
#' comparison functions. It computes per-subject SSE by reconstructing each
#' subject's data using the assigned cluster model.
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' sim <- simulate_dynamic_synergy_data(N=30, K=3, r=3, M=8, T_each=200)
#' 
#' # Fit all methods
#' results <- fit_all_methods_on_data(
#'   sim$list_of_data, z_true=sim$z_true,
#'   N=30, K=3, r=3, M=8, seed=123
#' )
#' 
#' # Create SSE heatmap
#' plot_sse_heatmap(results, sim$list_of_data)
#' 
#' # Use linear scale instead of log scale
#' plot_sse_heatmap(results, sim$list_of_data, log_scale=FALSE)
#' }
#'
#' @export
plot_sse_heatmap <- function(results_list, list_of_data, log_scale = TRUE,
                             theme_base_size = 12) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function.")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required for this function.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for this function.")
  }
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("Package 'viridis' is required for this function.")
  }
  
  N <- length(list_of_data)
  methods <- names(results_list)
  
  # Compute per-subject SSE for each method
  sse_matrix <- matrix(0, nrow = N, ncol = length(methods))
  colnames(sse_matrix) <- methods
  rownames(sse_matrix) <- paste0("S", seq_len(N))
  
  for (m in seq_along(methods)) {
    method_name <- methods[m]
    result <- results_list[[method_name]]
    
    # Compute per-subject SSE based on method type
    if (grepl("TwoStep", method_name) || grepl("Mixture", method_name)) {
      # Methods with cluster assignments
      z_vec <- result$z_est
      
      for (i in seq_len(N)) {
        X_i <- list_of_data[[i]]
        k_i <- z_vec[i]
        
        # Reconstruct based on method type
        if (grepl("PCA", method_name)) {
          # PCA-based reconstruction
          if (!is.null(result$cluster_models)) {
            # Two-step with cluster models
            cluster_model <- result$cluster_models[[k_i]]
            W_k <- cluster_model$W
            mu_i <- result$mu_list[[i]]
            sd_i <- result$sd_list[[i]]
            
            # Standardize
            X_std <- scale(X_i, center = mu_i, scale = sd_i)
            # Project and reconstruct
            scores <- X_std %*% W_k
            X_recon_std <- scores %*% t(W_k)
            # Unstandardize
            X_recon <- sweep(X_recon_std, 2, sd_i, "*")
            X_recon <- sweep(X_recon, 2, mu_i, "+")
          } else if (!is.null(result$W)) {
            # Mixture PCA
            W_k <- result$W[[k_i]]
            mu_k <- result$mu[[k_i]]
            X_centered <- sweep(X_i, 2, mu_k, "-")
            scores <- X_centered %*% W_k
            X_recon <- scores %*% t(W_k)
            X_recon <- sweep(X_recon, 2, mu_k, "+")
          }
        } else {
          # FA-based reconstruction
          if (!is.null(result$cluster_models)) {
            # Two-step FA (not typically used, but handle it)
            cluster_model <- result$cluster_models[[k_i]]
            Lambda_k <- cluster_model$W  # Using W as Lambda
            mu_i <- result$mu_list[[i]]
            sd_i <- result$sd_list[[i]]
            
            X_std <- scale(X_i, center = mu_i, scale = sd_i)
            scores <- X_std %*% Lambda_k %*% solve(t(Lambda_k) %*% Lambda_k)
            X_recon_std <- scores %*% t(Lambda_k)
            X_recon <- sweep(X_recon_std, 2, sd_i, "*")
            X_recon <- sweep(X_recon, 2, mu_i, "+")
          } else if (!is.null(result$Lambda)) {
            # Mixture FA
            Lambda_k <- result$Lambda[[k_i]]
            mu_k <- result$mu[[k_i]]
            Psi_k <- result$Psi[[k_i]]
            
            X_centered <- sweep(X_i, 2, mu_k, "-")
            # Factor scores: E[z|x] = (I + Lambda' Psi^-1 Lambda)^-1 Lambda' Psi^-1 x
            Psi_inv <- diag(1 / diag(Psi_k))
            M_inv <- solve(diag(ncol(Lambda_k)) + t(Lambda_k) %*% Psi_inv %*% Lambda_k)
            scores <- t(M_inv %*% t(Lambda_k) %*% Psi_inv %*% t(X_centered))
            X_recon <- scores %*% t(Lambda_k)
            X_recon <- sweep(X_recon, 2, mu_k, "+")
          }
        }
        
        sse_matrix[i, m] <- sum((X_i - X_recon)^2)
      }
      
    } else if (grepl("Single", method_name)) {
      # Single model methods (no clustering)
      for (i in seq_len(N)) {
        X_i <- list_of_data[[i]]
        
        if (grepl("PCA", method_name)) {
          # Single PCA
          P <- result$P
          mu <- result$mu
          X_centered <- sweep(X_i, 2, mu, "-")
          scores <- X_centered %*% P
          X_recon <- scores %*% t(P)
          X_recon <- sweep(X_recon, 2, mu, "+")
        } else {
          # Single FA
          Lambda <- result$loadings
          mu <- result$mu
          diagPsi <- result$diagPsi
          
          X_centered <- sweep(X_i, 2, mu, "-")
          Psi_inv <- diag(1 / diagPsi)
          M_inv <- solve(diag(ncol(Lambda)) + t(Lambda) %*% Psi_inv %*% Lambda)
          scores <- t(M_inv %*% t(Lambda) %*% Psi_inv %*% t(X_centered))
          X_recon <- scores %*% t(Lambda)
          X_recon <- sweep(X_recon, 2, mu, "+")
        }
        
        sse_matrix[i, m] <- sum((X_i - X_recon)^2)
      }
    }
  }
  
  # Convert to long format for ggplot
  sse_df <- as.data.frame(sse_matrix)
  sse_df$subject <- rownames(sse_matrix)
  sse_long <- tidyr::pivot_longer(sse_df, 
                                   cols = -subject,
                                   names_to = "method",
                                   values_to = "SSE")
  
  # Apply log scale if requested
  if (log_scale) {
    sse_long$SSE_plot <- log10(sse_long$SSE + 1e-10)  # Add small constant to avoid log(0)
    fill_label <- "log10(SSE)"
  } else {
    sse_long$SSE_plot <- sse_long$SSE
    fill_label <- "SSE"
  }
  
  # Order subjects and methods
  sse_long$subject <- factor(sse_long$subject, levels = paste0("S", seq_len(N)))
  sse_long$method <- factor(sse_long$method, levels = methods)
  
  # Create the plot
  p <- ggplot2::ggplot(sse_long, ggplot2::aes(x = method, y = subject, fill = SSE_plot)) +
    ggplot2::geom_tile() +
    viridis::scale_fill_viridis(option = "magma", name = fill_label) +
    ggplot2::labs(
      title = "Reconstruction Error (SSE) Heatmap",
      subtitle = sprintf("N=%d subjects across %d methods", N, length(methods)),
      x = "Method",
      y = "Subject"
    ) +
    ggplot2::theme_minimal(base_size = theme_base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.7)),
      panel.grid = ggplot2::element_blank()
    )
  
  return(p)
}


#' Plot Confusion Matrix for Cluster Assignments
#'
#' Creates a confusion matrix heatmap comparing true cluster assignments with
#' estimated cluster assignments. This is the standard way to visualize
#' clustering accuracy and identify systematic misclassification patterns.
#'
#' @param z_true Vector of true cluster assignments (length N, values 1 to K).
#' @param z_est Vector of estimated cluster assignments (length N, values 1 to K).
#' @param normalize Character string specifying normalization: "none" (raw counts),
#'   "true" (normalize by true cluster size), or "pred" (normalize by predicted
#'   cluster size). Default is "true".
#' @param theme_base_size Base font size for the plot theme. Default is 14.
#'
#' @return A ggplot object showing the confusion matrix as a heatmap.
#'
#' @details
#' The confusion matrix shows the correspondence between true and estimated
#' cluster assignments. Rows represent true clusters, columns represent
#' estimated clusters. Cell values indicate the number (or proportion) of
#' subjects with that true-estimated cluster pair.
#'
#' Perfect clustering produces a diagonal matrix (all subjects correctly
#' assigned). Off-diagonal entries indicate misclassifications. This
#' visualization complements the ARI metric by showing:
#' \itemize{
#'   \item Which clusters are confused with each other
#'   \item Whether misclassifications are systematic or random
#'   \item Cluster size imbalances
#'   \item Label permutation issues (if estimated clusters are numbered differently)
#' }
#'
#' Normalization options:
#' \itemize{
#'   \item \code{"none"}: Raw counts (useful for seeing absolute numbers)
#'   \item \code{"true"}: Divide by true cluster size (shows recall/sensitivity)
#'   \item \code{"pred"}: Divide by predicted cluster size (shows precision)
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' sim <- simulate_dynamic_synergy_data(N=30, K=3, r=3, M=8, T_each=200)
#' 
#' # Fit Two-step baseline
#' fit <- fit_two_step_baseline(
#'   sim$list_of_data, r=3, K=3, method="PCA",
#'   z_true=sim$z_true, align_basis=TRUE
#' )
#' 
#' # Plot confusion matrix
#' plot_confusion_matrix(sim$z_true, fit$z_est)
#' 
#' # Show raw counts instead of proportions
#' plot_confusion_matrix(sim$z_true, fit$z_est, normalize="none")
#' 
#' # Normalize by predicted cluster size
#' plot_confusion_matrix(sim$z_true, fit$z_est, normalize="pred")
#' }
#'
#' @export
plot_confusion_matrix <- function(z_true, z_est, normalize = "true",
                                  theme_base_size = 14) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function.")
  }
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("Package 'viridis' is required for this function.")
  }
  
  normalize <- match.arg(normalize, c("none", "true", "pred"))
  
  if (length(z_true) != length(z_est)) {
    stop("z_true and z_est must have the same length")
  }
  
  K_true <- length(unique(z_true))
  K_est <- length(unique(z_est))
  K <- max(K_true, K_est)
  
  # Create confusion matrix
  conf_mat <- table(True = z_true, Estimated = z_est)
  
  # Convert to matrix and ensure all clusters are represented
  conf_mat <- as.matrix(conf_mat)
  
  # Normalize if requested
  if (normalize == "true") {
    # Normalize by row (true cluster size)
    row_sums <- rowSums(conf_mat)
    conf_mat <- sweep(conf_mat, 1, row_sums, "/")
    value_label <- "Proportion"
  } else if (normalize == "pred") {
    # Normalize by column (predicted cluster size)
    col_sums <- colSums(conf_mat)
    conf_mat <- sweep(conf_mat, 2, col_sums, "/")
    value_label <- "Proportion"
  } else {
    value_label <- "Count"
  }
  
  # Convert to long format for ggplot
  conf_df <- as.data.frame(as.table(conf_mat))
  colnames(conf_df) <- c("True", "Estimated", "Value")
  
  # Ensure factors have all levels
  conf_df$True <- factor(conf_df$True, levels = seq_len(K_true))
  conf_df$Estimated <- factor(conf_df$Estimated, levels = seq_len(K_est))
  
  # Create the plot
  p <- ggplot2::ggplot(conf_df, ggplot2::aes(x = Estimated, y = True, fill = Value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(normalize == "none", 
                                                     sprintf("%d", round(Value)),
                                                     sprintf("%.2f", Value))),
                       size = theme_base_size / 3) +
    viridis::scale_fill_viridis(option = "plasma", name = value_label) +
    ggplot2::labs(
      title = "Confusion Matrix: True vs Estimated Clusters",
      subtitle = sprintf("ARI = %.3f", mclust::adjustedRandIndex(z_true, z_est)),
      x = "Estimated Cluster",
      y = "True Cluster"
    ) +
    ggplot2::theme_minimal(base_size = theme_base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1.0))
    ) +
    ggplot2::coord_fixed()
  
  return(p)
}


#' Plot Synergy Correlation Comparison
#'
#' Compares true synergies with estimated synergies by computing and visualizing
#' the correlation matrix between them. This extends the existing
#' \code{plot_synergy_correlation()} function to compare true vs estimated
#' synergies rather than comparing between clusters.
#'
#' @param W_true List of true basis matrices (e.g., \code{sim$Lambda_list}).
#'   Each matrix should be M x r.
#' @param W_est List of estimated basis matrices (e.g., from fitted model).
#'   Each matrix should be M x r.
#' @param z_true Vector of true cluster assignments (length N).
#' @param z_est Vector of estimated cluster assignments (length N).
#' @param cluster_id Integer specifying which cluster to analyze. If NULL,
#'   analyzes all clusters separately. Default is NULL.
#' @param method Correlation method: "pearson", "kendall", or "spearman".
#'   Default is "pearson".
#' @param theme_base_size Base font size for the plot theme. Default is 12.
#'
#' @return A ggplot object showing correlation heatmap(s) between true and
#'   estimated synergies.
#'
#' @details
#' This function computes the correlation between true synergy patterns
#' (from simulation data) and estimated synergy patterns (from fitted models).
#' For each cluster, it creates a correlation matrix where:
#' \itemize{
#'   \item Rows represent true synergies (components)
#'   \item Columns represent estimated synergies (components)
#'   \item Cell values show correlation between true component i and estimated component j
#' }
#'
#' Perfect recovery produces a diagonal matrix with high correlations (close to 1)
#' on the diagonal and low correlations off-diagonal. This indicates that each
#' estimated component correctly matches its corresponding true component.
#'
#' Common failure patterns:
#' \itemize{
#'   \item Permuted components: High correlations off-diagonal (alignment failure)
#'   \item Low correlations everywhere: Poor decomposition quality
#'   \item Mixed patterns: Some components recovered, others not
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' sim <- simulate_dynamic_synergy_data(N=30, K=3, r=3, M=8, T_each=200)
#' 
#' # Fit Two-step baseline
#' fit <- fit_two_step_baseline(
#'   sim$list_of_data, r=3, K=3, method="PCA",
#'   z_true=sim$z_true, align_basis=TRUE
#' )
#' 
#' # Compare true and estimated synergies
#' plot_synergy_correlation_comparison(
#'   sim$Lambda_list, fit$W_list,
#'   sim$z_true, fit$z_est
#' )
#' 
#' # Analyze specific cluster
#' plot_synergy_correlation_comparison(
#'   sim$Lambda_list, fit$W_list,
#'   sim$z_true, fit$z_est,
#'   cluster_id = 1
#' )
#' }
#'
#' @export
plot_synergy_correlation_comparison <- function(W_true, W_est, z_true, z_est,
                                                cluster_id = NULL,
                                                method = "pearson",
                                                theme_base_size = 12) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function.")
  }
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("Package 'viridis' is required for this function.")
  }
  
  method <- match.arg(method, c("pearson", "kendall", "spearman"))
  
  N <- length(W_true)
  M <- nrow(W_true[[1]])
  r <- ncol(W_true[[1]])
  K <- length(unique(z_true))
  
  # Determine which clusters to analyze
  if (is.null(cluster_id)) {
    clusters_to_plot <- seq_len(K)
  } else {
    clusters_to_plot <- cluster_id
  }
  
  # Compute correlation matrices for each cluster
  cor_list <- list()
  
  for (k in clusters_to_plot) {
    # Get subjects in this cluster (true and estimated)
    subjects_true <- which(z_true == k)
    subjects_est <- which(z_est == k)
    
    # Use subjects that are in both true and estimated cluster k
    subjects_both <- intersect(subjects_true, subjects_est)
    
    if (length(subjects_both) == 0) {
      warning(sprintf("No subjects in both true and estimated cluster %d, skipping", k))
      next
    }
    
    # Average basis vectors within cluster
    W_true_k <- Reduce("+", W_true[subjects_true]) / length(subjects_true)
    W_est_k <- Reduce("+", W_est[subjects_est]) / length(subjects_est)
    
    # Compute correlation matrix between true and estimated components
    cor_mat <- matrix(0, nrow = r, ncol = r)
    for (i in seq_len(r)) {
      for (j in seq_len(r)) {
        cor_mat[i, j] <- cor(W_true_k[, i], W_est_k[, j], method = method)
      }
    }
    
    rownames(cor_mat) <- paste0("True_", seq_len(r))
    colnames(cor_mat) <- paste0("Est_", seq_len(r))
    
    cor_list[[k]] <- cor_mat
  }
  
  # Convert to long format for ggplot
  cor_df_list <- lapply(seq_along(cor_list), function(idx) {
    k <- clusters_to_plot[idx]
    cor_mat <- cor_list[[idx]]
    
    if (is.null(cor_mat)) return(NULL)
    
    cor_df <- expand.grid(
      True_Component = rownames(cor_mat),
      Est_Component = colnames(cor_mat),
      stringsAsFactors = FALSE
    )
    cor_df$Correlation <- as.vector(cor_mat)
    cor_df$Cluster <- k
    cor_df
  })
  
  # Remove NULL entries
  cor_df_list <- cor_df_list[!sapply(cor_df_list, is.null)]
  
  if (length(cor_df_list) == 0) {
    stop("No valid clusters to plot")
  }
  
  cor_df <- do.call(rbind, cor_df_list)
  
  # Convert to factors
  cor_df$True_Component <- factor(cor_df$True_Component, 
                                   levels = paste0("True_", seq_len(r)))
  cor_df$Est_Component <- factor(cor_df$Est_Component,
                                  levels = paste0("Est_", seq_len(r)))
  cor_df$Cluster <- factor(cor_df$Cluster, levels = clusters_to_plot)
  
  # Create the plot
  p <- ggplot2::ggplot(cor_df, ggplot2::aes(x = Est_Component, y = True_Component, 
                                              fill = Correlation)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", Correlation)),
                       size = theme_base_size / 3.5) +
    ggplot2::scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, limits = c(-1, 1),
      name = paste0(stringr::str_to_title(method), "\nCorrelation")
    ) +
    ggplot2::labs(
      title = "True vs Estimated Synergy Correlation",
      subtitle = sprintf("%s correlation between true and estimated basis vectors", 
                        stringr::str_to_title(method)),
      x = "Estimated Component",
      y = "True Component"
    ) +
    ggplot2::theme_minimal(base_size = theme_base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::coord_fixed()
  
  # Add faceting if multiple clusters
  if (length(clusters_to_plot) > 1) {
    p <- p + ggplot2::facet_wrap(~ Cluster, nrow = 1,
                                  labeller = ggplot2::labeller(
                                    Cluster = function(x) paste0("Cluster ", x)
                                  ))
  }
  
  return(p)
}
