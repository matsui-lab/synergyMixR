#' Plot Synergy Bases (Lambda Matrices)
#'
#' Creates heatmaps of synergy matrices (Lambda_k) for each cluster.
#' Each heatmap shows the muscle-synergy structure with muscles on the y-axis
#' and synergies on the x-axis.
#'
#' @param sim A list returned by \code{simulate_dynamic_synergy_data()}.
#'   Must contain \code{Lambda_list} component.
#' @param theme_base_size Base font size for the plot theme. Default is 14.
#'
#' @return A ggplot object showing heatmaps of synergy matrices per cluster.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_dynamic_synergy_data(N=30, K=3, r=3, M=8, T_each=200)
#' plot_synergy_bases(sim)
#' }
#'
#' @export
plot_synergy_bases <- function(sim, theme_base_size = 14) {
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
  
  if (is.null(sim$Lambda_list)) {
    stop("sim must contain Lambda_list component")
  }
  
  K <- length(sim$Lambda_list)
  M <- nrow(sim$Lambda_list[[1]])
  r <- ncol(sim$Lambda_list[[1]])
  
  # Convert Lambda matrices to long format
  lambda_df <- do.call(rbind, lapply(seq_len(K), function(k) {
    lambda_k <- sim$Lambda_list[[k]]
    data.frame(
      cluster_id = k,
      muscle = rep(seq_len(M), r),
      synergy = rep(seq_len(r), each = M),
      value = as.vector(lambda_k)
    )
  }))
  
  lambda_df$cluster_id <- factor(lambda_df$cluster_id, levels = seq_len(K))
  lambda_df$muscle <- factor(lambda_df$muscle, levels = seq_len(M))
  lambda_df$synergy <- factor(lambda_df$synergy, levels = seq_len(r))
  
  p <- ggplot2::ggplot(lambda_df, ggplot2::aes(x = synergy, y = muscle, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~ cluster_id, nrow = 1, labeller = ggplot2::labeller(
      cluster_id = function(x) paste0("Cluster ", x)
    )) +
    viridis::scale_fill_viridis(option = "viridis", name = "Weight") +
    ggplot2::labs(
      title = "Synergy Bases (Lambda Matrices)",
      x = "Synergy",
      y = "Muscle"
    ) +
    ggplot2::theme_minimal(base_size = theme_base_size) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.8)),
      panel.grid = ggplot2::element_blank()
    )
  
  return(p)
}


#' Plot Activation Patterns
#'
#' Visualizes temporal activation curves per cluster by showing representative
#' activation patterns from subjects in each cluster.
#'
#' @param sim A list returned by \code{simulate_dynamic_synergy_data()}.
#'   Must contain \code{C_list} and \code{z_true} components.
#' @param n_subjects_per_cluster Number of subjects to display per cluster.
#'   Default is 3.
#' @param theme_base_size Base font size for the plot theme. Default is 14.
#'
#' @return A ggplot object showing temporal activation curves per cluster.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_dynamic_synergy_data(N=30, K=3, r=3, M=8, T_each=200)
#' plot_activation_patterns(sim)
#' }
#'
#' @export
plot_activation_patterns <- function(sim, n_subjects_per_cluster = 3, 
                                     theme_base_size = 14) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function.")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required for this function.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for this function.")
  }
  
  if (is.null(sim$C_list) || is.null(sim$z_true)) {
    stop("sim must contain C_list and z_true components")
  }
  
  K <- length(unique(sim$z_true))
  r <- ncol(sim$C_list[[1]])
  T_each <- nrow(sim$C_list[[1]])
  
  # Select representative subjects from each cluster
  activation_df <- do.call(rbind, lapply(seq_len(K), function(k) {
    cluster_subjects <- which(sim$z_true == k)
    n_select <- min(n_subjects_per_cluster, length(cluster_subjects))
    selected_subjects <- cluster_subjects[seq_len(n_select)]
    
    do.call(rbind, lapply(selected_subjects, function(i) {
      C_i <- sim$C_list[[i]]
      data.frame(
        cluster_id = k,
        subject_id = i,
        time = rep(seq_len(T_each), r),
        synergy = rep(seq_len(r), each = T_each),
        activation = as.vector(C_i)
      )
    }))
  }))
  
  activation_df$cluster_id <- factor(activation_df$cluster_id, levels = seq_len(K))
  activation_df$synergy <- factor(activation_df$synergy, levels = seq_len(r))
  
  p <- ggplot2::ggplot(activation_df, ggplot2::aes(x = time, y = activation, 
                                                     color = synergy, 
                                                     group = interaction(subject_id, synergy))) +
    ggplot2::geom_line(alpha = 0.6, linewidth = 0.5) +
    ggplot2::facet_wrap(~ cluster_id, nrow = 1, labeller = ggplot2::labeller(
      cluster_id = function(x) paste0("Cluster ", x)
    )) +
    viridis::scale_color_viridis(discrete = TRUE, option = "plasma", name = "Synergy") +
    ggplot2::labs(
      title = "Temporal Activation Patterns",
      x = "Time",
      y = "Activation"
    ) +
    ggplot2::theme_minimal(base_size = theme_base_size)
  
  return(p)
}


#' Plot Reconstructed EMG Signals
#'
#' Visualizes EMG time-series per cluster and subject. When cluster_id is NULL,
#' shows all clusters. For each cluster, individual subject EMG traces are shown
#' as thin semi-transparent lines, with the cluster-mean waveform overlaid as a
#' bold black line.
#'
#' @param sim A list returned by \code{simulate_dynamic_synergy_data()}.
#'   Must contain \code{X_list} and \code{z_true} components.
#' @param cluster_id Integer specifying which cluster to plot. If NULL, plots
#'   all clusters. Default is NULL.
#' @param muscles Integer vector specifying which muscles to plot. If NULL,
#'   plots all muscles. Default is NULL.
#' @param mean_ci Logical; if TRUE, display +/-1 SD ribbons (not yet implemented).
#'   Default is FALSE.
#' @param theme_base_size Base font size for the plot theme. Default is 13.
#'
#' @return A ggplot object showing EMG time-series with individual subjects and
#'   cluster means.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_dynamic_synergy_data(N=30, K=3, r=3, M=8, T_each=200)
#' plot_reconstructed_emg(sim)
#' plot_reconstructed_emg(sim, cluster_id = 1)
#' plot_reconstructed_emg(sim, muscles = c(1, 2, 3))
#' }
#'
#' @export
plot_reconstructed_emg <- function(sim, cluster_id = NULL, muscles = NULL, 
                                   mean_ci = FALSE, theme_base_size = 13) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function.")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required for this function.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for this function.")
  }
  
  if (is.null(sim$X_list) || is.null(sim$z_true)) {
    stop("sim must contain X_list and z_true components")
  }
  
  K <- length(unique(sim$z_true))
  M <- ncol(sim$X_list[[1]])
  T_each <- nrow(sim$X_list[[1]])
  
  # Determine which clusters to plot
  if (is.null(cluster_id)) {
    clusters_to_plot <- seq_len(K)
  } else {
    clusters_to_plot <- cluster_id
  }
  
  # Determine which muscles to plot
  if (is.null(muscles)) {
    muscles_to_plot <- seq_len(M)
  } else {
    muscles_to_plot <- muscles
  }
  
  # Build data frame with individual subject traces
  emg_df <- do.call(rbind, lapply(clusters_to_plot, function(k) {
    cluster_subjects <- which(sim$z_true == k)
    
    do.call(rbind, lapply(cluster_subjects, function(i) {
      X_i <- sim$X_list[[i]]
      data.frame(
        cluster_id = k,
        subject_id = i,
        time = rep(seq_len(T_each), length(muscles_to_plot)),
        muscle = rep(muscles_to_plot, each = T_each),
        value = as.vector(X_i[, muscles_to_plot])
      )
    }))
  }))
  
  # Compute cluster means
  emg_mean_df <- emg_df %>%
    dplyr::group_by(cluster_id, time, muscle) %>%
    dplyr::summarise(mean_value = mean(value), .groups = "drop")
  
  emg_df$cluster_id <- factor(emg_df$cluster_id, levels = clusters_to_plot)
  emg_df$muscle <- factor(emg_df$muscle, levels = muscles_to_plot)
  emg_mean_df$cluster_id <- factor(emg_mean_df$cluster_id, levels = clusters_to_plot)
  emg_mean_df$muscle <- factor(emg_mean_df$muscle, levels = muscles_to_plot)
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = emg_df, 
                       ggplot2::aes(x = time, y = value, group = subject_id),
                       color = "gray40", alpha = 0.3, linewidth = 0.3) +
    ggplot2::geom_line(data = emg_mean_df,
                       ggplot2::aes(x = time, y = mean_value),
                       color = "black", linewidth = 1.2) +
    ggplot2::facet_grid(muscle ~ cluster_id, 
                        scales = "free_y",
                        labeller = ggplot2::labeller(
                          cluster_id = function(x) paste0("Cluster ", x),
                          muscle = function(x) paste0("M", x)
                        )) +
    ggplot2::labs(
      title = "Reconstructed EMG Signals",
      x = "Time",
      y = "EMG Amplitude"
    ) +
    ggplot2::theme_minimal(base_size = theme_base_size) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = ggplot2::rel(0.9))
    )
  
  return(p)
}


#' Plot Noise Profiles
#'
#' Visualizes the distribution of noise parameters across clusters, showing
#' how stability-adaptability varies between clusters.
#'
#' @param sim A list returned by \code{simulate_dynamic_synergy_data()}.
#'   Must contain \code{stability_params} component.
#' @param theme_base_size Base font size for the plot theme. Default is 14.
#'
#' @return A ggplot object showing noise parameter distributions per cluster.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_dynamic_synergy_data(N=30, K=3, r=3, M=8, T_each=200)
#' plot_noise_profiles(sim)
#' }
#'
#' @export
plot_noise_profiles <- function(sim, theme_base_size = 14) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function.")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required for this function.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for this function.")
  }
  
  if (is.null(sim$stability_params)) {
    stop("sim must contain stability_params component")
  }
  
  K <- length(sim$stability_params)
  
  # Extract noise parameters
  noise_df <- do.call(rbind, lapply(seq_len(K), function(k) {
    params <- sim$stability_params[[k]]
    data.frame(
      cluster_id = k,
      control_profile = params$control_profile,
      sigma_white = params$sigma_white,
      sigma_low = params$sigma_low,
      filter_length = params$filter_length
    )
  }))
  
  # Reshape to long format for plotting
  noise_long <- tidyr::pivot_longer(
    noise_df,
    cols = c(sigma_white, sigma_low, filter_length),
    names_to = "parameter",
    values_to = "value"
  )
  
  noise_long$cluster_id <- factor(noise_long$cluster_id, levels = seq_len(K))
  noise_long$parameter <- factor(
    noise_long$parameter,
    levels = c("sigma_white", "sigma_low", "filter_length"),
    labels = c("White Noise (sigma)", "Low-Freq Noise (sigma)", "Filter Length")
  )
  
  p <- ggplot2::ggplot(noise_long, ggplot2::aes(x = cluster_id, y = value, 
                                                  fill = cluster_id)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~ parameter, scales = "free_y", nrow = 1) +
    viridis::scale_fill_viridis(discrete = TRUE, option = "mako", 
                                 name = "Cluster") +
    ggplot2::labs(
      title = "Noise Profiles (Stability-Adaptability)",
      x = "Cluster",
      y = "Parameter Value"
    ) +
    ggplot2::theme_minimal(base_size = theme_base_size)
  
  return(p)
}


#' Plot Cluster Embedding
#'
#' Creates a low-dimensional embedding (UMAP or PCA) to visualize cluster
#' separability. Each point represents a subject, colored by true cluster
#' assignment.
#'
#' @param sim A list returned by \code{simulate_dynamic_synergy_data()}.
#'   Must contain \code{X_list} and \code{z_true} components.
#' @param fit Optional fitted model object (e.g., from MFA/MPCA). If provided,
#'   predicted cluster assignments can be overlaid. Default is NULL.
#' @param method Character string specifying the embedding method: "umap" or
#'   "pca". Default is "umap".
#' @param theme_base_size Base font size for the plot theme. Default is 14.
#'
#' @return A ggplot object showing low-dimensional embedding of subjects.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_dynamic_synergy_data(N=30, K=3, r=3, M=8, T_each=200)
#' plot_cluster_embedding(sim, method = "umap")
#' plot_cluster_embedding(sim, method = "pca")
#' }
#'
#' @export
plot_cluster_embedding <- function(sim, fit = NULL, method = "umap", 
                                   theme_base_size = 14) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function.")
  }
  
  if (is.null(sim$X_list) || is.null(sim$z_true)) {
    stop("sim must contain X_list and z_true components")
  }
  
  N <- length(sim$X_list)
  K <- length(unique(sim$z_true))
  
  # Extract features: compute mean and variance per subject
  feature_mat <- do.call(rbind, lapply(seq_len(N), function(i) {
    X_i <- sim$X_list[[i]]
    c(colMeans(X_i), apply(X_i, 2, stats::var))
  }))
  
  # Apply dimensionality reduction
  if (method == "umap") {
    if (!requireNamespace("uwot", quietly = TRUE)) {
      stop("Package 'uwot' is required for UMAP embedding. Install it with install.packages('uwot')")
    }
    
    # Set seed for reproducibility
    set.seed(123)
    
    # Use uwot for UMAP embedding
    embedding <- uwot::umap(feature_mat,
                           n_neighbors = min(15, N - 1),
                           min_dist = 0.1,
                           metric = "euclidean")
    embed_coords <- embedding
    colnames(embed_coords) <- c("Dim1", "Dim2")
    
    x_label <- "UMAP 1"
    y_label <- "UMAP 2"
    
  } else if (method == "pca") {
    pca_result <- stats::prcomp(feature_mat, center = TRUE, scale. = TRUE)
    embed_coords <- pca_result$x[, 1:2]
    colnames(embed_coords) <- c("Dim1", "Dim2")
    
    var_explained <- summary(pca_result)$importance[2, 1:2] * 100
    x_label <- sprintf("PC1 (%.1f%%)", var_explained[1])
    y_label <- sprintf("PC2 (%.1f%%)", var_explained[2])
    
  } else {
    stop("method must be either 'umap' or 'pca'")
  }
  
  # Create data frame for plotting
  embed_df <- data.frame(
    Dim1 = embed_coords[, 1],
    Dim2 = embed_coords[, 2],
    cluster_id = factor(sim$z_true, levels = seq_len(K))
  )
  
  p <- ggplot2::ggplot(embed_df, ggplot2::aes(x = Dim1, y = Dim2, 
                                                color = cluster_id)) +
    ggplot2::geom_point(size = 3, alpha = 0.7) +
    viridis::scale_color_viridis(discrete = TRUE, option = "turbo", 
                                  name = "Cluster") +
    ggplot2::labs(
      title = sprintf("Cluster Embedding (%s)", toupper(method)),
      x = x_label,
      y = y_label
    ) +
    ggplot2::theme_minimal(base_size = theme_base_size) +
    ggplot2::theme(
      legend.position = "right"
    )
  
  return(p)
}


#' Plot Neural Command Trajectories
#'
#' Visualizes time-varying neural command trajectories estimated from synergy
#' activation patterns. This function creates a low-dimensional representation
#' of the temporal evolution of neural control signals, revealing how motor
#' commands change over time within and across clusters.
#'
#' @param sim A list returned by \code{simulate_dynamic_synergy_data()}.
#'   Must contain \code{C_list} and \code{z_true} components.
#' @param method Character string specifying the embedding method: "pca" or
#'   "umap". Default is "pca".
#' @param color_by Character string specifying how to color points: "cluster"
#'   (by cluster assignment) or "time" (by temporal progression). Default is
#'   "cluster".
#' @param smooth Logical; if TRUE, connect time points with smooth paths using
#'   geom_path with alpha gradient. Default is TRUE.
#' @param theme_base_size Base font size for the plot theme. Default is 14.
#'
#' @return A ggplot object showing neural command trajectories in 2D space.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_dynamic_synergy_data(N=30, K=3, r=3, M=8, T_each=200)
#' plot_neural_command(sim, method = "pca", color_by = "cluster")
#' plot_neural_command(sim, method = "umap", color_by = "time", smooth = TRUE)
#' }
#'
#' @export
plot_neural_command <- function(sim, method = "pca", color_by = "cluster", 
                                smooth = TRUE, theme_base_size = 14) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function.")
  }
  
  if (is.null(sim$C_list) || is.null(sim$z_true)) {
    stop("sim must contain C_list and z_true components")
  }
  
  if (!method %in% c("pca", "umap")) {
    stop("method must be either 'pca' or 'umap'")
  }
  
  if (!color_by %in% c("cluster", "time")) {
    stop("color_by must be either 'cluster' or 'time'")
  }
  
  N <- length(sim$C_list)
  K <- length(unique(sim$z_true))
  r <- ncol(sim$C_list[[1]])
  T_each <- nrow(sim$C_list[[1]])
  
  # Extract and normalize activation patterns
  # Stack all time points from all subjects into a single matrix
  all_activations <- do.call(rbind, lapply(seq_len(N), function(i) {
    C_i <- sim$C_list[[i]]
    # Normalize each synergy component to [0, 1]
    C_norm <- apply(C_i, 2, function(col) {
      col_range <- max(col) - min(col)
      if (col_range > 1e-10) {
        (col - min(col)) / col_range
      } else {
        col
      }
    })
    C_norm
  }))
  
  # Apply dimensionality reduction
  if (method == "umap") {
    if (!requireNamespace("uwot", quietly = TRUE)) {
      stop("Package 'uwot' is required for UMAP embedding. Install it with install.packages('uwot')")
    }
    
    # Set seed for reproducibility
    set.seed(123)
    embedding <- uwot::umap(all_activations, 
                           n_neighbors = min(15, nrow(all_activations) - 1),
                           min_dist = 0.1,
                           metric = "euclidean")
    embed_coords <- embedding
    colnames(embed_coords) <- c("Dim1", "Dim2")
    
    x_label <- "UMAP 1"
    y_label <- "UMAP 2"
    
  } else if (method == "pca") {
    pca_result <- stats::prcomp(all_activations, center = TRUE, scale. = TRUE)
    embed_coords <- pca_result$x[, 1:2]
    colnames(embed_coords) <- c("Dim1", "Dim2")
    
    var_explained <- summary(pca_result)$importance[2, 1:2] * 100
    x_label <- sprintf("PC1 (%.1f%%)", var_explained[1])
    y_label <- sprintf("PC2 (%.1f%%)", var_explained[2])
    
  }
  
  # Create data frame with subject and time information
  neural_df <- data.frame(
    Dim1 = embed_coords[, 1],
    Dim2 = embed_coords[, 2],
    subject_id = rep(seq_len(N), each = T_each),
    time = rep(seq_len(T_each), N),
    cluster_id = rep(sim$z_true, each = T_each)
  )
  
  neural_df$cluster_id <- factor(neural_df$cluster_id, levels = seq_len(K))
  
  # Normalize time to [0, 1] for coloring
  neural_df$time_norm <- (neural_df$time - 1) / (T_each - 1)
  
  # Create plot based on color_by parameter
  if (color_by == "cluster") {
    if (smooth) {
      p <- ggplot2::ggplot(neural_df, ggplot2::aes(x = Dim1, y = Dim2, 
                                                     color = cluster_id,
                                                     group = subject_id)) +
        ggplot2::geom_path(alpha = 0.3, linewidth = 0.5) +
        ggplot2::geom_point(alpha = 0.5, size = 0.8) +
        viridis::scale_color_viridis(discrete = TRUE, option = "turbo", 
                                      name = "Cluster")
    } else {
      p <- ggplot2::ggplot(neural_df, ggplot2::aes(x = Dim1, y = Dim2, 
                                                     color = cluster_id)) +
        ggplot2::geom_point(alpha = 0.5, size = 1.5) +
        viridis::scale_color_viridis(discrete = TRUE, option = "turbo", 
                                      name = "Cluster")
    }
  } else {  # color_by == "time"
    if (smooth) {
      p <- ggplot2::ggplot(neural_df, ggplot2::aes(x = Dim1, y = Dim2, 
                                                     color = time_norm,
                                                     group = subject_id)) +
        ggplot2::geom_path(alpha = 0.3, linewidth = 0.5) +
        ggplot2::geom_point(alpha = 0.5, size = 0.8) +
        viridis::scale_color_viridis(option = "plasma", name = "Time\n(normalized)")
    } else {
      p <- ggplot2::ggplot(neural_df, ggplot2::aes(x = Dim1, y = Dim2, 
                                                     color = time_norm)) +
        ggplot2::geom_point(alpha = 0.5, size = 1.5) +
        viridis::scale_color_viridis(option = "plasma", name = "Time\n(normalized)")
    }
  }
  
  # Add faceting if coloring by time (to show clusters separately)
  if (color_by == "time") {
    p <- p + ggplot2::facet_wrap(~ cluster_id, nrow = 1, labeller = ggplot2::labeller(
      cluster_id = function(x) paste0("Cluster ", x)
    ))
  }
  
  p <- p +
    ggplot2::labs(
      title = sprintf("Neural Command Trajectories (%s)", toupper(method)),
      x = x_label,
      y = y_label
    ) +
    ggplot2::theme_minimal(base_size = theme_base_size) +
    ggplot2::theme(
      legend.position = "right"
    )
  
  return(p)
}


#' Plot Synergy Summary
#'
#' Creates a comprehensive summary plot combining multiple visualization panels
#' using patchwork. This provides an overview of spatial synergy structures,
#' temporal activation dynamics, and cluster separability.
#'
#' @param sim A list returned by \code{simulate_dynamic_synergy_data()}.
#' @param fit Optional fitted model object (e.g., from MFA/MPCA). Default is NULL.
#' @param theme_base_size Base font size for the plot theme. Default is 12.
#'
#' @return A patchwork object combining multiple ggplot panels.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_dynamic_synergy_data(N=30, K=3, r=3, M=8, T_each=200)
#' plot_synergy_summary(sim)
#' }
#'
#' @export
plot_synergy_summary <- function(sim, fit = NULL, theme_base_size = 12) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for this function. Install it with install.packages('patchwork')")
  }
  
  # Generate individual plots
  p1 <- plot_synergy_bases(sim, theme_base_size = theme_base_size)
  p2 <- plot_activation_patterns(sim, theme_base_size = theme_base_size)
  p3 <- plot_noise_profiles(sim, theme_base_size = theme_base_size)
  p4 <- plot_cluster_embedding(sim, fit = fit, method = "umap", 
                               theme_base_size = theme_base_size)
  
  # Combine plots using patchwork
  layout <- (p1 / p2) | (p3 / p4)
  
  combined <- layout + 
    patchwork::plot_annotation(
      title = "Synergy Data Summary",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = theme_base_size + 2, 
                                                                  face = "bold"))
    )
  
  return(combined)
}
