#' Generic Internal Helper for Plotting Cluster Synergy Loadings
#'
#' Internal function that provides shared logic for both MFA and MPCA
#' synergy loading plots.
#'
#' @param loading_list A list of loading matrices (M x r each)
#' @param model_type Character, either "MFA" or "MPCA"
#' @param cluster_ids A numeric vector of cluster IDs to plot
#' @param plot_type One of "heatmap" or "bar"
#' @param variable_names Optional character vector of variable/muscle names
#' @param variable_label Label for the variable axis (e.g., "Muscle" or "Variable")
#'
#' @return A ggplot object
#' @keywords internal
.plot_cluster_synergy_loadings_generic <- function(loading_list,
                                                    model_type = c("MFA", "MPCA"),
                                                    cluster_ids = NULL,
                                                    plot_type = c("heatmap", "bar"),
                                                    variable_names = NULL,
                                                    variable_label = "Muscle") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function. Please install it.")
  }
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("Package 'reshape2' is required for this function. Please install it.")
  }

  model_type <- match.arg(model_type)
  plot_type <- match.arg(plot_type)

  K <- length(loading_list)

  if (is.null(cluster_ids)) {
    cluster_ids <- seq_len(K)
  }

  # Collect all loadings into one data.frame for ggplot
  df_all <- list()
  for (k in cluster_ids) {
    L_k <- loading_list[[k]]
    if (is.null(L_k)) next  # skip empty clusters
    M <- nrow(L_k)
    r <- ncol(L_k)

    df_L <- data.frame(L_k)
    if (!is.null(variable_names) && length(variable_names) == M) {
      df_L$Variable <- factor(variable_names, levels = variable_names)
    } else {
      df_L$Variable <- factor(seq_len(M))
    }
    df_L <- reshape2::melt(df_L, id.vars = "Variable",
                           variable.name = "Factor", value.name = "Loading")
    df_L$Cluster <- factor(k)
    df_all[[k]] <- df_L
  }
  df_plot <- do.call(rbind, df_all)

  # Create plot title
  title_heatmap <- sprintf("%s Synergy Loadings (by Cluster)", model_type)
  title_bar <- sprintf("%s Synergy Loadings (by Cluster, Factor)", model_type)

  p <- NULL
  if (plot_type == "heatmap") {
    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Factor, y = Variable, fill = Loading)) +
      ggplot2::geom_tile() +
      scale_fill_mixsynergy_diverging() +
      ggplot2::facet_grid(. ~ Cluster, scales = "free_x") +
      ggplot2::labs(title = title_heatmap,
                    x = "Factor", y = variable_label) +
      theme_mixsynergy_heatmap()

  } else if (plot_type == "bar") {
    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Variable, y = Loading, fill = Factor)) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::facet_grid(Cluster ~ Factor, scales = "free_y") +
      ggplot2::labs(title = title_bar,
                    x = variable_label, y = "Loading") +
      theme_mixsynergy()
  }

  return(p)
}


#' Plot Cluster-Specific Synergy Coefficients (Factor Loadings) for MFA
#'
#' This function takes a fitted MFA model (e.g., from \code{mfa_em_fit()}) and
#' produces a faceted heatmap or bar plot of each cluster's factor loadings.
#'
#' @param mfa_fit A fitted MFA model object containing \code{Lambda}, \code{z}, etc.
#' @param cluster_ids A numeric vector of cluster IDs to plot. Defaults to all.
#' @param plot_type One of \code{"heatmap"} or \code{"bar"}. Specifies how to visualize.
#' @param muscle_names Optional character vector of muscle names for y-axis labels.
#'   If NULL, default channel names will be used.
#'
#' @return A \code{ggplot} object (which you can print or further modify).
#' @examples
#' \dontrun{
#'   # Suppose we have a best MFA model in 'best_mfa_model'
#'   p <- plot_cluster_synergy_loadings_mfa(best_mfa_model, cluster_ids = 1:2)
#'   print(p)
#' }
#' @export
plot_cluster_synergy_loadings_mfa <- function(mfa_fit,
                                              cluster_ids = NULL,
                                              plot_type = c("heatmap", "bar"),
                                              muscle_names = NULL) {
  .plot_cluster_synergy_loadings_generic(
    loading_list = mfa_fit$Lambda,
    model_type = "MFA",
    cluster_ids = cluster_ids,
    plot_type = plot_type,
    variable_names = muscle_names,
    variable_label = "Muscle"
  )
}

#' Plot All Factor Scores (Z) for MFA, Faceted by Cluster and Factor
#'
#' This function computes the factor scores Z_i(t) for each subject i,
#' identifies that subject's cluster k, and then creates a large facet plot
#' of the time evolution of each factor (1..r) in each cluster (1..K).
#' Optionally overlays all subjects in the same cluster.
#'
#' @param mfa_fit A fitted MFA model (with z, Lambda, Psi, mu, etc.).
#' @param list_of_data A list of (T_i x M) data matrices.
#' @param overlay_subjects Logical; if TRUE, all subjects in the same cluster
#'   are overlaid in the same facet. If FALSE, you could average them, etc.
#' @param muscle_names Optional character vector of muscle names for y-axis labels.
#'   If NULL, default channel names will be used.
#'
#' @return A ggplot object that shows the factor scores over time,
#'   with facets by (Cluster, Factor).
#'
#' @examples
#' \dontrun{
#'   p <- plot_all_factor_scores_mfa(my_mfa_fit, list_of_data)
#'   print(p)
#' }
#' @export
plot_all_factor_scores_mfa <- function(mfa_fit, list_of_data, overlay_subjects=TRUE, muscle_names = NULL) {

  # 必要: compute_factor_scores() という関数がある前提
  # (X_i, mu_k, Lambda_k, Psi_k) => Z_i(t)  (T_i x r)
  # 例: EZ = compute_factor_scores(X, mu, Lambda, Psi)

  N <- length(list_of_data)
  z_vec <- mfa_fit$z
  df_all <- list()

  for (i in seq_len(N)) {
    k_i <- z_vec[i]  # subject i のクラスタ
    X_i <- list_of_data[[i]]
    T_i <- nrow(X_i)

    # クラスタ k_i のパラメータ
    Lambda_k <- mfa_fit$Lambda[[k_i]]
    Psi_k    <- mfa_fit$Psi[[k_i]]
    mu_k     <- mfa_fit$mu[[k_i]]

    # 因子スコア Z_i(t,f)
    EZ_i <- compute_factor_scores(X_i, mu_k, Lambda_k, Psi_k)  # (T_i x r)

    # データフレーム化
    df_i <- as.data.frame(EZ_i)
    r <- ncol(EZ_i)
    # 列名を Factor1...r に
    factor_names <- paste0("Factor", seq_len(r))
    names(df_i) <- factor_names
    df_i$Time <- seq_len(T_i)
    df_i$Subject <- i
    df_i$Cluster <- k_i

    df_all[[i]] <- df_i
  }

  # 結合
  df_long <- do.call(rbind, df_all) %>%
    tidyr::pivot_longer(cols = tidyselect::starts_with("Factor"),
                 names_to = "Factor",
                 values_to = "Score")

  if (overlay_subjects) {
    p <- ggplot2::ggplot(df_long, ggplot2::aes(x=Time, y=Score, group=Subject)) +
      ggplot2::geom_line(alpha=0.5) +
      ggplot2::facet_grid(Cluster ~ Factor, scales="free_y") +
      ggplot2::labs(title="Factor Score Waveforms by Cluster and Factor",
           x="Time", y="Factor Score") +
      theme_mixsynergy()
  } else {
    df_avg <- df_long %>%
      dplyr::group_by(Cluster, Factor, Time) %>%
      dplyr::summarize(Score = mean(Score), .groups="drop")
    p <- ggplot2::ggplot(df_avg, ggplot2::aes(x=Time, y=Score)) +
      ggplot2::geom_line() +
      ggplot2::facet_grid(Cluster ~ Factor, scales="free_y") +
      ggplot2::labs(title="Average Factor Score Waveforms by Cluster and Factor",
           x="Time", y="Factor Score") +
      theme_mixsynergy()
  }

  return(p)
}

#' Plot Reconstructed Waveforms for All Clusters and Channels (MFA)
#'
#' This function plots the reconstructed waveforms (Xhat_i) for each subject
#' across all channels, grouped by cluster. By default, it overlays all subjects
#' in the same cluster, then uses facet to separate by (Cluster, Channel).
#'
#' @param mfa_fit A fitted MFA model (with z, Lambda, Psi, mu, etc.).
#' @param list_of_data A list of (T_i x M) data matrices used in fitting.
#' @param overlay_subjects Logical; if TRUE, each subject in the same cluster is drawn
#'   in the same facet. If FALSE, the function computes the average waveforms across
#'   subjects in each cluster.
#' @param muscle_names Optional character vector of muscle names for y-axis labels.
#'   If NULL, default channel names will be used.
#'
#' @return A ggplot object with line plots of reconstructed waveforms, facet by cluster & channel.
#' @examples
#' \dontrun{
#'   p_recon <- plot_all_reconstructed_waveforms_mfa(mfa_fit, list_of_data)
#'   print(p_recon)
#' }
#' @export
plot_all_reconstructed_waveforms_mfa <- function(mfa_fit, list_of_data, overlay_subjects = TRUE, muscle_names = NULL) {

  z_vec <- mfa_fit$z
  N <- length(list_of_data)
  df_list <- vector("list", N)

  for (i in seq_len(N)) {
    k_i <- z_vec[i]
    X_i <- list_of_data[[i]]
    T_i <- nrow(X_i)
    M_i <- ncol(X_i)

    Lambda_k <- mfa_fit$Lambda[[k_i]]
    Psi_k    <- mfa_fit$Psi[[k_i]]
    mu_k     <- mfa_fit$mu[[k_i]]

    # diagPsi_k
    diagPsi_k <- diag(Psi_k)
    # 再構成
    Xhat_i <- reconstruct_singleFA(X_i, mu_k, Lambda_k, diagPsi_k)

    # データフレーム化
    df_i <- as.data.frame(Xhat_i)
    # 列名を "Channel1"..."ChannelM" に
    if (!is.null(muscle_names) && length(muscle_names) == M_i) {
      channel_names <- muscle_names
    } else {
      channel_names <- paste0("Channel", seq_len(M_i))
    }
    colnames(df_i) <- channel_names
    df_i$Time <- seq_len(T_i)
    df_i$Subject <- i
    df_i$Cluster <- factor(k_i)

    df_list[[i]] <- df_i
  }

  # 結合
  df_plot <- do.call(rbind, df_list)
  # ロング形式へ
  df_long <- df_plot %>%
    tidyr::pivot_longer(cols = tidyselect::starts_with("Channel"),
                 names_to = "Channel",
                 values_to = "Amplitude")

  if (overlay_subjects) {
    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Time, y = Amplitude, group = Subject)) +
      ggplot2::geom_line(alpha = 0.6) +
      ggplot2::facet_grid(Cluster ~ Channel, scales = "free_y") +
      ggplot2::labs(title="Reconstructed Waveforms (All Clusters & Channels)",
           x="Time", y="Amplitude") +
      theme_mixsynergy()
  } else {
    df_avg <- df_long %>%
      dplyr::group_by(Cluster, Channel, Time) %>%
      dplyr::summarize(Amplitude = mean(Amplitude), .groups="drop")
    p <- ggplot2::ggplot(df_avg, ggplot2::aes(x = Time, y = Amplitude)) +
      ggplot2::geom_line() +
      ggplot2::facet_grid(Cluster ~ Channel, scales = "free_y") +
      ggplot2::labs(title="Average Reconstructed Waveforms (All Clusters & Channels)",
           x="Time", y="Amplitude") +
      theme_mixsynergy()
  }

  return(p)
}
#' Plot Original vs. Reconstructed Waveform for an MFA Model
#'
#' @param mfa_fit A fitted MFA model with \code{Lambda, Psi, mu, z}.
#' @param list_of_data The original data list used to fit the model.
#' @param subject_id Index of the subject to plot.
#' @param channel Which channel (muscle) to plot. Defaults to 1.
#' @param muscle_names Optional character vector of muscle names for y-axis labels.
#'   If NULL, default channel names will be used.
#'
#' @return A \code{ggplot} object showing original and reconstructed waveforms.
#' @export
plot_waveform_comparison_mfa <- function(mfa_fit,
                                         list_of_data,
                                         subject_id,
                                         channel = 1,
                                         muscle_names = NULL) {

  # Identify the cluster for this subject
  z_vec <- mfa_fit$z
  k_i <- z_vec[subject_id]

  # Extract parameters
  Lambda_k <- mfa_fit$Lambda[[k_i]]
  Psi_k    <- mfa_fit$Psi[[k_i]]
  mu_k     <- mfa_fit$mu[[k_i]]

  # Original data
  X_i <- list_of_data[[subject_id]]
  T_i <- nrow(X_i)

  # Reconstruct the data using the synergy-based function
  diagPsi_k <- diag(Psi_k)
  Xhat_i <- reconstruct_singleFA(X_i, mu_k, Lambda_k, diagPsi_k)

  # Prepare for ggplot
  df_plot <- data.frame(
    Time = seq_len(T_i),
    Original = X_i[, channel],
    Reconstructed = Xhat_i[, channel]
  )
  df_long <- tidyr::pivot_longer(df_plot,
                                 cols = c("Original","Reconstructed"),
                                 names_to = "Series",
                                 values_to = "Amplitude")

  channel_label <- if (!is.null(muscle_names) && length(muscle_names) >= channel) {
    muscle_names[channel]
  } else {
    paste0("Channel ", channel)
  }
  
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x=Time, y=Amplitude, color=Series)) +
    ggplot2::geom_line() +
    ggplot2::labs(title=paste0("Subject ", subject_id,
                      ", Cluster ", k_i,
                      ", ", channel_label),
         x="Time", y="Amplitude") +
    theme_mixsynergy()

  return(p)
}


#' Summarize VAF and SSE for an MFA Model
#'
#' Computes total SSE and VAF across all subjects, and optionally by cluster.
#'
#' @param mfa_fit A fitted MFA model (with z, Lambda, Psi, mu, etc.).
#' @param list_of_data The data list (a list of (T_i x M) matrices).
#' @param by_cluster Logical; if TRUE, returns a data.frame with SSE/VAF per cluster.
#'
#' @return If \code{by_cluster=FALSE}, returns a named vector with \code{SSE} and \code{VAF}.
#' If \code{by_cluster=TRUE}, returns a data.frame with columns
#'   \code{Cluster, SSE, T, SSE_per_time, VAF}.
#'
#' @export
summarize_vaf_sse_mfa <- function(mfa_fit, list_of_data, by_cluster=FALSE) {

  # 全データに対する SSE と SST, VAF を計算
  SSE_total <- calc_reconstruction_error_mixtureFA(list_of_data, mfa_fit)
  SST_total <- calc_total_SST(list_of_data)
  VAF_total <- 1 - SSE_total / SST_total

  # クラスタ別でなく、全体の統合値だけの場合
  if (!by_cluster) {
    return(c(SSE = SSE_total, VAF = VAF_total))
  }

  # クラスタ別集計
  z_vec <- mfa_fit$z
  K <- length(mfa_fit$Lambda)
  out_list <- vector("list", K)

  for (k in seq_len(K)) {
    idx_k <- which(z_vec == k)
    list_k <- list_of_data[idx_k]

    # もし該当クラスタに属する被験者がいなければ NA
    if (length(list_k) == 0) {
      out_list[[k]] <- data.frame(
        Cluster = k, SSE = NA, T = 0, SSE_per_time = NA, VAF = NA
      )
      next
    }

    #-- ここがポイント: サブモデル内でクラスタ番号をすべて1に変える --#
    # クラスタk向けの "1クラスタだけ" のサブモデルを作る
    sub_model <- list(
      z      = rep(1, length(list_k)),         # 全サンプルが "1" として扱われる
      Lambda = list(mfa_fit$Lambda[[k]]),      # 要素数1だけのリスト
      Psi    = list(mfa_fit$Psi[[k]]),         # 同様
      mu     = list(mfa_fit$mu[[k]])           # 同様
      # pi等が必要なら適宜追加
    )

    # サブモデルを使って該当クラスタデータだけのSSEを計算
    SSE_k <- calc_reconstruction_error_mixtureFA(list_k, sub_model)

    # T_k: このクラスタに属する全サブジェクトの合計行数
    T_k <- sum(sapply(list_k, nrow))

    # クラスタkの部分SSTを計算
    bigX <- do.call(rbind, list_k)
    mu_global <- colMeans(bigX)
    SST_k <- sum((sweep(bigX, 2, mu_global)^2))

    # VAF_k
    VAF_k <- 1 - SSE_k / SST_k

    out_list[[k]] <- data.frame(
      Cluster = k,
      SSE = SSE_k,
      T = T_k,
      SSE_per_time = SSE_k / T_k,
      VAF = VAF_k
    )
  }

  df_out <- do.call(rbind, out_list)
  return(df_out)
}


#' Line Plot of BIC vs. K or r
#'
#' Takes the summary data.frame from \code{select_optimal_K_r_mfa()} (which has columns
#' \code{K, r, logLik, BIC}) and produces a simple line plot of BIC. You can choose whether
#' the x-axis is \code{K} or \code{r}, and color/line by the other.
#'
#' @param df_summary A data.frame with columns \code{K, r, BIC}.
#' @param x_axis One of \code{"K"} or \code{"r"}. That variable is placed on the x-axis.
#'
#' @return A \code{ggplot} object.
#' @export
plot_bic_line <- function(df_summary, x_axis = c("K","r")) {

  x_axis <- match.arg(x_axis)

  # Decide which to treat as the grouping/color variable:
  if (x_axis == "K") {
    # x=K, color=factor(r)
    p <- ggplot2::ggplot(df_summary, ggplot2::aes(x = K, y = BIC, color=factor(r))) +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::labs(color="r") +
      ggplot2::scale_x_continuous(breaks = sort(unique(df_summary$K)))
  } else {
    # x=r, color=factor(K)
    p <- ggplot2::ggplot(df_summary, ggplot2::aes(x = r, y = BIC, color=factor(K))) +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::labs(color="K") +
      ggplot2::scale_x_continuous(breaks = sort(unique(df_summary$r)))
  }

  p <- p + ggplot2::labs(title="BIC vs. (K, r)", x=x_axis, y="BIC") +
    theme_mixsynergy()

  return(p)
}

#' Plot Cluster-Specific Synergy Coefficients (Factor Loadings) for MPCA
#'
#' This function takes a fitted MPCA model (e.g. from \code{mixture_pca_em_fit()})
#' and produces a faceted heatmap or bar plot of each cluster's factor loadings
#' \code{W[[k]]} (an \code{M x r} matrix).
#'
#' @param mpca_fit A fitted MPCA model object containing \code{W}, \code{z}, etc.
#' @param cluster_ids A numeric vector of cluster IDs to plot. Defaults to all.
#' @param plot_type One of \code{"heatmap"} or \code{"bar"}. Specifies how to visualize.
#' @param muscle_names Optional character vector of muscle names for y-axis labels.
#'   If NULL, default channel names will be used.
#'
#' @return A \code{ggplot} object (which you can print or further modify).
#'
#' @examples
#' \dontrun{
#'   # Suppose we have a best MPCA model in 'best_mpca_model'
#'   # that has best_mpca_model$W[[k]] for each cluster
#'   p <- plot_cluster_synergy_loadings_mpca(best_mpca_model, cluster_ids = 1:2)
#'   print(p)
#' }
#' @export
plot_cluster_synergy_loadings_mpca <- function(mpca_fit,
                                               cluster_ids = NULL,
                                               plot_type = c("heatmap", "bar"),
                                               muscle_names = NULL) {
  .plot_cluster_synergy_loadings_generic(
    loading_list = mpca_fit$W,
    model_type = "MPCA",
    cluster_ids = cluster_ids,
    plot_type = plot_type,
    variable_names = muscle_names,
    variable_label = "Variable"
  )
}

#' Plot All Factor Scores (Z) for MPCA, Faceted by Cluster and Factor
#'
#' This function uses precomputed factor scores \code{mpca_fit$factor_scores[[i]]}
#' (each \code{T_i x r}) and plots them over time, grouped by cluster \code{z[i]}.
#' Optionally overlays all subjects in the same cluster on one facet, or can average them.
#'
#' @param mpca_fit A fitted MPCA model (with z, W, mu, etc.) plus
#'   \code{factor_scores[[i]]} if available.
#' @param list_of_data A list of (T_i x M) data matrices (just for dimension reference).
#' @param overlay_subjects Logical; if TRUE, all subjects in the same cluster
#'   are overlaid in the same facet. If FALSE, we can average them etc.
#' @param muscle_names Optional character vector of muscle names for y-axis labels.
#'   If NULL, default channel names will be used.
#'
#' @return A ggplot object that shows the factor scores over time,
#'   with facets by (Cluster, Factor).
#'
#' @examples
#' \dontrun{
#'   mpca_fit <- compute_factor_scores_mpca(mpca_fit, list_of_data)
#'   p <- plot_all_factor_scores_mpca(mpca_fit, list_of_data, overlay_subjects=TRUE)
#'   print(p)
#' }
#' @export
plot_all_factor_scores_mpca <- function(mpca_fit, list_of_data, overlay_subjects=TRUE, muscle_names = NULL) {

  if (is.null(mpca_fit$factor_scores)) {
    stop("mpca_fit$factor_scores not found. Please run compute_factor_scores_mpca() first.")
  }

  N <- length(list_of_data)
  z_vec <- mpca_fit$z
  df_all <- list()

  for (i in seq_len(N)) {
    k_i <- z_vec[i]
    # factor scores for subject i => (T_i x r)
    scores_i <- mpca_fit$factor_scores[[i]]
    if (is.null(scores_i)) next  # skip if missing

    T_i <- nrow(scores_i)
    r_  <- ncol(scores_i)

    df_i <- as.data.frame(scores_i)
    factor_names <- paste0("Factor", seq_len(r_))
    colnames(df_i) <- factor_names
    df_i$Time    <- seq_len(T_i)
    df_i$Subject <- i
    df_i$Cluster <- k_i

    df_all[[i]] <- df_i
  }

  df_long <- do.call(rbind, df_all) %>%
    tidyr::pivot_longer(cols = tidyselect::starts_with("Factor"),
                 names_to = "Factor",
                 values_to = "Score")

  if (overlay_subjects) {
    p <- ggplot2::ggplot(df_long, ggplot2::aes(x=Time, y=Score, group=Subject)) +
      ggplot2::geom_line(alpha=0.5) +
      ggplot2::facet_grid(Cluster ~ Factor, scales="free_y") +
      ggplot2::labs(title="MPCA Factor Score Waveforms by Cluster and Factor",
           x="Time", y="Score") +
      theme_mixsynergy()
  } else {
    df_avg <- df_long %>%
      dplyr::group_by(Cluster, Factor, Time) %>%
      dplyr::summarize(Score = mean(Score), .groups="drop")

    p <- ggplot2::ggplot(df_avg, ggplot2::aes(x=Time, y=Score)) +
      ggplot2::geom_line() +
      ggplot2::facet_grid(Cluster ~ Factor, scales="free_y") +
      ggplot2::labs(title="Average MPCA Factor Score Waveforms by Cluster and Factor",
           x="Time", y="Score") +
      theme_mixsynergy()
  }

  return(p)
}

#' Reconstruct a single subject's data from a Mixture PCA model
#'
#' @param X Original \code{(T x M)} data
#' @param mu An \code{(M)} mean vector
#' @param W An \code{(M x r)} loading matrix
#' @param Z A \code{(T x r)} factor score matrix
#'
#' @return A \code{(T x M)} reconstructed matrix Xhat
#'
#' @export
reconstruct_singleMPCA <- function(X, mu, W, Z) {
  # X : (T x M) original (just used for dimension reference if needed)
  # Z : (T x r)
  # W : (M x r)
  # mu: length M
  T_i <- nrow(X)
  # Xhat = Z * W^T + mu
  # but be mindful of dimension => (T x r) %*% (r x M) => (T x M)
  Xhat <- Z %*% t(W)
  # add mu
  for(i in seq_len(T_i)){
    Xhat[i,] <- Xhat[i,] + mu
  }
  return(Xhat)
}


#' Plot Reconstructed Waveforms for All Clusters and Channels (MPCA)
#'
#' This function plots the reconstructed waveforms (Xhat_i) for each subject
#' across all channels, grouped by cluster. By default, it overlays all subjects
#' in the same cluster, then uses facet to separate by (Cluster, Channel).
#'
#' @param mpca_fit A fitted MPCA model (with z, W, mu, etc.).
#'   Also \code{factor_scores[[i]]} must be available if we want to reconstruct Xhat.
#' @param list_of_data A list of (T_i x M) data matrices used in fitting.
#' @param overlay_subjects Logical; if TRUE, each subject in the same cluster is drawn
#'   in the same facet. If FALSE, the function computes the average waveforms across
#'   subjects in each cluster.
#' @param muscle_names Optional character vector of muscle names for y-axis labels.
#'   If NULL, default channel names will be used.
#'
#' @return A ggplot object with line plots of reconstructed waveforms,
#'   facet by cluster & channel.
#'
#' @examples
#' \dontrun{
#'   mpca_fit <- compute_factor_scores_mpca(mpca_fit, list_of_data)
#'   p_recon <- plot_all_reconstructed_waveforms_mpca(mpca_fit, list_of_data)
#'   print(p_recon)
#' }
#' @export
plot_all_reconstructed_waveforms_mpca <- function(mpca_fit, list_of_data, overlay_subjects = TRUE, muscle_names = NULL) {

  if (is.null(mpca_fit$factor_scores)) {
    stop("mpca_fit$factor_scores not found. Please run compute_factor_scores_mpca() first.")
  }

  z_vec <- mpca_fit$z
  N <- length(list_of_data)
  df_list <- vector("list", N)

  for (i in seq_len(N)) {
    k_i <- z_vec[i]
    X_i <- list_of_data[[i]]
    T_i <- nrow(X_i)
    M_i <- ncol(X_i)

    W_k   <- mpca_fit$W[[k_i]]
    mu_k  <- mpca_fit$mu[[k_i]]
    Z_i   <- mpca_fit$factor_scores[[i]]  # (T_i x r)

    # Reconstruct
    Xhat_i <- reconstruct_singleMPCA(X_i, mu_k, W_k, Z_i)

    # convert to data.frame
    df_i <- as.data.frame(Xhat_i)
    if (!is.null(muscle_names) && length(muscle_names) == M_i) {
      channel_names <- muscle_names
    } else {
      channel_names <- paste0("Channel", seq_len(M_i))
    }
    colnames(df_i) <- channel_names
    df_i$Time    <- seq_len(T_i)
    df_i$Subject <- i
    df_i$Cluster <- factor(k_i)

    df_list[[i]] <- df_i
  }

  df_plot <- do.call(rbind, df_list)
  df_long <- df_plot %>%
    tidyr::pivot_longer(cols = tidyselect::starts_with("Channel"),
                 names_to = "Channel",
                 values_to = "Amplitude")

  if (overlay_subjects) {
    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Time, y = Amplitude, group = Subject)) +
      ggplot2::geom_line(alpha = 0.6) +
      ggplot2::facet_grid(Cluster ~ Channel, scales = "free_y") +
      ggplot2::labs(title="Reconstructed Waveforms (All Clusters & Channels, MPCA)",
           x="Time", y="Amplitude") +
      theme_mixsynergy()
  } else {
    df_avg <- df_long %>%
      dplyr::group_by(Cluster, Channel, Time) %>%
      dplyr::summarize(Amplitude = mean(Amplitude), .groups="drop")
    p <- ggplot2::ggplot(df_avg, ggplot2::aes(x = Time, y = Amplitude)) +
      ggplot2::geom_line() +
      ggplot2::facet_grid(Cluster ~ Channel, scales = "free_y") +
      ggplot2::labs(title="Avg Reconstructed Waveforms (All Clusters & Channels, MPCA)",
           x="Time", y="Amplitude") +
      theme_mixsynergy()
  }

  return(p)
}

#' Plot Original vs. Reconstructed Waveform for an MPCA Model
#'
#' @param mpca_fit A fitted MPCA model with \code{W, mu, z}, plus factor_scores.
#' @param list_of_data The original data list used to fit the model.
#' @param subject_id Index of the subject to plot (1-based).
#' @param channel Which channel (1-based) to plot. Defaults to 1.
#' @param muscle_names Optional character vector of muscle names for y-axis labels.
#'   If NULL, default channel names will be used.
#'
#' @return A \code{ggplot} object showing original and reconstructed waveforms.
#'
#' @export
plot_waveform_comparison_mpca <- function(mpca_fit,
                                          list_of_data,
                                          subject_id,
                                          channel = 1,
                                          muscle_names = NULL) {

  if (is.null(mpca_fit$factor_scores)) {
    stop("mpca_fit$factor_scores not found. Please run compute_factor_scores_mpca() first.")
  }

  z_vec <- mpca_fit$z
  k_i <- z_vec[subject_id]

  # Extract parameters
  W_k   <- mpca_fit$W[[k_i]]
  mu_k  <- mpca_fit$mu[[k_i]]
  Z_i   <- mpca_fit$factor_scores[[subject_id]]

  # Original data
  X_i <- list_of_data[[subject_id]]
  T_i <- nrow(X_i)

  # Reconstruct
  Xhat_i <- reconstruct_singleMPCA(X_i, mu_k, W_k, Z_i)

  # Prepare for ggplot
  df_plot <- data.frame(
    Time = seq_len(T_i),
    Original = X_i[, channel],
    Reconstructed = Xhat_i[, channel]
  )
  df_long <- tidyr::pivot_longer(df_plot,
                                 cols = c("Original","Reconstructed"),
                                 names_to = "Series",
                                 values_to = "Amplitude")

  channel_label <- if (!is.null(muscle_names) && length(muscle_names) >= channel) {
    muscle_names[channel]
  } else {
    paste0("Channel ", channel)
  }
  
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x=Time, y=Amplitude, color=Series)) +
    ggplot2::geom_line() +
    ggplot2::labs(title=paste0("Subject ", subject_id,
                      ", Cluster ", k_i,
                      ", ", channel_label),
         x="Time", y="Amplitude") +
    theme_mixsynergy()

  return(p)
}
