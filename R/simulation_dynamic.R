#' Simulate Dynamic Synergy Data with Three Orthogonal Control Axes
#'
#' This function extends \code{simulate_dynamic_synergy_data()} to systematically
#' evaluate how different motor control strategies affect muscle synergy structures.
#' The simulation model includes three orthogonal control axes:
#' (1) Spatial structure, (2) Temporal coordination, and (3) Stability-Adaptability.
#' These axes represent complementary aspects of motor control diversity observed
#' in both rehabilitation and sports science contexts.
#'
#' Each cluster (k = 1, ..., K) is assigned a unique combination of parameters
#' along these three axes, emulating distinct movement control strategies
#' (e.g., proximal-dominant, alternating coordination, or adaptive noisy control).
#'
#' @section Spatial Structure Axis:
#' The spatial axis defines how strongly each muscle group contributes to each
#' synergy component. For each cluster k, the synergy matrix (Lambda_k) is
#' generated with cluster-specific bias patterns (B_k) that shift the activation
#' balance among proximal, distal, or antagonistic muscle groups.
#'
#' The parameter \code{cluster_sep_spatial} controls the degree of differentiation
#' between clusters, ranging from 0 (identical synergy weights) to 2 (strongly
#' distinct spatial modules).
#'
#' @section Temporal Coordination Axis:
#' The temporal axis governs how synergies are activated over time, including
#' their onset timing, overlap, and rhythmic frequency. For each cluster, the
#' activation matrix C_i^(k)(t) is generated with cluster-specific frequencies
#' (f_k) and phase offsets (phi_ij^(k)).
#'
#' Three representative temporal coordination modes:
#' \itemize{
#'   \item Cluster 1: Co-contraction/rigid (all synergies in-phase, strong overlap)
#'   \item Cluster 2: Alternating/efficient (phases evenly distributed, minimal overlap)
#'   \item Cluster 3: Adaptive/reactive (random phases, heterogeneous frequencies)
#' }
#'
#' The parameter \code{cluster_sep_temporal} (0-2) scales the magnitude of
#' inter-cluster differences in phase and frequency.
#'
#' @section Stability-Adaptability Axis:
#' The stability axis introduces physiological and measurement variability into
#' the EMG signal. Each cluster k is assigned distinct noise parameters to
#' simulate differences in control precision and adaptability.
#'
#' Three representative control profiles:
#' \itemize{
#'   \item Cluster 1: Stable control (minimal noise, elite or late rehab)
#'   \item Cluster 2: Moderate control (balanced smoothness vs flexibility)
#'   \item Cluster 3: Adaptive/exploratory (variable, early rehab or fatigue)
#' }
#'
#' @param N Number of subjects (time series). Default is 30.
#' @param K Number of clusters (subgroups). Default is 3.
#' @param r Number of synergies per cluster. Default is 3.
#' @param M Number of observed channels (muscles). Default is 8.
#' @param T_each Time-series length per subject. Default is 200.
#' @param cluster_sep_spatial Numeric scalar (0-2) controlling spatial structure
#'   differentiation between clusters. Higher values create more distinct spatial
#'   modules. Default is 1.0.
#' @param cluster_sep_temporal Numeric scalar (0-2) controlling temporal coordination
#'   differentiation between clusters. Higher values create more distinct temporal
#'   patterns. Default is 1.0.
#' @param cluster_sep_stability Numeric scalar (0-2) controlling stability-adaptability
#'   differentiation between clusters. Higher values create more distinct noise
#'   profiles. Default is 1.0.
#' @param freq_range Numeric vector of length 2 specifying the range of base
#'   frequencies (in Hz) for the periodic activation signals. Default is c(0.5, 3.0).
#' @param amp_range Numeric vector of length 2 specifying the range of amplitudes
#'   for the periodic activation signals. Default is c(0.5, 1.5).
#' @param base_level Numeric scalar specifying the baseline activation level
#'   added to C(t) to avoid all-zero activations. Default is 0.05.
#' @param rectify_output Logical; if TRUE, apply half-wave rectification to
#'   X(t) to produce non-negative EMG-like signals. Default is TRUE.
#' @param seed Random seed for reproducibility. Default is 123.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{X_list}}{List of length N; each element is a (T_each x M) matrix
#'     of generated EMG signals.}
#'   \item{\code{C_list}}{List of length N; each element is a (T_each x r) matrix
#'     of synergy activation signals.}
#'   \item{\code{Lambda_list}}{List of length K; each element is a (M x r) matrix
#'     of cluster-specific synergy loadings (non-negative).}
#'   \item{\code{z_true}}{Integer vector of length N indicating true cluster
#'     assignments (1 to K).}
#'   \item{\code{spatial_params}}{List containing spatial structure parameters
#'     for each cluster (bias patterns).}
#'   \item{\code{temporal_params}}{List containing temporal coordination parameters
#'     for each cluster (frequencies, phase offsets).}
#'   \item{\code{stability_params}}{List containing stability-adaptability parameters
#'     for each cluster (noise levels, filter lengths).}
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' sim_v2 <- simulate_dynamic_synergy_data(
#'   N = 30, K = 3, r = 3, M = 8, T_each = 200
#' )
#' 
#' # High spatial separation, low temporal separation
#' sim_spatial <- simulate_dynamic_synergy_data(
#'   N = 30, K = 3, r = 3, M = 8, T_each = 200,
#'   cluster_sep_spatial = 2.0,
#'   cluster_sep_temporal = 0.5,
#'   cluster_sep_stability = 1.0,
#'   seed = 456
#' )
#' 
#' # Factorial design: 3x3x3 combinations
#' spatial_levels <- c(0.5, 1.0, 2.0)
#' temporal_levels <- c(0.5, 1.0, 2.0)
#' stability_levels <- c(0.5, 1.0, 2.0)
#' 
#' results <- list()
#' idx <- 1
#' for (sp in spatial_levels) {
#'   for (tp in temporal_levels) {
#'     for (st in stability_levels) {
#'       results[[idx]] <- simulate_dynamic_synergy_data(
#'         N = 30, K = 3, r = 3, M = 8, T_each = 200,
#'         cluster_sep_spatial = sp,
#'         cluster_sep_temporal = tp,
#'         cluster_sep_stability = st,
#'         seed = 100 + idx
#'       )
#'       idx <- idx + 1
#'     }
#'   }
#' }
#' }
#'
#' @export
simulate_dynamic_synergy_data <- function(
  N = 30,
  K = 3,
  r = 3,
  M = 8,
  T_each = 200,
  cluster_sep_spatial = 1.0,
  cluster_sep_temporal = 1.0,
  cluster_sep_stability = 1.0,
  freq_range = c(0.5, 3.0),
  amp_range = c(0.5, 1.5),
  base_level = 0.05,
  rectify_output = TRUE,
  seed = 123
) {
  set.seed(seed)
  
  # Input validation
  if (N <= 0 || K <= 0 || r <= 0 || M <= 0 || T_each <= 0) {
    stop("N, K, r, M, and T_each must be positive integers")
  }
  if (length(freq_range) != 2 || freq_range[1] >= freq_range[2]) {
    stop("freq_range must be a vector of length 2 with freq_range[1] < freq_range[2]")
  }
  if (length(amp_range) != 2 || amp_range[1] >= amp_range[2]) {
    stop("amp_range must be a vector of length 2 with amp_range[1] < amp_range[2]")
  }
  if (base_level < 0) {
    stop("base_level must be non-negative")
  }
  if (cluster_sep_spatial < 0 || cluster_sep_temporal < 0 || cluster_sep_stability < 0) {
    stop("cluster separation parameters must be non-negative")
  }
  
  # Cluster assignments (uniform probability)
  z_true <- sample.int(K, N, replace = TRUE)
  
  # ========================================================================
  # AXIS 1: SPATIAL STRUCTURE - Cluster-Specific Synergy Matrices (Lambda_k)
  # ========================================================================
  
  # Generate base synergy matrix
  Lambda_base <- matrix(runif(M * r, 0.3, 0.7), M, r)
  
  # Define cluster-specific bias patterns for spatial structure
  spatial_params <- vector("list", K)
  Lambda_list <- vector("list", K)
  
  for (k in seq_len(K)) {
    # Create cluster-specific bias pattern
    # Cluster 1: Proximal-dominant (emphasize first M/2 muscles)
    # Cluster 2: Distal-dominant (emphasize last M/2 muscles)
    # Cluster 3: Balanced/mixed pattern
    
    if (K == 1) {
      # Single cluster: no bias
      B_k <- matrix(0, M, r)
    } else if (K == 2) {
      # Two clusters: proximal vs distal
      if (k == 1) {
        # Proximal-dominant
        B_k <- matrix(0, M, r)
        B_k[1:ceiling(M/2), ] <- cluster_sep_spatial * 0.3
        B_k[(ceiling(M/2)+1):M, ] <- -cluster_sep_spatial * 0.2
      } else {
        # Distal-dominant
        B_k <- matrix(0, M, r)
        B_k[1:ceiling(M/2), ] <- -cluster_sep_spatial * 0.2
        B_k[(ceiling(M/2)+1):M, ] <- cluster_sep_spatial * 0.3
      }
    } else {
      # Three or more clusters
      if (k == 1) {
        # Proximal-dominant (hip and thigh)
        B_k <- matrix(0, M, r)
        B_k[1:ceiling(M/3), ] <- cluster_sep_spatial * 0.4
        B_k[(ceiling(M/3)+1):M, ] <- -cluster_sep_spatial * 0.15
      } else if (k == 2) {
        # Distal-dominant (calf and shin)
        B_k <- matrix(0, M, r)
        B_k[1:ceiling(2*M/3), ] <- -cluster_sep_spatial * 0.15
        B_k[(ceiling(2*M/3)+1):M, ] <- cluster_sep_spatial * 0.4
      } else {
        # Balanced/mixed pattern with cluster-specific variation
        # Use cluster-specific seed to ensure unique patterns for each cluster
        set.seed(seed + k * 100)
        B_k <- matrix(rnorm(M * r, 0, cluster_sep_spatial * 0.1), M, r)
        
        # Add cluster-specific structured bias to ensure distinctness
        # Create a sinusoidal pattern that varies by cluster
        for (j in seq_len(r)) {
          muscle_idx <- seq_len(M)
          # Sinusoidal modulation with cluster-specific phase
          phase <- 2 * pi * (k - 3) / max(1, K - 2)
          modulation <- sin(2 * pi * muscle_idx / M + phase)
          B_k[, j] <- B_k[, j] + cluster_sep_spatial * 0.15 * modulation
        }
      }
    }
    
    # Generate Lambda_k with cluster-specific bias
    Lambda_k_raw <- Lambda_base + B_k + 
                    matrix(rnorm(M * r, 0, 0.1 * cluster_sep_spatial), M, r)
    
    # Apply non-negativity constraint and normalize
    Lambda_k <- pmax(Lambda_k_raw, 0)
    
    # Ensure at least some positive values
    if (sum(Lambda_k) < 1e-6) {
      Lambda_k <- Lambda_k + matrix(runif(M * r, 0.01, 0.1), M, r)
    }
    
    # Column-wise normalization to make synergies comparable
    for (j in seq_len(r)) {
      col_sum <- sum(Lambda_k[, j])
      if (col_sum > 1e-6) {
        Lambda_k[, j] <- Lambda_k[, j] / col_sum
      }
    }
    
    Lambda_list[[k]] <- Lambda_k
    spatial_params[[k]] <- list(bias_pattern = B_k)
  }
  
  # ========================================================================
  # AXIS 2: TEMPORAL COORDINATION - Phase and Frequency of Neural Drives
  # ========================================================================
  
  temporal_params <- vector("list", K)
  
  # Define coordination modes that will be cycled through for K>3
  mode_cycle <- c("rigid", "alternating", "adaptive")
  
  for (k in seq_len(K)) {
    if (K == 1) {
      # Single cluster: default moderate coordination
      coord_mode <- "moderate"
      phase_pattern <- "moderate"
      # Base frequency at midpoint
      base_freq <- mean(freq_range)
    } else if (K == 2) {
      # Two clusters: rigid vs alternating
      if (k == 1) {
        coord_mode <- "rigid"
        phase_pattern <- "cocontraction"
        base_freq <- freq_range[1]
      } else {
        coord_mode <- "alternating"
        phase_pattern <- "alternating"
        base_freq <- freq_range[2]
      }
    } else {
      # Three or more clusters: cycle through modes and distribute frequencies
      # Cycle coordination modes across clusters
      mode_idx <- ((k - 1) %% 3) + 1
      coord_mode <- mode_cycle[mode_idx]
      
      # Assign phase pattern based on mode
      if (coord_mode == "rigid") {
        phase_pattern <- "cocontraction"
      } else if (coord_mode == "alternating") {
        phase_pattern <- "alternating"
      } else {
        phase_pattern <- "random"
      }
      
      # Distribute base frequencies evenly across freq_range for all K clusters
      # This ensures each cluster has a unique base frequency
      base_freq <- freq_range[1] + (k - 1) * diff(freq_range) / (K - 1)
    }
    
    # Generate cluster-specific frequencies with variation
    freq_variation <- cluster_sep_temporal * 0.3 * diff(freq_range)
    # Use cluster-specific seed for reproducibility
    set.seed(seed + k * 1000)
    freqs_k <- rep(base_freq, r) + rnorm(r, 0, freq_variation)
    freqs_k <- pmax(freqs_k, freq_range[1])
    freqs_k <- pmin(freqs_k, freq_range[2])
    
    # Generate cluster-specific phase offsets
    # Add cluster-specific base phase offset to ensure uniqueness
    cluster_phase_offset <- 2 * pi * (k - 1) / K
    
    if (phase_pattern == "cocontraction") {
      # All synergies in-phase (phi approximately 0) + cluster offset
      set.seed(seed + k * 1000 + 1)
      phases_k <- cluster_phase_offset + rnorm(r, 0, cluster_sep_temporal * 0.1)
    } else if (phase_pattern == "alternating") {
      # Evenly distributed phases over 2*pi + cluster offset
      phases_k <- seq(0, 2*pi, length.out = r + 1)[1:r] + cluster_phase_offset
      # Add small variation based on cluster_sep_temporal
      set.seed(seed + k * 1000 + 1)
      phases_k <- phases_k + rnorm(r, 0, cluster_sep_temporal * 0.2)
    } else {
      # Random phases with cluster-specific seed
      set.seed(seed + k * 1000 + 1)
      phases_k <- runif(r, 0, 2*pi)
    }
    
    # Wrap phases to [0, 2*pi)
    phases_k <- phases_k %% (2 * pi)
    
    temporal_params[[k]] <- list(
      coordination_mode = coord_mode,
      base_frequency = base_freq,
      frequencies = freqs_k,
      phase_pattern = phase_pattern,
      phases = phases_k
    )
  }
  
  # Reset seed to original value for subsequent operations
  set.seed(seed)
  
  # ========================================================================
  # AXIS 3: STABILITY-ADAPTABILITY - Noise and Variability Profiles
  # ========================================================================
  
  stability_params <- vector("list", K)
  
  # Define stability profiles that will be cycled through for K>3
  profile_cycle <- c("stable", "moderate", "adaptive")
  
  # Define parameter ranges for interpolation
  sigma_white_range <- c(0.02, 0.10)
  sigma_low_range <- c(0.01, 0.05)
  filter_length_range <- c(3, 8)
  
  for (k in seq_len(K)) {
    if (K == 1) {
      # Single cluster: moderate control
      control_profile <- "moderate"
      sigma_white <- 0.05
      sigma_low <- 0.02
      filter_length <- 5
    } else if (K == 2) {
      # Two clusters: stable vs adaptive
      if (k == 1) {
        control_profile <- "stable"
        sigma_white <- sigma_white_range[1]
        sigma_low <- sigma_low_range[1]
        filter_length <- filter_length_range[1]
      } else {
        control_profile <- "adaptive"
        sigma_white <- sigma_white_range[2]
        sigma_low <- sigma_low_range[2]
        filter_length <- filter_length_range[2]
      }
    } else {
      # Three or more clusters: cycle through profiles and distribute parameters
      # Cycle stability profiles across clusters
      profile_idx <- ((k - 1) %% 3) + 1
      control_profile <- profile_cycle[profile_idx]
      
      # Distribute noise parameters evenly across ranges for all K clusters
      # This ensures each cluster has unique noise parameters
      sigma_white <- sigma_white_range[1] + (k - 1) * diff(sigma_white_range) / (K - 1)
      sigma_low <- sigma_low_range[1] + (k - 1) * diff(sigma_low_range) / (K - 1)
      filter_length <- filter_length_range[1] + (k - 1) * diff(filter_length_range) / (K - 1)
    }
    
    # Scale noise parameters by cluster_sep_stability
    # Higher cluster_sep_stability increases the spread between clusters
    sigma_white_k <- sigma_white * (1 + cluster_sep_stability * 0.5)
    sigma_low_k <- sigma_low * (1 + cluster_sep_stability * 0.5)
    filter_length_k <- max(3L, min(as.integer(filter_length * (1 + cluster_sep_stability * 0.3)), 
                                    as.integer(T_each / 10)))
    
    stability_params[[k]] <- list(
      control_profile = control_profile,
      sigma_white = sigma_white_k,
      sigma_low = sigma_low_k,
      filter_length = filter_length_k
    )
  }
  
  # ========================================================================
  # GENERATE ACTIVATION SIGNALS C(t) AND EMG SIGNALS X(t)
  # ========================================================================
  
  C_list <- vector("list", N)
  X_list <- vector("list", N)
  
  for (i in seq_len(N)) {
    k_i <- z_true[i]
    Lambda <- Lambda_list[[k_i]]
    
    # Get cluster-specific temporal parameters
    freqs <- temporal_params[[k_i]]$frequencies
    phases <- temporal_params[[k_i]]$phases
    
    # Random amplitudes for this subject
    amps <- runif(r, amp_range[1], amp_range[2])
    
    # Time vector (normalized to [0, 1])
    time_vec <- seq(0, 1, length.out = T_each)
    
    # Generate periodic activation signals with noise
    C <- matrix(0, nrow = T_each, ncol = r)
    for (j in seq_len(r)) {
      # Sine wave + Gaussian noise
      signal <- amps[j] * sin(2 * pi * freqs[j] * time_vec + phases[j])
      signal <- signal + rnorm(T_each, mean = 0, sd = 0.1 * amps[j])
      
      # Rectification (non-negativity) + baseline activation
      C[, j] <- base_level + pmax(0, signal)
    }
    
    C_list[[i]] <- C
    
    # Generate EMG signals X(t) = Lambda * C(t)^T + noise
    # X(t) = Lambda * C(t)^T, then transpose to get (T_each x M)
    X_clean <- C %*% t(Lambda)  # (T_each x r) %*% (r x M) = (T_each x M)
    
    # Get cluster-specific stability parameters
    sigma_white_k <- stability_params[[k_i]]$sigma_white
    sigma_low_k <- stability_params[[k_i]]$sigma_low
    filter_length_k <- stability_params[[k_i]]$filter_length
    
    # Add white noise
    X_noisy <- X_clean + sigma_white_k * matrix(rnorm(T_each * M),
                                                 nrow = T_each, ncol = M)
    
    # Add low-frequency noise (temporal smoothing per channel)
    if (sigma_low_k > 0) {
      # Generate white noise matrix
      W <- matrix(rnorm(T_each * M), nrow = T_each, ncol = M)
      
      # Apply moving average filter per channel to create low-frequency noise
      kernel <- rep(1 / filter_length_k, filter_length_k)
      
      for (j in seq_len(M)) {
        # Apply filter with sides = 1 to avoid NAs at both ends
        lf <- stats::filter(W[, j], kernel, sides = 1)
        # Replace leading NAs with 0
        lf[is.na(lf)] <- 0
        W[, j] <- as.numeric(lf)
      }
      
      X_noisy <- X_noisy + sigma_low_k * W
    }
    
    # Apply rectification if requested (half-wave rectification)
    if (rectify_output) {
      X_noisy[X_noisy < 0] <- 0
    }
    
    X_list[[i]] <- X_noisy
  }
  
  # Return all components including parameter specifications
  list(
    X_list = X_list,
    C_list = C_list,
    Lambda_list = Lambda_list,
    z_true = z_true,
    spatial_params = spatial_params,
    temporal_params = temporal_params,
    stability_params = stability_params
  )
}
