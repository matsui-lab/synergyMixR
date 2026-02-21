# EMG Preprocessing - Normalization Functions
# Part of the synergyMixR package

#' Normalize EMG Signal
#'
#' Normalizes EMG signals using various methods to allow comparison
#' across channels, subjects, or sessions.
#'
#' @param data EMG signal. Can be a numeric vector (single channel),
#'   a matrix (T x M, time x channels), or a list of matrices.
#' @param method Normalization method:
#'   \itemize{
#'     \item \code{"max"}: Divide by maximum value (default)
#'     \item \code{"mvc"}: Divide by Maximum Voluntary Contraction reference
#'     \item \code{"mean"}: Divide by mean value
#'     \item \code{"zscore"}: Z-score standardization (subtract mean, divide by SD)
#'     \item \code{"range"}: Min-max normalization to [0, 1]
#'     \item \code{"median"}: Divide by median value (robust to outliers)
#'   }
#' @param reference Reference values for MVC normalization.
#'   Can be a single value (applied to all channels) or a vector
#'   (one value per channel). Only used when \code{method = "mvc"}.
#' @param per_channel Logical. If TRUE (default), normalize each channel
#'   independently. If FALSE, use global statistics across all channels.
#'
#' @return Normalized EMG signal, same structure as input.
#'
#' @details
#' Normalization is essential for:
#' \itemize{
#'   \item Comparing EMG across different muscles
#'   \item Comparing EMG across different subjects
#'   \item Comparing EMG across different recording sessions
#' }
#'
#' \strong{MVC normalization} (\code{method = "mvc"}) is the gold standard
#' in EMG research, expressing activation as a percentage of maximum
#' voluntary contraction. Requires separate MVC trials for reference.
#'
#' \strong{Max normalization} (\code{method = "max"}) is useful when MVC
#' data is not available, normalizing to the peak value in the data.
#'
#' @examples
#' # Max normalization (0 to 1 range)
#' emg <- matrix(rnorm(1000), ncol = 4)
#' norm_emg <- normalize_emg(emg, method = "max")
#'
#' # MVC normalization (% MVC)
#' mvc_values <- c(100, 150, 80, 120)  # MVC for each channel
#' norm_emg <- normalize_emg(emg, method = "mvc", reference = mvc_values)
#'
#' # Z-score normalization
#' norm_emg <- normalize_emg(emg, method = "zscore")
#'
#' @export
normalize_emg <- function(
  data,
  method = "max",
  reference = NULL,
  per_channel = TRUE
) {
  method <- match.arg(method, c("max", "mvc", "mean", "zscore", "range", "median"))

  # Validate MVC reference
  if (method == "mvc" && is.null(reference)) {
    stop("reference values required for MVC normalization")
  }

  # Handle list input
  if (is.list(data) && !is.data.frame(data)) {
    return(lapply(data, normalize_emg,
      method = method,
      reference = reference,
      per_channel = per_channel
    ))
  }

  # Handle vector input
  if (is.vector(data) && !is.matrix(data)) {
    ref <- if (!is.null(reference)) reference[1] else NULL
    return(normalize_vector(data, method, ref))
  }

  # Handle matrix input
  if (is.matrix(data) || is.data.frame(data)) {
    data <- as.matrix(data)
    return(normalize_matrix(data, method, reference, per_channel))
  }

  stop("Input must be a vector, matrix, or list of matrices")
}


#' Normalize a single vector
#' @keywords internal
normalize_vector <- function(x, method, reference = NULL) {
  switch(method,
    "max" = x / max(abs(x), na.rm = TRUE),
    "mvc" = {
      if (is.null(reference) || reference == 0) {
        warning("Invalid MVC reference, using max normalization")
        x / max(abs(x), na.rm = TRUE)
      } else {
        x / reference * 100  # Express as percentage
      }
    },
    "mean" = x / mean(abs(x), na.rm = TRUE),
    "zscore" = (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE),
    "range" = {
      min_val <- min(x, na.rm = TRUE)
      max_val <- max(x, na.rm = TRUE)
      if (max_val == min_val) {
        rep(0, length(x))
      } else {
        (x - min_val) / (max_val - min_val)
      }
    },
    "median" = x / stats::median(abs(x), na.rm = TRUE)
  )
}


#' Normalize a matrix
#' @keywords internal
normalize_matrix <- function(data, method, reference, per_channel) {
  n_channels <- ncol(data)

  if (per_channel) {
    # Normalize each channel independently
    result <- matrix(NA, nrow = nrow(data), ncol = ncol(data))

    for (j in seq_len(n_channels)) {
      ref <- if (!is.null(reference)) {
        if (length(reference) == 1) reference else reference[j]
      } else {
        NULL
      }
      result[, j] <- normalize_vector(data[, j], method, ref)
    }

    # Preserve column names
    colnames(result) <- colnames(data)
    return(result)

  } else {
    # Use global statistics
    switch(method,
      "max" = data / max(abs(data), na.rm = TRUE),
      "mvc" = {
        if (is.null(reference) || all(reference == 0)) {
          warning("Invalid MVC reference, using max normalization")
          data / max(abs(data), na.rm = TRUE)
        } else if (length(reference) == 1) {
          data / reference * 100
        } else {
          # Apply channel-specific reference even in global mode
          t(t(data) / reference) * 100
        }
      },
      "mean" = data / mean(abs(data), na.rm = TRUE),
      "zscore" = {
        global_mean <- mean(data, na.rm = TRUE)
        global_sd <- stats::sd(as.vector(data), na.rm = TRUE)
        (data - global_mean) / global_sd
      },
      "range" = {
        min_val <- min(data, na.rm = TRUE)
        max_val <- max(data, na.rm = TRUE)
        if (max_val == min_val) {
          matrix(0, nrow = nrow(data), ncol = ncol(data))
        } else {
          (data - min_val) / (max_val - min_val)
        }
      },
      "median" = data / stats::median(abs(data), na.rm = TRUE)
    )
  }
}


#' Compute MVC Reference Values
#'
#' Computes Maximum Voluntary Contraction (MVC) reference values from
#' MVC trial data for use in EMG normalization.
#'
#' @param mvc_data MVC trial data. Can be a matrix (T x M) or a list
#'   of matrices (multiple MVC trials).
#' @param method Method to compute MVC:
#'   \itemize{
#'     \item \code{"peak"}: Maximum value (default)
#'     \item \code{"mean_peak"}: Mean of top N% values
#'     \item \code{"rms_peak"}: RMS of the peak window
#'   }
#' @param window_ms Window size in ms for peak detection (used with
#'   \code{"rms_peak"}). Default is 500 ms.
#' @param top_percent Percentage of top values to average (used with
#'   \code{"mean_peak"}). Default is 5.
#' @param sampling_rate Sampling rate in Hz (required for \code{"rms_peak"}).
#'
#' @return Numeric vector of MVC values, one per channel.
#'
#' @examples
#' # From single MVC trial
#' mvc_trial <- matrix(abs(rnorm(1000 * 4)), ncol = 4)
#' mvc_ref <- compute_mvc_reference(mvc_trial)
#'
#' # From multiple MVC trials (takes maximum across trials)
#' mvc_trials <- list(
#'   trial1 = matrix(abs(rnorm(1000 * 4)), ncol = 4),
#'   trial2 = matrix(abs(rnorm(1000 * 4)), ncol = 4)
#' )
#' mvc_ref <- compute_mvc_reference(mvc_trials)
#'
#' @export
compute_mvc_reference <- function(
  mvc_data,
  method = "peak",
  window_ms = 500,
  top_percent = 5,
  sampling_rate = NULL
) {
  method <- match.arg(method, c("peak", "mean_peak", "rms_peak"))

  if (method == "rms_peak" && is.null(sampling_rate)) {
    stop("sampling_rate required for rms_peak method")
  }

  # Handle list of MVC trials
  if (is.list(mvc_data) && !is.data.frame(mvc_data)) {
    mvc_values_list <- lapply(mvc_data, compute_mvc_reference,
      method = method,
      window_ms = window_ms,
      top_percent = top_percent,
      sampling_rate = sampling_rate
    )
    # Take maximum across trials
    mvc_matrix <- do.call(rbind, mvc_values_list)
    return(apply(mvc_matrix, 2, max))
  }

  # Handle single matrix
  mvc_data <- as.matrix(mvc_data)
  n_channels <- ncol(mvc_data)

  result <- numeric(n_channels)

  for (j in seq_len(n_channels)) {
    x <- mvc_data[, j]

    result[j] <- switch(method,
      "peak" = max(abs(x), na.rm = TRUE),
      "mean_peak" = {
        n_top <- max(1, floor(length(x) * top_percent / 100))
        sorted <- sort(abs(x), decreasing = TRUE)
        mean(sorted[seq_len(n_top)])
      },
      "rms_peak" = {
        window_samples <- round(window_ms * sampling_rate / 1000)
        if (window_samples > length(x)) window_samples <- length(x)

        # Find the window with maximum RMS
        max_rms <- 0
        for (i in seq_len(length(x) - window_samples + 1)) {
          window_rms <- sqrt(mean(x[i:(i + window_samples - 1)]^2))
          if (window_rms > max_rms) max_rms <- window_rms
        }
        max_rms
      }
    )
  }

  result
}
