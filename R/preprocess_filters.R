# EMG Preprocessing - Filter Functions
# Part of the synergyMixR package

#' Remove DC Offset from EMG Signal
#'
#' Removes the DC component (baseline offset) from EMG signals using
#' various methods.
#'
#' @param data EMG signal. Can be a numeric vector (single channel),
#'   a matrix (T x M, time x channels), or a list of matrices.
#' @param method Method for DC removal:
#'   \itemize{
#'     \item \code{"mean"}: Subtract the mean (default)
#'     \item \code{"median"}: Subtract the median (robust to outliers
#'     \item \code{"linear"}: Remove linear trend (detrend)
#'     \item \code{"polynomial"}: Remove polynomial trend
#'   }
#' @param poly_order Polynomial order when \code{method = "polynomial"}.
#'   Default is 1 (equivalent to linear detrending).
#'
#' @return Data with DC offset removed, same structure as input.
#'
#' @examples
#' # Single channel
#' signal <- rnorm(1000) + 0.5
#' clean <- remove_dc_offset(signal)
#'
#' # Multi-channel matrix
#' emg <- matrix(rnorm(1000 * 8), nrow = 1000, ncol = 8)
#' emg_clean <- remove_dc_offset(emg, method = "linear")
#'
#' @export
remove_dc_offset <- function(data, method = "mean", poly_order = 1) {
  # Validate method

method <- match.arg(method, c("mean", "median", "linear", "polynomial"))

  # Handle different input types
  if (is.list(data) && !is.data.frame(data)) {
    return(lapply(data, remove_dc_offset, method = method, poly_order = poly_order))
  }

  if (is.vector(data) && !is.matrix(data)) {
    return(remove_dc_offset_vector(data, method, poly_order))
  }

  if (is.matrix(data) || is.data.frame(data)) {
    data <- as.matrix(data)
    result <- apply(data, 2, remove_dc_offset_vector, method = method, poly_order = poly_order)
    return(result)
  }

  stop("Input must be a vector, matrix, or list of matrices")
}

#' Remove DC offset from a single vector
#' @keywords internal
remove_dc_offset_vector <- function(x, method, poly_order) {
  n <- length(x)
  t_idx <- seq_len(n)

  switch(method,
    "mean" = x - mean(x, na.rm = TRUE),
    "median" = x - stats::median(x, na.rm = TRUE),
    "linear" = {
      fit <- stats::lm(x ~ t_idx)
      x - stats::fitted(fit)
    },
    "polynomial" = {
      fit <- stats::lm(x ~ stats::poly(t_idx, degree = poly_order, raw = TRUE))
      x - stats::fitted(fit)
    }
  )
}


#' Bandpass Filter for EMG Signals
#'
#' Applies a bandpass filter to EMG signals to extract the frequency band
#' of interest (typically 20-450 Hz for surface EMG).
#'
#' @param data EMG signal. Can be a numeric vector (single channel),
#'   a matrix (T x M, time x channels), or a list of matrices.
#' @param sampling_rate Sampling frequency in Hz.
#' @param low_cutoff Lower cutoff frequency in Hz. Default is 20 Hz.
#' @param high_cutoff Upper cutoff frequency in Hz. Default is 450 Hz.
#' @param filter_order Filter order. Default is 4 (creates 4th order filter,
#'   8th order effective with zero-phase filtering).
#' @param filter_type Type of filter: \code{"butterworth"} (default),
#'   \code{"chebyshev1"}, or \code{"chebyshev2"}.
#' @param zero_phase Logical. If TRUE (default), applies zero-phase filtering
#'   (forward-backward filtering) to avoid phase distortion.
#'
#' @return Filtered EMG signal, same structure as input.
#'
#' @details
#' The bandpass filter removes:
#' \itemize{
#'   \item Low-frequency noise and motion artifacts (below \code{low_cutoff})
#'   \item High-frequency noise (above \code{high_cutoff})
#' }
#'
#' Standard EMG bandwidth is 20-450 Hz for surface EMG. The lower cutoff
#' helps remove motion artifacts, while the upper cutoff prevents aliasing
#' and removes high-frequency noise.
#'
#' @examples
#' \dontrun{
#' # Generate sample EMG-like signal
#' fs <- 1000  # 1000 Hz sampling rate
#' t <- seq(0, 1, by = 1/fs)
#' emg <- rnorm(length(t)) + 0.1 * sin(2 * pi * 10 * t)  # with low-freq noise
#'
#' # Apply bandpass filter
#' filtered <- bandpass_filter(emg, sampling_rate = fs)
#' }
#'
#' @export
bandpass_filter <- function(
  data,
  sampling_rate,
  low_cutoff = 20,
  high_cutoff = 450,
  filter_order = 4,
  filter_type = "butterworth",
  zero_phase = TRUE
) {
  # Validate inputs
  if (missing(sampling_rate)) {
    stop("sampling_rate is required")
  }

  nyquist <- sampling_rate / 2

  if (low_cutoff <= 0 || low_cutoff >= nyquist) {
    stop("low_cutoff must be between 0 and Nyquist frequency (", nyquist, " Hz)")
  }
  if (high_cutoff <= low_cutoff || high_cutoff >= nyquist) {
    stop("high_cutoff must be between low_cutoff and Nyquist frequency (", nyquist, " Hz)")
  }

  filter_type <- match.arg(filter_type, c("butterworth", "chebyshev1", "chebyshev2"))

  # Design filter
  W <- c(low_cutoff, high_cutoff) / nyquist

  bf <- switch(filter_type,
    "butterworth" = signal::butter(filter_order, W, type = "pass"),
    "chebyshev1" = signal::cheby1(filter_order, 0.5, W, type = "pass"),
    "chebyshev2" = signal::cheby2(filter_order, 20, W, type = "pass")
  )

  # Apply filter
  apply_filter(data, bf, zero_phase)
}


#' Notch Filter for Power Line Noise Removal
#'
#' Applies a notch (band-stop) filter to remove power line interference
#' (50 Hz or 60 Hz) and optionally its harmonics.
#'
#' @param data EMG signal. Can be a numeric vector (single channel),
#'   a matrix (T x M, time x channels), or a list of matrices.
#' @param sampling_rate Sampling frequency in Hz.
#' @param notch_freq Frequency to remove in Hz. Default is 50 Hz.
#'   Use 60 Hz for regions with 60 Hz power lines (e.g., Americas, Japan).
#' @param Q Quality factor determining the filter bandwidth.
#'   Higher Q = narrower notch. Default is 30.
#' @param harmonics Logical. If TRUE (default), also removes harmonics
#'   (2x, 3x, etc. of notch_freq up to Nyquist frequency).
#' @param zero_phase Logical. If TRUE (default), applies zero-phase filtering.
#'
#' @return Filtered EMG signal with power line noise removed.
#'
#' @details
#' Power line interference (50/60 Hz) is a common artifact in EMG recordings.
#' The notch filter creates a narrow band-stop at the specified frequency.
#' When \code{harmonics = TRUE}, harmonics (100/120 Hz, 150/180 Hz, etc.)
#' are also removed.
#'
#' @examples
#' \dontrun{
#' # Remove 50 Hz power line noise
#' clean_emg <- notch_filter(emg, sampling_rate = 1000, notch_freq = 50)
#'
#' # Remove 60 Hz (Americas) with harmonics
#' clean_emg <- notch_filter(emg, sampling_rate = 1000, notch_freq = 60, harmonics = TRUE)
#' }
#'
#' @export
notch_filter <- function(
  data,
  sampling_rate,
  notch_freq = 50,
  Q = 30,
  harmonics = TRUE,
  zero_phase = TRUE
) {
  # Validate inputs
  if (missing(sampling_rate)) {
    stop("sampling_rate is required")
  }

  nyquist <- sampling_rate / 2

  if (notch_freq <= 0 || notch_freq >= nyquist) {
    stop("notch_freq must be between 0 and Nyquist frequency (", nyquist, " Hz)")
  }

  # Determine frequencies to filter
  if (harmonics) {
    freqs <- notch_freq * seq_len(floor(nyquist / notch_freq))
    freqs <- freqs[freqs < nyquist * 0.95]  # Leave margin from Nyquist
  } else {
    freqs <- notch_freq
  }

  result <- data

  # Apply notch filter for each frequency
  for (freq in freqs) {
    bandwidth <- freq / Q
    low <- (freq - bandwidth / 2) / nyquist
    high <- (freq + bandwidth / 2) / nyquist

    # Ensure valid range
    if (low > 0 && high < 1) {
      bf <- signal::butter(2, c(low, high), type = "stop")
      result <- apply_filter(result, bf, zero_phase)
    }
  }

  result
}


#' Apply Filter to EMG Data
#'
#' Internal function to apply a filter to various data structures.
#'
#' @param data EMG data (vector, matrix, or list)
#' @param bf Filter object from signal package
#' @param zero_phase Whether to use zero-phase filtering
#'
#' @return Filtered data
#' @keywords internal
apply_filter <- function(data, bf, zero_phase) {
  filter_func <- if (zero_phase) signal::filtfilt else signal::filter

  if (is.list(data) && !is.data.frame(data)) {
    return(lapply(data, apply_filter, bf = bf, zero_phase = zero_phase))
  }

  if (is.vector(data) && !is.matrix(data)) {
    return(as.numeric(filter_func(bf, data)))
  }

  if (is.matrix(data) || is.data.frame(data)) {
    data <- as.matrix(data)
    result <- apply(data, 2, function(x) as.numeric(filter_func(bf, x)))
    return(result)
  }

  stop("Input must be a vector, matrix, or list of matrices")
}


#' Rectify EMG Signal
#'
#' Applies rectification to EMG signals, converting all values to positive.
#'
#' @param data EMG signal. Can be a numeric vector (single channel),
#'   a matrix (T x M, time x channels), or a list of matrices.
#' @param method Rectification method:
#'   \itemize{
#'     \item \code{"full_wave"}: Take absolute value (default)
#'     \item \code{"half_wave"}: Set negative values to zero
#'   }
#'
#' @return Rectified EMG signal.
#'
#' @details
#' Rectification is a standard step in EMG processing that prepares the
#' signal for envelope extraction. Full-wave rectification (absolute value)
#' is most commonly used.
#'
#' @examples
#' signal <- rnorm(1000)
#' rectified <- rectify(signal)  # full-wave rectification
#' half_rect <- rectify(signal, method = "half_wave")
#'
#' @export
rectify <- function(data, method = "full_wave") {
  method <- match.arg(method, c("full_wave", "half_wave"))

  rectify_func <- switch(method,
    "full_wave" = abs,
    "half_wave" = function(x) pmax(x, 0)
  )

  if (is.list(data) && !is.data.frame(data)) {
    return(lapply(data, rectify, method = method))
  }

  if (is.vector(data) && !is.matrix(data)) {
    return(rectify_func(data))
  }

  if (is.matrix(data) || is.data.frame(data)) {
    return(rectify_func(as.matrix(data)))
  }

  stop("Input must be a vector, matrix, or list of matrices")
}
