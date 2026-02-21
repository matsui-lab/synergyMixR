# EMG Preprocessing - Envelope Extraction Functions
# Part of the synergyMixR package

#' Extract Envelope from EMG Signal
#'
#' Extracts the amplitude envelope from rectified EMG signals using
#' various methods.
#'
#' @param data Rectified EMG signal. Can be a numeric vector (single channel),
#'   a matrix (T x M, time x channels), or a list of matrices.
#'   Should typically be rectified before envelope extraction.
#' @param sampling_rate Sampling frequency in Hz.
#' @param method Envelope extraction method:
#'   \itemize{
#'     \item \code{"lowpass"}: Low-pass filter (default, most common)
#'     \item \code{"rms"}: Root Mean Square with sliding window
#'     \item \code{"moving_average"}: Simple moving average
#'   }
#' @param cutoff Cutoff frequency for low-pass filter in Hz.
#'   Only used when \code{method = "lowpass"}. Default is 6 Hz.
#' @param window_ms Window size in milliseconds for RMS or moving average.
#'   Only used when \code{method} is \code{"rms"} or \code{"moving_average"}.
#'   Default is 50 ms.
#' @param filter_order Filter order for low-pass method. Default is 4.
#' @param zero_phase Logical. If TRUE (default), applies zero-phase filtering
#'   for the low-pass method.
#'
#' @return EMG envelope signal, same structure as input.
#'
#' @details
#' The envelope represents the amplitude modulation of the EMG signal,
#' showing muscle activation patterns over time. Common approaches:
#'
#' \itemize{
#'   \item \strong{Low-pass filter}: Applies a Butterworth low-pass filter
#'     to the rectified signal. Cutoff frequency typically 4-10 Hz.
#'   \item \strong{RMS}: Computes root mean square in a sliding window.
#'     Good for preserving amplitude information.
#'   \item \strong{Moving average}: Simple smoothing with a sliding window.
#'     Fast but may introduce artifacts at window edges.
#' }
#'
#' @examples
#' \dontrun{
#' # Generate and process EMG
#' fs <- 1000
#' raw_emg <- rnorm(1000)
#' rect_emg <- rectify(raw_emg)
#'
#' # Low-pass envelope (most common)
#' env_lp <- extract_envelope(rect_emg, sampling_rate = fs, method = "lowpass")
#'
#' # RMS envelope
#' env_rms <- extract_envelope(rect_emg, sampling_rate = fs,
#'                             method = "rms", window_ms = 100)
#' }
#'
#' @export
extract_envelope <- function(
  data,
  sampling_rate,
  method = "lowpass",
  cutoff = 6,
  window_ms = 50,
  filter_order = 4,
  zero_phase = TRUE
) {
  # Validate inputs
  if (missing(sampling_rate)) {
    stop("sampling_rate is required")
  }

  method <- match.arg(method, c("lowpass", "rms", "moving_average"))

  # Handle different input types
  if (is.list(data) && !is.data.frame(data)) {
    return(lapply(data, extract_envelope,
      sampling_rate = sampling_rate,
      method = method,
      cutoff = cutoff,
      window_ms = window_ms,
      filter_order = filter_order,
      zero_phase = zero_phase
    ))
  }

  if (is.vector(data) && !is.matrix(data)) {
    return(extract_envelope_vector(data, sampling_rate, method,
      cutoff, window_ms, filter_order, zero_phase))
  }

  if (is.matrix(data) || is.data.frame(data)) {
    data <- as.matrix(data)
    result <- apply(data, 2, extract_envelope_vector,
      sampling_rate = sampling_rate,
      method = method,
      cutoff = cutoff,
      window_ms = window_ms,
      filter_order = filter_order,
      zero_phase = zero_phase
    )
    return(result)
  }

  stop("Input must be a vector, matrix, or list of matrices")
}


#' Extract envelope from a single vector
#' @keywords internal
extract_envelope_vector <- function(
  x, sampling_rate, method, cutoff,
  window_ms, filter_order, zero_phase
) {
  switch(method,
    "lowpass" = envelope_lowpass(x, sampling_rate, cutoff, filter_order, zero_phase),
    "rms" = envelope_rms(x, sampling_rate, window_ms),
    "moving_average" = envelope_moving_average(x, sampling_rate, window_ms)
  )
}


#' Low-pass filter envelope extraction
#' @keywords internal
envelope_lowpass <- function(x, sampling_rate, cutoff, filter_order, zero_phase) {
  nyquist <- sampling_rate / 2

  if (cutoff >= nyquist) {
    warning("cutoff frequency adjusted to 95% of Nyquist frequency")
    cutoff <- nyquist * 0.95
  }

  W <- cutoff / nyquist
  bf <- signal::butter(filter_order, W, type = "low")

  filter_func <- if (zero_phase) signal::filtfilt else signal::filter
  as.numeric(filter_func(bf, x))
}


#' RMS envelope extraction
#' @keywords internal
envelope_rms <- function(x, sampling_rate, window_ms) {
  window_samples <- round(window_ms * sampling_rate / 1000)
  if (window_samples < 1) window_samples <- 1
  if (window_samples > length(x)) window_samples <- length(x)

  n <- length(x)
  result <- numeric(n)

  # Half window for centering
  half_win <- floor(window_samples / 2)

  for (i in seq_len(n)) {
    start_idx <- max(1, i - half_win)
    end_idx <- min(n, i + half_win)
    result[i] <- sqrt(mean(x[start_idx:end_idx]^2))
  }

  result
}


#' Moving average envelope extraction
#' @keywords internal
envelope_moving_average <- function(x, sampling_rate, window_ms) {
  window_samples <- round(window_ms * sampling_rate / 1000)
  if (window_samples < 1) window_samples <- 1
  if (window_samples > length(x)) window_samples <- length(x)

  # Use stats::filter for moving average
  weights <- rep(1 / window_samples, window_samples)
  result <- stats::filter(x, weights, sides = 2)

  # Handle NA values at edges (use partial window)
  na_idx <- which(is.na(result))
  for (i in na_idx) {
    start_idx <- max(1, i - floor(window_samples / 2))
    end_idx <- min(length(x), i + floor(window_samples / 2))
    result[i] <- mean(x[start_idx:end_idx])
  }

  as.numeric(result)
}
