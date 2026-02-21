# EMG Preprocessing - Main Pipeline Functions
# Part of the synergyMixR package

#' EMG Preprocessing Pipeline
#'
#' Applies a complete preprocessing pipeline to raw EMG data, transforming
#' it into the format required for MFA/MPCA analysis.
#'
#' @param raw_data Raw EMG data. Can be:
#'   \itemize{
#'     \item A matrix (T x M) for a single trial
#'     \item A list of matrices for multiple trials/subjects
#'   }
#' @param sampling_rate Sampling frequency in Hz.
#' @param pipeline_config Preprocessing configuration. Can be:
#'   \itemize{
#'     \item A character string: preset name ("default", "gait", "isometric", "fine_motor")
#'     \item A list: custom configuration (see Details)
#'   }
#' @param mvc_data Optional MVC reference data for MVC normalization.
#'   Can be a matrix (T x M) or a list of matrices.
#' @param channel_names Optional character vector of channel/muscle names.
#' @param verbose Logical. If TRUE (default), prints progress information.
#'
#' @return An object of class \code{"preprocessed_emg"} containing:
#'   \itemize{
#'     \item \code{data}: Preprocessed data in list_of_data format (ready for MFA/MPCA)
#'     \item \code{processing_info}: Record of all preprocessing steps and parameters
#'     \item \code{channel_names}: Channel names
#'     \item \code{sampling_rate}: Original sampling rate
#'   }
#'
#' @details
#' \subsection{Processing Flow}{
#' The preprocessing pipeline transforms raw EMG signals into activation
#' patterns suitable for muscle synergy analysis:
#'
#' \preformatted{
#' Raw EMG --> [DC Removal] --> [Bandpass] --> [Notch] --> [Rectify]
#'                 |               |             |            |
#'            Remove bias    Extract signal  Remove power  Convert to
#'            and drift      bandwidth       line noise    positive
#'                                                             |
#'                                                             v
#'         Preprocessed <-- [Normalize] <-- [Envelope] <-------+
#'         (list_of_data)        |              |
#'                          Scale for      Extract amplitude
#'                          comparison     modulation
#' }
#' }
#'
#' \subsection{Processing Steps}{
#' \describe{
#'   \item{\strong{1. DC Offset Removal}}{
#'     Removes electrode DC bias and signal drift.
#'     Methods: "mean" (subtract mean), "linear" (detrend), "polynomial".
#'     Typical use: Always recommended as first step.
#'   }
#'   \item{\strong{2. Bandpass Filter}}{
#'     Extracts the EMG frequency band (typically 20-450 Hz).
#'     Removes low-frequency motion artifacts (<20 Hz) and
#'     high-frequency noise (>450 Hz).
#'     Parameters: low_cutoff=20, high_cutoff=450, order=4 (Butterworth).
#'   }
#'   \item{\strong{3. Notch Filter}}{
#'     Removes power line interference (50 Hz in Japan/Europe, 60 Hz in Americas).
#'     Can also remove harmonics (100, 150, 200 Hz, etc.).
#'     Parameters: freq=50 or 60, Q=30, harmonics=TRUE.
#'   }
#'   \item{\strong{4. Rectification}}{
#'     Converts bipolar EMG to unipolar by taking absolute value (full-wave)
#'     or setting negatives to zero (half-wave).
#'     Full-wave rectification is standard for envelope extraction.
#'   }
#'   \item{\strong{5. Envelope Extraction}}{
#'     Extracts the amplitude envelope representing muscle activation pattern.
#'     Methods: "lowpass" (4-10 Hz cutoff), "rms" (root mean square),
#'     "moving_average". Lowpass with 6 Hz cutoff is most common.
#'   }
#'   \item{\strong{6. Normalization}}{
#'     Scales signals for comparison across muscles/subjects/sessions.
#'     Methods: "mvc" (percentage of max voluntary contraction - gold standard),
#'     "max" (peak value), "zscore" (standardization), "range" (0-1 scaling).
#'   }
#' }
#' }
#'
#' \subsection{Preset Configurations}{
#' \describe{
#'   \item{\strong{default}}{
#'     General-purpose settings. Bandpass 20-450 Hz, 50 Hz notch,
#'     lowpass envelope (6 Hz), max normalization.
#'   }
#'   \item{\strong{gait}}{
#'     Optimized for walking/gait analysis. Linear detrend,
#'     slower envelope (4 Hz cutoff), MVC normalization.
#'   }
#'   \item{\strong{isometric}}{
#'     For sustained contractions. RMS envelope (100 ms window),
#'     MVC normalization. No time normalization needed.
#'   }
#'   \item{\strong{fine_motor}}{
#'     For precise movements. Wider bandwidth (10-500 Hz),
#'     faster envelope (10 Hz cutoff).
#'   }
#' }
#' }
#'
#' \subsection{Custom Configuration}{
#' \preformatted{
#' config <- list(
#'   dc_removal = list(method = "mean"),
#'   bandpass = list(low = 20, high = 450, order = 4),
#'   notch = list(freq = 50, harmonics = TRUE),
#'   rectify = list(method = "full_wave"),
#'   envelope = list(method = "lowpass", cutoff = 6),
#'   normalize = list(method = "max")
#' )
#' }
#' Set any step to NULL or FALSE to skip it.
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with default settings
#' raw_emg <- matrix(rnorm(1000 * 8), ncol = 8)
#' processed <- preprocess_emg(raw_emg, sampling_rate = 1000)
#'
#' # Use with MFA analysis
#' fit <- mfa_em_fit(processed$data, K = 2, r = 3)
#'
#' # Custom configuration
#' config <- list(
#'   bandpass = list(low = 30, high = 400),
#'   notch = list(freq = 60),  # 60 Hz for Americas
#'   envelope = list(method = "rms", window_ms = 50),
#'   normalize = list(method = "mvc")
#' )
#' processed <- preprocess_emg(raw_emg, sampling_rate = 1000,
#'                             pipeline_config = config,
#'                             mvc_data = mvc_trials)
#' }
#'
#' @seealso
#' \code{\link{emg_preset_config}} for preset configurations,
#' \code{\link{bandpass_filter}}, \code{\link{notch_filter}},
#' \code{\link{rectify}}, \code{\link{extract_envelope}},
#' \code{\link{normalize_emg}} for individual processing functions.
#'
#' @export
preprocess_emg <- function(
  raw_data,
  sampling_rate,
  pipeline_config = "default",
  mvc_data = NULL,
  channel_names = NULL,
  verbose = TRUE
) {
  # Validate inputs
  if (missing(sampling_rate)) {
    stop("sampling_rate is required")
  }

  # Get configuration
  if (is.character(pipeline_config)) {
    config <- emg_preset_config(pipeline_config)
    if (is.null(config)) {
      stop("Unknown preset: ", pipeline_config, ". Use 'default', 'gait', 'isometric', or 'fine_motor'")
    }
  } else if (is.list(pipeline_config)) {
    config <- pipeline_config
  } else {
    stop("pipeline_config must be a preset name or a list")
  }

  # Convert single matrix to list format
  is_single_trial <- is.matrix(raw_data) || is.data.frame(raw_data)
  if (is_single_trial) {
    raw_data <- list(raw_data)
  }

  # Ensure all elements are matrices
  raw_data <- lapply(raw_data, as.matrix)

  # Get channel names
  if (is.null(channel_names)) {
    n_channels <- ncol(raw_data[[1]])
    channel_names <- if (!is.null(colnames(raw_data[[1]]))) {
      colnames(raw_data[[1]])
    } else {
      paste0("Ch", seq_len(n_channels))
    }
  }

  # Compute MVC reference if needed
  mvc_ref <- NULL
  if (!is.null(mvc_data) && !is.null(config$normalize) && config$normalize$method == "mvc") {
    if (verbose) message("Computing MVC reference values...")
    mvc_ref <- compute_mvc_reference(mvc_data, sampling_rate = sampling_rate)
  }

  # Record processing steps
  steps_applied <- character()
  parameters <- list()

  # Apply preprocessing pipeline
  data <- raw_data

  # Step 1: DC offset removal
  if (!is.null(config$dc_removal) && !isFALSE(config$dc_removal)) {
    if (verbose) message("Removing DC offset...")
    dc_params <- config$dc_removal
    data <- lapply(data, remove_dc_offset,
      method = dc_params$method %||% "mean",
      poly_order = dc_params$poly_order %||% 1
    )
    steps_applied <- c(steps_applied, "dc_removal")
    parameters$dc_removal <- dc_params
  }

  # Step 2: Bandpass filter
  if (!is.null(config$bandpass) && !isFALSE(config$bandpass)) {
    if (verbose) message("Applying bandpass filter...")
    bp_params <- config$bandpass
    data <- lapply(data, bandpass_filter,
      sampling_rate = sampling_rate,
      low_cutoff = bp_params$low %||% 20,
      high_cutoff = bp_params$high %||% 450,
      filter_order = bp_params$order %||% 4,
      zero_phase = bp_params$zero_phase %||% TRUE
    )
    steps_applied <- c(steps_applied, "bandpass")
    parameters$bandpass <- bp_params
  }

  # Step 3: Notch filter
  if (!is.null(config$notch) && !isFALSE(config$notch)) {
    if (verbose) message("Applying notch filter...")
    notch_params <- config$notch
    data <- lapply(data, notch_filter,
      sampling_rate = sampling_rate,
      notch_freq = notch_params$freq %||% 50,
      Q = notch_params$Q %||% 30,
      harmonics = notch_params$harmonics %||% TRUE,
      zero_phase = notch_params$zero_phase %||% TRUE
    )
    steps_applied <- c(steps_applied, "notch")
    parameters$notch <- notch_params
  }

  # Step 4: Rectification
  if (!is.null(config$rectify) && !isFALSE(config$rectify)) {
    if (verbose) message("Rectifying signal...")
    rect_params <- config$rectify
    data <- lapply(data, rectify,
      method = rect_params$method %||% "full_wave"
    )
    steps_applied <- c(steps_applied, "rectify")
    parameters$rectify <- rect_params
  }

  # Step 5: Envelope extraction
  if (!is.null(config$envelope) && !isFALSE(config$envelope)) {
    if (verbose) message("Extracting envelope...")
    env_params <- config$envelope
    data <- lapply(data, extract_envelope,
      sampling_rate = sampling_rate,
      method = env_params$method %||% "lowpass",
      cutoff = env_params$cutoff %||% 6,
      window_ms = env_params$window_ms %||% 50,
      zero_phase = env_params$zero_phase %||% TRUE
    )
    steps_applied <- c(steps_applied, "envelope")
    parameters$envelope <- env_params
  }

  # Step 6: Normalization
  if (!is.null(config$normalize) && !isFALSE(config$normalize)) {
    if (verbose) message("Normalizing signal...")
    norm_params <- config$normalize
    data <- lapply(data, normalize_emg,
      method = norm_params$method %||% "max",
      reference = mvc_ref,
      per_channel = norm_params$per_channel %||% TRUE
    )
    steps_applied <- c(steps_applied, "normalize")
    parameters$normalize <- norm_params
    if (!is.null(mvc_ref)) {
      parameters$normalize$mvc_reference <- mvc_ref
    }
  }

  # Set channel names on all matrices
  data <- lapply(data, function(x) {
    if (ncol(x) == length(channel_names)) {
      colnames(x) <- channel_names
    }
    x
  })

  # Unwrap if single trial
  if (is_single_trial) {
    data <- data[[1]]
  }

  if (verbose) message("Preprocessing complete!")

  # Create result object
  result <- structure(
    list(
      data = data,
      processing_info = list(
        steps_applied = steps_applied,
        parameters = parameters,
        sampling_rate = sampling_rate,
        timestamp = Sys.time()
      ),
      channel_names = channel_names,
      sampling_rate = sampling_rate
    ),
    class = c("preprocessed_emg", "list")
  )

  result
}


#' EMG Preprocessing Preset Configurations
#'
#' Returns preset configurations for common EMG preprocessing scenarios.
#'
#' @param preset Preset name:
#'   \itemize{
#'     \item \code{"default"}: General-purpose EMG preprocessing
#'     \item \code{"gait"}: Optimized for gait/walking analysis
#'     \item \code{"isometric"}: Optimized for isometric contractions
#'     \item \code{"fine_motor"}: Optimized for fine motor tasks
#'   }
#'
#' @return A list containing preprocessing configuration for each step.
#'
#' @examples
#' # Get default configuration
#' config <- emg_preset_config("default")
#'
#' # Modify for custom needs
#' config$notch$freq <- 60  # Change to 60 Hz
#' config$envelope$cutoff <- 10  # Higher cutoff for faster movements
#'
#' @export
emg_preset_config <- function(preset = "default") {
  configs <- list(
    default = list(
      dc_removal = list(method = "mean"),
      bandpass = list(low = 20, high = 450, order = 4),
      notch = list(freq = 50, harmonics = TRUE),
      rectify = list(method = "full_wave"),
      envelope = list(method = "lowpass", cutoff = 6),
      normalize = list(method = "max")
    ),

    gait = list(
      dc_removal = list(method = "linear"),
      bandpass = list(low = 20, high = 450, order = 4),
      notch = list(freq = 50, harmonics = TRUE),
      rectify = list(method = "full_wave"),
      envelope = list(method = "lowpass", cutoff = 4),
      normalize = list(method = "mvc")
    ),

    isometric = list(
      dc_removal = list(method = "mean"),
      bandpass = list(low = 20, high = 450, order = 4),
      notch = list(freq = 50, harmonics = TRUE),
      rectify = list(method = "full_wave"),
      envelope = list(method = "rms", window_ms = 100),
      normalize = list(method = "mvc")
    ),

    fine_motor = list(
      dc_removal = list(method = "mean"),
      bandpass = list(low = 10, high = 500, order = 4),
      notch = list(freq = 50, harmonics = TRUE),
      rectify = list(method = "full_wave"),
      envelope = list(method = "lowpass", cutoff = 10),
      normalize = list(method = "max")
    )
  )

  configs[[preset]]
}


#' Print Method for preprocessed_emg Objects
#'
#' @param x A preprocessed_emg object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.preprocessed_emg <- function(x, ...) {
  cat("Preprocessed EMG Data\n")
  cat("=====================\n\n")

  # Data summary
  if (is.list(x$data) && !is.matrix(x$data)) {
    cat("Trials/Subjects:", length(x$data), "\n")
    dims <- sapply(x$data, dim)
    cat("Time points:", paste(unique(dims[1, ]), collapse = ", "), "\n")
    cat("Channels:", dims[2, 1], "\n")
  } else {
    cat("Time points:", nrow(x$data), "\n")
    cat("Channels:", ncol(x$data), "\n")
  }

  cat("Channel names:", paste(x$channel_names, collapse = ", "), "\n\n")

  # Processing info
  cat("Processing steps:\n")
  for (step in x$processing_info$steps_applied) {
    cat("  -", step, "\n")
  }
  cat("\nOriginal sampling rate:", x$sampling_rate, "Hz\n")
  cat("Processed at:", format(x$processing_info$timestamp), "\n")

  invisible(x)
}


#' Summary Method for preprocessed_emg Objects
#'
#' @param object A preprocessed_emg object
#' @param ... Additional arguments (ignored)
#'
#' @export
summary.preprocessed_emg <- function(object, ...) {
  cat("Preprocessed EMG Data Summary\n")
  cat("=============================\n\n")

  # Data summary
  if (is.list(object$data) && !is.matrix(object$data)) {
    cat("Number of trials/subjects:", length(object$data), "\n")

    # Compute statistics across all data
    all_data <- do.call(rbind, object$data)
  } else {
    cat("Single trial/subject\n")
    all_data <- object$data
  }

  cat("Channels:", length(object$channel_names), "\n")
  cat("Channel names:", paste(object$channel_names, collapse = ", "), "\n\n")

  # Per-channel statistics
  cat("Per-channel statistics:\n")
  stats_df <- data.frame(
    Channel = object$channel_names,
    Mean = apply(all_data, 2, mean, na.rm = TRUE),
    SD = apply(all_data, 2, stats::sd, na.rm = TRUE),
    Min = apply(all_data, 2, min, na.rm = TRUE),
    Max = apply(all_data, 2, max, na.rm = TRUE)
  )
  print(stats_df, row.names = FALSE)

  cat("\nProcessing parameters:\n")
  for (step in names(object$processing_info$parameters)) {
    cat("  ", step, ":\n", sep = "")
    params <- object$processing_info$parameters[[step]]
    for (param_name in names(params)) {
      cat("    ", param_name, ": ", format(params[[param_name]]), "\n", sep = "")
    }
  }

  invisible(object)
}


# Null coalescing operator (if not already defined)
`%||%` <- function(x, y) if (is.null(x)) y else x
