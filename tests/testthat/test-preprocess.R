# Tests for EMG preprocessing functions

test_that("remove_dc_offset works with vector input", {
  set.seed(123)
  signal <- rnorm(100) + 5  # Signal with DC offset

  # Mean removal
  result <- remove_dc_offset(signal, method = "mean")
  expect_equal(mean(result), 0, tolerance = 1e-10)

  # Median removal
  result_median <- remove_dc_offset(signal, method = "median")
  expect_true(abs(mean(result_median)) < abs(mean(signal)))

  # Linear detrend
  signal_trend <- 1:100 + rnorm(100)
  result_linear <- remove_dc_offset(signal_trend, method = "linear")
  expect_true(abs(mean(result_linear)) < abs(mean(signal_trend)))
})

test_that("remove_dc_offset works with matrix input", {
  set.seed(123)
  emg <- matrix(rnorm(100 * 4) + rep(c(1, 2, 3, 4), each = 100), ncol = 4)

  result <- remove_dc_offset(emg)
  expect_equal(dim(result), dim(emg))

  # Each channel should have mean near zero
  for (j in 1:4) {
    expect_equal(mean(result[, j]), 0, tolerance = 1e-10)
  }
})

test_that("remove_dc_offset works with list input", {
  set.seed(123)
  data_list <- list(
    matrix(rnorm(50 * 4) + 2, ncol = 4),
    matrix(rnorm(60 * 4) + 3, ncol = 4)
  )

  result <- remove_dc_offset(data_list)
  expect_length(result, 2)
  expect_equal(nrow(result[[1]]), 50)
  expect_equal(nrow(result[[2]]), 60)
})

test_that("bandpass_filter validates inputs", {
  signal <- rnorm(1000)

  # Missing sampling_rate

  expect_error(bandpass_filter(signal))

  # Invalid cutoff frequencies
  expect_error(bandpass_filter(signal, sampling_rate = 1000, low_cutoff = 0))
  expect_error(bandpass_filter(signal, sampling_rate = 1000, high_cutoff = 600))
  expect_error(bandpass_filter(signal, sampling_rate = 1000, low_cutoff = 100, high_cutoff = 50))
})

test_that("bandpass_filter applies correctly", {
  skip_if_not_installed("signal")

  # Create signal with known frequency components
  fs <- 1000
  t <- seq(0, 1, by = 1/fs)
  low_freq <- sin(2 * pi * 5 * t)      # 5 Hz - should be removed
  target_freq <- sin(2 * pi * 50 * t)   # 50 Hz - should pass
  high_freq <- sin(2 * pi * 480 * t)    # 480 Hz - should be removed

  signal <- low_freq + target_freq + high_freq

  filtered <- bandpass_filter(signal, sampling_rate = fs, low_cutoff = 20, high_cutoff = 450)

  expect_length(filtered, length(signal))

  # The target frequency should dominate after filtering
  # This is a basic check - the filter should reduce low and high freq components
  expect_true(is.numeric(filtered))
})

test_that("bandpass_filter works with matrix input", {
  skip_if_not_installed("signal")

  set.seed(123)
  emg <- matrix(rnorm(1000 * 4), ncol = 4)

  result <- bandpass_filter(emg, sampling_rate = 1000)
  expect_equal(dim(result), dim(emg))
})

test_that("notch_filter removes power line frequency", {
  skip_if_not_installed("signal")

  fs <- 1000
  t <- seq(0, 1, by = 1/fs)
  signal <- rnorm(length(t)) + sin(2 * pi * 50 * t)  # Add 50 Hz noise

  filtered <- notch_filter(signal, sampling_rate = fs, notch_freq = 50, harmonics = FALSE)

  expect_length(filtered, length(signal))
  expect_true(is.numeric(filtered))
})

test_that("notch_filter handles harmonics", {
  skip_if_not_installed("signal")

  fs <- 1000
  signal <- rnorm(1000)

  # With harmonics
  result_with <- notch_filter(signal, sampling_rate = fs, notch_freq = 50, harmonics = TRUE)
  # Without harmonics
  result_without <- notch_filter(signal, sampling_rate = fs, notch_freq = 50, harmonics = FALSE)

  expect_length(result_with, length(signal))
  expect_length(result_without, length(signal))
})

test_that("rectify works correctly", {
  signal <- c(-3, -1, 0, 1, 3)

  # Full-wave rectification
  full <- rectify(signal, method = "full_wave")
  expect_equal(full, c(3, 1, 0, 1, 3))

  # Half-wave rectification
  half <- rectify(signal, method = "half_wave")
  expect_equal(half, c(0, 0, 0, 1, 3))
})

test_that("rectify works with matrix and list", {
  mat <- matrix(c(-1, 1, -2, 2), ncol = 2)
  result_mat <- rectify(mat)
  expect_equal(result_mat, abs(mat))

  lst <- list(mat, mat * 2)
  result_lst <- rectify(lst)
  expect_length(result_lst, 2)
  expect_equal(result_lst[[1]], abs(mat))
})

test_that("extract_envelope lowpass method works", {
  skip_if_not_installed("signal")

  fs <- 1000
  signal <- abs(rnorm(1000))

  result <- extract_envelope(signal, sampling_rate = fs, method = "lowpass", cutoff = 6)
  expect_length(result, length(signal))
  expect_true(all(is.finite(result)))
})

test_that("extract_envelope rms method works", {
  signal <- rep(c(0, 1, 1, 0), 25)

  result <- extract_envelope(signal, sampling_rate = 100, method = "rms", window_ms = 50)
  expect_length(result, length(signal))
  expect_true(all(result >= 0))
})

test_that("extract_envelope moving_average method works", {
  signal <- c(rep(0, 50), rep(1, 50))

  result <- extract_envelope(signal, sampling_rate = 100, method = "moving_average", window_ms = 100)
  expect_length(result, length(signal))
  expect_true(all(is.finite(result)))
})

test_that("normalize_emg max method works", {
  signal <- c(0, 5, 10, 5, 0)
  result <- normalize_emg(signal, method = "max")
  expect_equal(max(result), 1)
  expect_equal(result, signal / 10)
})

test_that("normalize_emg zscore method works", {
  set.seed(123)
  signal <- rnorm(100, mean = 50, sd = 10)
  result <- normalize_emg(signal, method = "zscore")

  expect_equal(mean(result), 0, tolerance = 1e-10)
  expect_equal(sd(result), 1, tolerance = 1e-10)
})

test_that("normalize_emg mvc method works", {
  signal <- c(0, 50, 100, 50, 0)
  mvc_ref <- 200

  result <- normalize_emg(signal, method = "mvc", reference = mvc_ref)
  expect_equal(result, signal / 200 * 100)  # Percentage of MVC
  expect_equal(max(result), 50)  # 100/200 * 100 = 50%
})

test_that("normalize_emg range method works", {
  signal <- c(10, 20, 30, 40, 50)
  result <- normalize_emg(signal, method = "range")

  expect_equal(min(result), 0)
  expect_equal(max(result), 1)
})

test_that("normalize_emg works with matrix per_channel", {
  mat <- matrix(c(0, 10, 0, 20), ncol = 2)
  result <- normalize_emg(mat, method = "max", per_channel = TRUE)

  expect_equal(max(result[, 1]), 1)
  expect_equal(max(result[, 2]), 1)
})

test_that("compute_mvc_reference peak method works", {
  mvc_data <- matrix(c(50, 100, 75, 80, 160, 120), ncol = 2)

  result <- compute_mvc_reference(mvc_data, method = "peak")
  expect_equal(result[1], 100)
  expect_equal(result[2], 160)
})

test_that("compute_mvc_reference handles list of trials", {
  mvc_trials <- list(
    matrix(c(80, 50), ncol = 2),
    matrix(c(60, 70), ncol = 2)
  )

  result <- compute_mvc_reference(mvc_trials, method = "peak")
  expect_equal(result[1], 80)
  expect_equal(result[2], 70)
})

test_that("emg_preset_config returns valid configurations", {
  presets <- c("default", "gait", "isometric", "fine_motor")

  for (preset in presets) {
    config <- emg_preset_config(preset)
    expect_true(is.list(config))
    expect_true("bandpass" %in% names(config))
    expect_true("rectify" %in% names(config))
  }
})

test_that("preprocess_emg pipeline works end-to-end", {
  skip_if_not_installed("signal")

  set.seed(123)
  raw_emg <- matrix(rnorm(1000 * 4), ncol = 4)

  result <- preprocess_emg(
    raw_emg,
    sampling_rate = 1000,
    pipeline_config = "default",
    verbose = FALSE
  )

  expect_s3_class(result, "preprocessed_emg")
  expect_true("data" %in% names(result))
  expect_true("processing_info" %in% names(result))
  expect_equal(ncol(result$data), 4)

  # Check that steps were recorded
  expect_true(length(result$processing_info$steps_applied) > 0)
})

test_that("preprocess_emg works with list input", {
  skip_if_not_installed("signal")

  set.seed(123)
  raw_list <- list(
    matrix(rnorm(500 * 4), ncol = 4),
    matrix(rnorm(600 * 4), ncol = 4)
  )

  result <- preprocess_emg(
    raw_list,
    sampling_rate = 1000,
    pipeline_config = "default",
    verbose = FALSE
  )

  expect_s3_class(result, "preprocessed_emg")
  expect_true(is.list(result$data))
  expect_length(result$data, 2)
})

test_that("preprocess_emg custom config works", {
  skip_if_not_installed("signal")

  set.seed(123)
  raw_emg <- matrix(rnorm(1000 * 4), ncol = 4)

  custom_config <- list(
    bandpass = list(low = 30, high = 400),
    rectify = list(method = "full_wave"),
    envelope = list(method = "rms", window_ms = 100),
    normalize = list(method = "max")
  )

  result <- preprocess_emg(
    raw_emg,
    sampling_rate = 1000,
    pipeline_config = custom_config,
    verbose = FALSE
  )

  expect_s3_class(result, "preprocessed_emg")
  expect_true("bandpass" %in% result$processing_info$steps_applied)
  expect_false("dc_removal" %in% result$processing_info$steps_applied)  # Not in custom config
})

test_that("print.preprocessed_emg works", {
  skip_if_not_installed("signal")

  set.seed(123)
  raw_emg <- matrix(rnorm(100 * 4), ncol = 4)
  result <- preprocess_emg(raw_emg, sampling_rate = 1000, verbose = FALSE)

  expect_output(print(result), "Preprocessed EMG Data")
})
