#!/usr/bin/env Rscript
# ==============================================================================
# SRCV DEBUGGING SCRIPT
# ==============================================================================
# This script tests each component of the SRCV implementation to identify
# why the validation study is failing
# ==============================================================================

library(tidyverse)
library(signal)

# Source the main functions (adjust path as needed)
# source("srcv_validation_study.R")  # Uncomment and adjust path

# ==============================================================================
# TEST 1: SIMPLE SYNTHETIC SIGNAL TEST
# ==============================================================================

cat("\n")
cat("=" %.% rep("", 80), "\n", sep="")
cat("TEST 1: SIMPLE SYNTHETIC SIGNAL WITH KNOWN LATENCY\n")
cat("=" %.% rep("", 80), "\n", sep="")

# Create a very simple test case
test1_simple_case <- function() {
  # Parameters
  fs <- 100  # 100 Hz sampling
  duration <- 10  # 10 seconds
  true_latency_ms <- 500  # 500ms latency - should be easy to detect

  # Create time vector
  dt <- 1000 / fs  # 10ms per sample
  time_ms <- seq(0, duration * 1000, by = dt)
  n_samples <- length(time_ms)

  # Create simple event signal (3 events)
  event_signal <- rep(0, n_samples)
  event_times_ms <- c(1000, 4000, 7000)  # Events at 1s, 4s, 7s
  event_duration_ms <- 1000  # 1 second events

  for (et in event_times_ms) {
    idx <- which(time_ms >= et & time_ms < (et + event_duration_ms))
    event_signal[idx] <- 1
  }

  # Create gaze response (delayed version of events)
  gaze_signal <- rep(0, n_samples)
  for (et in event_times_ms) {
    # Response starts after latency
    response_start <- et + true_latency_ms
    response_end <- response_start + event_duration_ms
    idx <- which(time_ms >= response_start & time_ms < response_end)
    if (length(idx) > 0) {
      gaze_signal[idx] <- 1
    }
  }

  # Add small noise to gaze
  gaze_signal <- gaze_signal + rnorm(n_samples, 0, 0.05)

  # Plot the signals
  par(mfrow = c(3, 1), mar = c(3, 4, 2, 1))

  # Plot 1: Event signal
  plot(time_ms/1000, event_signal, type = "l", col = "red", lwd = 2,
       main = "Event Signal", xlab = "", ylab = "Signal", ylim = c(-0.1, 1.1))
  grid()

  # Plot 2: Gaze signal
  plot(time_ms/1000, gaze_signal, type = "l", col = "blue", lwd = 2,
       main = sprintf("Gaze Signal (True Latency = %d ms)", true_latency_ms),
       xlab = "", ylab = "Signal")
  grid()

  # Plot 3: Both overlaid
  plot(time_ms/1000, event_signal, type = "l", col = "red", lwd = 2,
       main = "Overlay: Red=Events, Blue=Gaze",
       xlab = "Time (seconds)", ylab = "Signal", ylim = range(c(event_signal, gaze_signal)))
  lines(time_ms/1000, gaze_signal, col = "blue", lwd = 2)
  grid()
  legend("topright", c("Events", "Gaze"), col = c("red", "blue"), lwd = 2)

  par(mfrow = c(1, 1))

  return(list(
    time_ms = time_ms,
    event_signal = event_signal,
    gaze_signal = gaze_signal,
    true_latency_ms = true_latency_ms,
    fs = fs
  ))
}

# Run Test 1
test1_data <- test1_simple_case()

# ==============================================================================
# TEST 2: CROSS-CORRELATION FUNCTION
# ==============================================================================

cat("\n")
cat("=" %.% rep("", 80), "\n", sep="")
cat("TEST 2: TESTING CROSS-CORRELATION FUNCTION\n")
cat("=" %.% rep("", 80), "\n", sep="")

# Simple cross-correlation implementation for comparison
simple_cross_correlation <- function(signal1, signal2, max_lag) {
  n <- length(signal1)
  lags <- -max_lag:max_lag
  correlations <- numeric(length(lags))

  # Normalize signals
  sig1_norm <- (signal1 - mean(signal1)) / sd(signal1)
  sig2_norm <- (signal2 - mean(signal2)) / sd(signal2)

  for (i in seq_along(lags)) {
    lag <- lags[i]
    if (lag < 0) {
      # signal2 leads signal1
      idx1 <- 1:(n + lag)
      idx2 <- (-lag + 1):n
    } else if (lag > 0) {
      # signal1 leads signal2
      idx1 <- (lag + 1):n
      idx2 <- 1:(n - lag)
    } else {
      idx1 <- 1:n
      idx2 <- 1:n
    }

    if (length(idx1) > 0 && length(idx2) > 0) {
      correlations[i] <- sum(sig1_norm[idx1] * sig2_norm[idx2]) / length(idx1)
    }
  }

  return(list(
    lags = lags,
    correlations = correlations
  ))
}

# Test with our simple signals
max_lag_samples <- 100  # 1 second at 100 Hz
xcorr_result <- simple_cross_correlation(
  test1_data$gaze_signal,
  test1_data$event_signal,
  max_lag_samples
)

# Convert lag to milliseconds
lag_ms <- xcorr_result$lags * (1000 / test1_data$fs)

# Find peak
peak_idx <- which.max(xcorr_result$correlations)
estimated_latency_ms <- lag_ms[peak_idx]

# Plot correlation function
plot(lag_ms, xcorr_result$correlations, type = "l", lwd = 2,
     main = "Cross-Correlation Function",
     xlab = "Lag (ms)", ylab = "Correlation")
abline(v = test1_data$true_latency_ms, col = "green", lty = 2, lwd = 2)
abline(v = estimated_latency_ms, col = "red", lty = 2, lwd = 2)
abline(h = 0, col = "gray")
grid()
legend("topright",
       c(sprintf("True Latency: %d ms", test1_data$true_latency_ms),
         sprintf("Estimated: %.1f ms", estimated_latency_ms),
         sprintf("Error: %.1f ms", abs(estimated_latency_ms - test1_data$true_latency_ms))),
       col = c("green", "red", "black"),
       lty = c(2, 2, NA),
       lwd = 2)

cat(sprintf("\nResults:\n"))
cat(sprintf("  True latency: %d ms\n", test1_data$true_latency_ms))
cat(sprintf("  Estimated latency: %.1f ms\n", estimated_latency_ms))
cat(sprintf("  Error: %.1f ms\n", abs(estimated_latency_ms - test1_data$true_latency_ms)))
cat(sprintf("  Peak correlation: %.3f\n", max(xcorr_result$correlations)))

# ==============================================================================
# TEST 3: TEST WITH ACTUAL SRCV FUNCTIONS
# ==============================================================================

cat("\n")
cat("=" %.% rep("", 80), "\n", sep="")
cat("TEST 3: TESTING WITH ACTUAL SRCV FUNCTIONS\n")
cat("=" %.% rep("", 80), "\n", sep="")

# If you have sourced the main script, uncomment this section:
if (exists("compute_cross_correlation")) {
  # Test with the actual function
  xcorr_actual <- compute_cross_correlation(
    test1_data$gaze_signal,
    test1_data$event_signal,
    max_lag_samples
  )

  cat("\nActual SRCV function results:\n")
  cat(sprintf("  Optimal lag: %d samples\n", xcorr_actual$optimal_lag))
  cat(sprintf("  Estimated latency: %.1f ms\n",
              xcorr_actual$optimal_lag * (1000 / test1_data$fs)))
  cat(sprintf("  Max correlation: %.3f\n", xcorr_actual$max_correlation))

  # Compare the two implementations
  plot(xcorr_result$lags, xcorr_result$correlations, type = "l", lwd = 2,
       main = "Comparison: Simple (black) vs SRCV Function (red)",
       xlab = "Lag (samples)", ylab = "Correlation")
  lines(xcorr_actual$lags, xcorr_actual$correlations, col = "red", lwd = 2)
  legend("topright", c("Simple Implementation", "SRCV Function"),
         col = c("black", "red"), lwd = 2)
} else {
  cat("\nNote: SRCV functions not loaded. Source the main script to test them.\n")
}

# ==============================================================================
# TEST 4: REALISTIC SIMULATION TEST
# ==============================================================================

cat("\n")
cat("=" %.% rep("", 80), "\n", sep="")
cat("TEST 4: MORE REALISTIC SIMULATION\n")
cat("=" %.% rep("", 80), "\n", sep="")

test4_realistic <- function() {
  if (!exists("generate_viewing_events") || !exists("generate_gaze_response")) {
    cat("Skipping Test 4: Functions not loaded\n")
    return(NULL)
  }

  # Generate events
  set.seed(123)
  events <- generate_viewing_events(duration = 30, n_events = 10)

  cat("\nGenerated events:\n")
  print(head(events))

  # Test with different SNRs and latencies
  test_cases <- expand.grid(
    snr_db = c(0, 10, 20),
    true_latency = c(0, 100, 500),
    stringsAsFactors = FALSE
  )

  results <- data.frame()

  for (i in 1:nrow(test_cases)) {
    # Generate gaze response
    gaze_data <- generate_gaze_response(
      events = events,
      snr_db = test_cases$snr_db[i],
      fs = 120,
      latency = test_cases$true_latency[i]
    )

    # Estimate latency
    if (exists("estimate_latency_srcv")) {
      result <- estimate_latency_srcv(gaze_data)

      results <- rbind(results, data.frame(
        snr_db = test_cases$snr_db[i],
        true_latency = test_cases$true_latency[i],
        estimated_latency = result$estimated_latency,
        error = result$error,
        correlation_peak = result$correlation_peak,
        prominence_ratio = result$prominence_ratio
      ))
    }
  }

  if (nrow(results) > 0) {
    cat("\nTest Results:\n")
    print(results)

    # Plot errors
    plot(results$true_latency, results$error,
         col = factor(results$snr_db),
         pch = 19, cex = 2,
         main = "Estimation Error vs True Latency",
         xlab = "True Latency (ms)",
         ylab = "Error (ms)")
    legend("topright",
           paste("SNR =", unique(results$snr_db), "dB"),
           col = 1:length(unique(results$snr_db)),
           pch = 19)
    abline(h = 50, col = "red", lty = 2)  # Success threshold
    grid()
  }

  return(results)
}

test4_results <- test4_realistic()

# ==============================================================================
# TEST 5: SIGNAL GENERATION DIAGNOSTIC
# ==============================================================================

cat("\n")
cat("=" %.% rep("", 80), "\n", sep="")
cat("TEST 5: SIGNAL GENERATION DIAGNOSTIC\n")
cat("=" %.% rep("", 80), "\n", sep="")

if (exists("generate_viewing_events") && exists("generate_gaze_response")) {
  # Generate one example with clear parameters
  set.seed(456)
  events_test <- generate_viewing_events(duration = 20, n_events = 5)
  gaze_test <- generate_gaze_response(
    events = events_test,
    snr_db = 20,
    fs = 120,
    latency = 200
  )

  # Plot and examine
  par(mfrow = c(2, 1), mar = c(3, 4, 2, 1))

  # Plot 1: Full signal
  plot(gaze_test$time/1000, gaze_test$gaze, type = "l",
       main = "Generated Gaze Signal (Full)",
       xlab = "", ylab = "Gaze Signal")
  for (i in 1:nrow(events_test)) {
    abline(v = events_test$start_time[i]/1000, col = "red", lty = 2, lwd = 0.5)
    abline(v = (events_test$start_time[i] + 200)/1000, col = "blue", lty = 2, lwd = 0.5)
  }
  legend("topright", c("Event Start", "Expected Response (200ms later)"),
         col = c("red", "blue"), lty = 2)

  # Plot 2: Zoomed on first event
  event1_time <- events_test$start_time[1]
  zoom_range <- which(gaze_test$time >= (event1_time - 500) &
                        gaze_test$time <= (event1_time + 2000))

  if (length(zoom_range) > 0) {
    plot(gaze_test$time[zoom_range], gaze_test$gaze[zoom_range], type = "l",
         main = "Zoomed: First Event Region",
         xlab = "Time (ms)", ylab = "Gaze Signal")
    abline(v = event1_time, col = "red", lty = 2, lwd = 2)
    abline(v = event1_time + 200, col = "blue", lty = 2, lwd = 2)
    grid()
  }

  par(mfrow = c(1, 1))

  # Check statistics
  cat("\nSignal Statistics:\n")
  cat(sprintf("  Time range: %.1f - %.1f ms\n",
              min(gaze_test$time), max(gaze_test$time)))
  cat(sprintf("  Number of samples: %d\n", length(gaze_test$time)))
  cat(sprintf("  Sampling rate: %.1f Hz\n", gaze_test$fs))
  cat(sprintf("  Mean: %.3f, SD: %.3f\n",
              mean(gaze_test$gaze), sd(gaze_test$gaze)))

  cat("\nEvent Statistics:\n")
  cat(sprintf("  Number of events: %d\n", nrow(events_test)))
  cat(sprintf("  Event time range: %.1f - %.1f ms\n",
              min(events_test$start_time), max(events_test$end_time)))
  cat(sprintf("  Average duration: %.1f ms\n", mean(events_test$duration)))
}

# ==============================================================================
# TEST 6: LAG DIRECTION TEST
# ==============================================================================

cat("\n")
cat("=" %.% rep("", 80), "\n", sep="")
cat("TEST 6: LAG DIRECTION TEST\n")
cat("=" %.% rep("", 80), "\n", sep="")

test_lag_direction <- function() {
  # Create two simple signals to test lag direction
  n <- 1000
  signal_a <- rep(0, n)
  signal_b <- rep(0, n)

  # Put a pulse in signal_a at position 200
  signal_a[200:250] <- 1

  # Put a pulse in signal_b at position 300 (lags by 100 samples)
  signal_b[300:350] <- 1

  # Test correlation
  xcorr <- simple_cross_correlation(signal_b, signal_a, 200)

  peak_lag <- xcorr$lags[which.max(xcorr$correlations)]

  par(mfrow = c(3, 1), mar = c(3, 4, 2, 1))

  # Signals
  plot(signal_a, type = "l", main = "Signal A (Reference)",
       xlab = "", ylab = "Signal", col = "blue", lwd = 2)
  plot(signal_b, type = "l", main = "Signal B (Lags A by 100 samples)",
       xlab = "", ylab = "Signal", col = "red", lwd = 2)

  # Correlation
  plot(xcorr$lags, xcorr$correlations, type = "l",
       main = sprintf("Cross-Correlation (Peak at lag = %d)", peak_lag),
       xlab = "Lag", ylab = "Correlation", lwd = 2)
  abline(v = peak_lag, col = "green", lty = 2)
  abline(v = 100, col = "red", lty = 2)  # Expected

  par(mfrow = c(1, 1))

  cat(sprintf("\nLag Direction Test:\n"))
  cat(sprintf("  Signal B lags Signal A by 100 samples\n"))
  cat(sprintf("  Detected lag: %d samples\n", peak_lag))
  cat(sprintf("  Interpretation: %s\n",
              ifelse(peak_lag > 0, "Positive lag = B lags A ✓",
                     "Negative lag = B leads A ✗")))

  return(peak_lag)
}

lag_test_result <- test_lag_direction()

# ==============================================================================
# TEST 7: EDGE CASES
# ==============================================================================

cat("\n")
cat("=" %.% rep("", 80), "\n", sep="")
cat("TEST 7: EDGE CASES\n")
cat("=" %.% rep("", 80), "\n", sep="")

# Test with no events
test_no_correlation <- function() {
  n <- 1000
  signal1 <- rnorm(n)
  signal2 <- rnorm(n)

  xcorr <- simple_cross_correlation(signal1, signal2, 50)

  cat("Test: Uncorrelated signals\n")
  cat(sprintf("  Max correlation: %.3f (should be near 0)\n",
              max(abs(xcorr$correlations))))

  # Test with perfect correlation
  signal3 <- signal1  # Perfect copy
  xcorr_perfect <- simple_cross_correlation(signal1, signal3, 50)

  cat("\nTest: Perfectly correlated signals\n")
  cat(sprintf("  Max correlation: %.3f (should be near 1)\n",
              max(xcorr_perfect$correlations)))
  cat(sprintf("  Lag at max: %d (should be 0)\n",
              xcorr_perfect$lags[which.max(xcorr_perfect$correlations)]))
}

test_no_correlation()

# ==============================================================================
# SUMMARY AND RECOMMENDATIONS
# ==============================================================================

cat("\n")
cat("=" %.% rep("", 80), "\n", sep="")
cat("DEBUGGING SUMMARY\n")
cat("=" %.% rep("", 80), "\n", sep="")

cat("\nKey Checks:\n")
cat("1. Does simple cross-correlation detect the 500ms latency in Test 1?\n")
cat("2. Do the SRCV functions match the simple implementation?\n")
cat("3. Are errors reasonable for different SNRs in Test 4?\n")
cat("4. Is the lag direction correct in Test 6?\n")
cat("5. Are signals properly normalized (non-zero variance)?\n")

cat("\nCommon Issues to Look For:\n")
cat("- Time unit mismatches (seconds vs milliseconds)\n")
cat("- Array indexing errors (off-by-one)\n")
cat("- Incorrect lag sign (positive vs negative)\n")
cat("- Signal normalization destroying timing\n")
cat("- Events outside the time range of gaze signal\n")

cat("\nNext Steps:\n")
cat("1. If Test 1 works but actual SRCV doesn't: Problem is in SRCV implementation\n")
cat("2. If Test 1 fails: Problem is fundamental to correlation approach\n")
cat("3. If only realistic tests fail: Problem is in signal generation\n")
cat("4. Check warnings() for any numerical issues\n")

cat("\n" %.% rep("=", 80), "\n", sep="")
