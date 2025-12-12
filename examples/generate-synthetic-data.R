# Synthetic Data Generator for gazeneuro Experiment 1
# Generates eye-tracking data for three modalities: research tracker, webcam, and reference

library(tidyverse)

#' Generate synthetic viewing session with proper timestamp alignment
#'
#' @param duration_seconds Total session duration in seconds (default: 30)
#' @param n_slices Total number of slices in the scan (default: 25)
#' @param n_scroll_events Number of scroll events during session (default: 20)
#' @param sampling_rate Gaze sampling rate in Hz (default: 60)
#' @param spatial_noise_sd Standard deviation of spatial noise in normalized coords
#' @param latency_seconds System latency in seconds (default: 0.05 = 50ms)
#' @param temporal_jitter Proportion of timing jitter (e.g., 0.05 = ±5%)
#' @param seed Random seed for reproducibility
#'
#' @return List with gaze_data and z_axis data frames (all timestamps in seconds)
generate_synthetic_session <- function(
    duration_seconds = 30,
    n_slices = 25,
    n_scroll_events = 20,
    sampling_rate = 60,
    spatial_noise_sd = 0.02,
    latency_seconds = 0.05,
    temporal_jitter = 0.02,
    seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)

  # Step 1: Generate z-axis events
  # Events occur: at start, during scrolling, and at end

  # Generate event times throughout the session
  # Start with beginning and end
  event_times <- c(0, duration_seconds)

  # Add scroll events in between
  if (n_scroll_events > 0) {
    # Generate random times when user scrolls
    scroll_times <- sort(runif(n_scroll_events,
                               min = 0.5,  # Don't scroll immediately
                               max = duration_seconds - 0.5))  # Leave time at end
    event_times <- sort(unique(c(event_times, scroll_times)))
  }

  n_events <- length(event_times)

  # For each event, determine which slice is being viewed
  # Simulate realistic scrolling behavior
  current_slice_indices <- numeric(n_events)
  current_slice_indices[1] <- runif(1, 0, 0.3)  # Start near beginning

  for (i in 2:n_events) {
    # Simulate different scrolling patterns
    scroll_type <- sample(c("next", "previous", "jump", "nearby"), 1,
                          prob = c(0.4, 0.2, 0.1, 0.3))

    prev_index <- current_slice_indices[i-1]

    if (scroll_type == "next") {
      # Scroll forward by 1-3 slices
      delta <- runif(1, 1, 3) / n_slices
      new_index <- prev_index + delta
    } else if (scroll_type == "previous") {
      # Scroll back by 1-3 slices
      delta <- runif(1, 1, 3) / n_slices
      new_index <- prev_index - delta
    } else if (scroll_type == "jump") {
      # Jump to a different part
      new_index <- runif(1, 0, 1)
    } else {  # nearby
      # Small adjustment
      delta <- runif(1, -2, 2) / n_slices
      new_index <- prev_index + delta
    }

    # Clamp to valid range
    current_slice_indices[i] <- pmax(0, pmin(1, new_index))
  }

  # Generate z_axis data frame
  z_axis <- data.frame(
    client_timestamp = event_times * 1000,  # Convert to milliseconds
    time_seconds = event_times,
    plane = "AXIAL",
    index = current_slice_indices,
    image_id = "synthetic_scan_001",
    slice_number = round(current_slice_indices * (n_slices - 1)) + 1,
    stringsAsFactors = FALSE
  )

  # Step 2: Generate gaze data
  # Calculate sampling times with optional jitter
  sampling_interval <- 1.0 / sampling_rate
  n_samples <- floor(duration_seconds * sampling_rate)

  # Generate gaze timestamps (in seconds)
  gaze_times_seconds <- seq(0, duration_seconds - sampling_interval, length.out = n_samples)

  # Add temporal jitter if requested
  if (temporal_jitter > 0) {
    jitter <- runif(n_samples, -temporal_jitter, temporal_jitter) * sampling_interval
    gaze_times_seconds <- gaze_times_seconds + jitter
    # Ensure times stay in bounds
    gaze_times_seconds <- pmax(0, pmin(duration_seconds - sampling_interval, gaze_times_seconds))
    # Sort to maintain temporal order
    gaze_times_seconds <- sort(gaze_times_seconds)
  }

  # Apply latency: gaze data arrives AFTER what's displayed
  gaze_times_with_latency <- gaze_times_seconds + latency_seconds

  # For each gaze point, determine which slice was being viewed
  # Use the z-axis events to determine current slice at each time
  slice_at_gaze_time <- numeric(n_samples)

  for (i in 1:n_samples) {
    # Find which z-axis event was active at this time (accounting for latency)
    actual_time <- gaze_times_seconds[i]  # Time when the slice was actually viewed
    event_idx <- findInterval(actual_time, event_times)
    event_idx <- max(1, min(event_idx, n_events))
    slice_at_gaze_time[i] <- current_slice_indices[event_idx]
  }

  # Generate gaze coordinates based on viewed content
  gaze_x <- numeric(n_samples)
  gaze_y <- numeric(n_samples)

  # Track fixation points for each unique slice
  slice_fixations <- list()

  for (i in 1:n_samples) {
    slice_key <- as.character(round(slice_at_gaze_time[i] * 1000))  # Key for this slice

    # Get or create fixation point for this slice
    if (!(slice_key %in% names(slice_fixations))) {
      # New slice - create a fixation point
      slice_fixations[[slice_key]] <- list(
        x = runif(1, 0.3, 0.7),
        y = runif(1, 0.3, 0.7),
        visit_count = 0
      )
    }

    fixation <- slice_fixations[[slice_key]]
    fixation$visit_count <- fixation$visit_count + 1

    # Add drift and microsaccades
    drift_amount <- 0.02 * sin(2 * pi * i / sampling_rate)  # Smooth drift

    gaze_x[i] <- fixation$x + drift_amount
    gaze_y[i] <- fixation$y + drift_amount * 0.7
  }

  # Add spatial noise (microsaccades and measurement error)
  if (spatial_noise_sd > 0) {
    gaze_x <- gaze_x + rnorm(n_samples, 0, spatial_noise_sd)
    gaze_y <- gaze_y + rnorm(n_samples, 0, spatial_noise_sd)
  }

  # Clamp to valid range [0, 1]
  gaze_x <- pmax(0, pmin(1, gaze_x))
  gaze_y <- pmax(0, pmin(1, gaze_y))

  # Create gaze data frame
  gaze_data <- data.frame(
    device_time_stamp = gaze_times_with_latency * 1e6,  # Convert to microseconds
    time_seconds = gaze_times_with_latency,
    gaze_point_on_display_area_x = gaze_x,
    gaze_point_on_display_area_y = gaze_y,
    tracking_session = paste0("synthetic_", sampling_rate, "Hz"),
    true_slice_index = slice_at_gaze_time,  # Which slice was actually being viewed
    stringsAsFactors = FALSE
  )

  # Keep only samples within session duration (some may exceed due to latency)
  gaze_data <- gaze_data[gaze_data$time_seconds <= duration_seconds, ]

  # Print summary for debugging
  cat("Synthetic Data Generation Summary:\n")
  cat("─────────────────────────────────\n")
  cat(sprintf("Duration: %.1f seconds\n", duration_seconds))
  cat(sprintf("Gaze samples: %d at %d Hz\n", nrow(gaze_data), sampling_rate))
  cat(sprintf("Z-axis events: %d events\n", nrow(z_axis)))
  cat(sprintf("Unique slices viewed: %d\n", length(unique(z_axis$slice_number))))
  cat(sprintf("Latency: %.0f ms\n", latency_seconds * 1000))
  cat(sprintf("Gaze time range: %.3f - %.3f seconds\n",
              min(gaze_data$time_seconds), max(gaze_data$time_seconds)))
  cat(sprintf("Z-axis time range: %.3f - %.3f seconds\n",
              min(z_axis$time_seconds), max(z_axis$time_seconds)))

  return(list(
    gaze_data = gaze_data,
    z_axis = z_axis,
    metadata = list(
      duration = duration_seconds,
      n_slices = n_slices,
      n_events = n_events,
      sampling_rate = sampling_rate,
      spatial_noise_sd = spatial_noise_sd,
      latency_seconds = latency_seconds,
      temporal_jitter = temporal_jitter
    )
  ))
}

# Visualization function to check the synthetic data
visualize_synthetic_session <- function(synthetic_data) {
  library(ggplot2)
  library(dplyr)

  # Extract data
  gaze <- synthetic_data$gaze_data
  z_axis <- synthetic_data$z_axis

  # Create a combined timeline plot
  # Sample gaze data to avoid overplotting
  gaze_sample <- gaze[seq(1, nrow(gaze), by = max(1, floor(nrow(gaze)/500))), ]

  # Plot 1: Timeline with both data sources
  p1 <- ggplot() +
    # Z-axis events as vertical lines
    geom_vline(data = z_axis,
               aes(xintercept = time_seconds),
               color = "red", alpha = 0.5, linetype = "dashed") +
    # Gaze points
    geom_point(data = gaze_sample,
               aes(x = time_seconds, y = 1),
               color = "blue", alpha = 0.3, size = 0.5) +
    # Z-axis points
    geom_point(data = z_axis,
               aes(x = time_seconds, y = 2),
               color = "red", size = 3) +
    scale_y_continuous(breaks = c(1, 2),
                       labels = c("Gaze Data", "Z-axis Events"),
                       limits = c(0.5, 2.5)) +
    labs(title = "Synthetic Session Timeline",
         subtitle = "Red lines show when slices changed",
         x = "Time (seconds)",
         y = "") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank())

  # Plot 2: Slice indices over time
  # Create a step function for slice viewing
  slice_steps <- data.frame(
    time = c(z_axis$time_seconds, synthetic_data$metadata$duration),
    slice = c(z_axis$slice_number, tail(z_axis$slice_number, 1))
  )

  p2 <- ggplot() +
    geom_step(data = slice_steps,
              aes(x = time, y = slice),
              color = "darkgreen", size = 1) +
    geom_point(data = z_axis,
               aes(x = time_seconds, y = slice_number),
               color = "red", size = 3) +
    labs(title = "Slice Viewing Pattern",
         subtitle = "Shows which slice was displayed over time",
         x = "Time (seconds)",
         y = "Slice Number") +
    theme_minimal()

  # Plot 3: Gaze heatmap
  p3 <- ggplot(gaze_sample,
               aes(x = gaze_point_on_display_area_x,
                   y = gaze_point_on_display_area_y)) +
    stat_density_2d(aes(fill = ..density..),
                    geom = "tile",
                    contour = FALSE,
                    n = 50) +
    scale_fill_viridis_c() +
    coord_fixed() +
    labs(title = "Gaze Distribution Heatmap",
         x = "Gaze X", y = "Gaze Y") +
    theme_minimal()

  # Combine plots
  library(gridExtra)
  grid.arrange(p1, p2, p3, layout_matrix = rbind(c(1, 1), c(2, 3)))
}

# Test function
test_synthetic_data <- function() {
  # Generate synthetic data
  synthetic <- generate_synthetic_session(
    duration_seconds = 30,
    n_slices = 25,
    n_scroll_events = 15,
    sampling_rate = 60,
    spatial_noise_sd = 0.02,
    latency_seconds = 0.050,
    seed = 123
  )

  # Visualize
  visualize_synthetic_session(synthetic)

  return(synthetic)
}

#' Generate multiple sessions for a specific modality
#'
#' @param modality One of "research", "webcam", or "reference"
#' @param n_sessions Number of sessions to generate
#' @param base_seed Base random seed
#'
#' @return List of synthetic sessions
generate_modality_sessions <- function(modality, n_sessions = 100, base_seed = 12345) {

  # Define modality parameters based on Methods section
  params <- switch(modality,
                   research = list(
                     sampling_rate = 120,
                     spatial_noise_sd = 0.5,  # 0.5° visual angle
                     latency_seconds = 0.05,
                     temporal_jitter = 0.05  # ±5%
                   ),
                   webcam = list(
                     sampling_rate = 20,
                     spatial_noise_sd = 2.0,  # 2.0° visual angle
                     latency_seconds = 0.15,
                     temporal_jitter = 0.05
                   ),
                   reference = list(
                     sampling_rate = 60,
                     spatial_noise_sd = 0,
                     latency_seconds = 0,
                     temporal_jitter = 0
                   ),
                   stop("Unknown modality")
  )

  # Generate sessions
  sessions <- list()

  for (i in 1:n_sessions) {
    message(sprintf("Generating %s session %d/%d", modality, i, n_sessions))

    sessions[[i]] <- generate_synthetic_session(
      duration_seconds = 300,  # 5 minutes
      n_slices = 25,
      sampling_rate = params$sampling_rate,
      spatial_noise_sd = params$spatial_noise_sd,
      latency_seconds = params$latency_seconds,
      temporal_jitter = params$temporal_jitter,
      seed = base_seed + i
    )

    # Add modality label
    sessions[[i]]$modality <- modality
    sessions[[i]]$session_id <- i
  }

  return(sessions)
}

#' Validate synthetic data characteristics
#'
#' @param session A synthetic session from generate_synthetic_session_fixed()
#' @return Data frame with validation metrics
validate_synthetic_data <- function(session) {
  gaze <- session$gaze_data
  z_axis <- session$z_axis

  # Calculate actual metrics
  # Account for the fact that timestamps don't start at 0
  gaze_start <- min(gaze$device_time_stamp)
  gaze_end <- max(gaze$device_time_stamp)
  duration_actual <- (gaze_end - gaze_start) / 1e6  # Convert from microseconds to seconds

  n_samples <- nrow(gaze)
  actual_rate <- n_samples / duration_actual

  # Calculate time differences to check sampling regularity
  time_diffs <- diff(gaze$device_time_stamp) / 1e6
  mean_interval <- mean(time_diffs)
  sd_interval <- sd(time_diffs)

  # Spatial spread
  spatial_spread_x <- sd(gaze$gaze_point_on_display_area_x)
  spatial_spread_y <- sd(gaze$gaze_point_on_display_area_y)

  # Number of slice events
  n_slices <- nrow(z_axis)

  # Calculate temporal alignment between gaze and z_axis
  gaze_start_sec <- gaze_start / 1e6
  z_start_sec <- min(z_axis$client_timestamp) / 1000
  apparent_latency_sec = (z_start_sec - gaze_start_sec)


  # Return validation metrics
  data.frame(
    expected_rate = session$metadata$sampling_rate,
    actual_rate = round(actual_rate, 1),
    rate_error = round(abs(actual_rate - session$metadata$sampling_rate), 2),
    mean_interval_ms = round(mean_interval * 1000, 2),
    interval_sd_ms = round(sd_interval * 1000, 2),
    spatial_spread_x = round(spatial_spread_x, 4),
    spatial_spread_y = round(spatial_spread_y, 4),
    n_samples = n_samples,
    n_slices = n_slices,
    duration_sec = round(duration_actual, 1),
    true_latency_sec = session$metadata$latency_seconds,
    apparent_latency_sec = round(apparent_latency_sec, 2),
    latency_check = round(abs(apparent_latency_sec - session$metadata$latency_seconds), 2)
  )
}
