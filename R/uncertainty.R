# =============================================================================
# Empirical Uncertainty Evaluation Functions
# gazeneuro package extension
#
# Purpose: Evaluate theoretical uncertainty predictions (P1-P5) against
#          empirical gaze-tracking recordings
#
# Author: Eric Prince
# Date: 2025-12-13
# =============================================================================

#' @importFrom zoo rollmean
#' @importFrom stats var sd cor lm coef predict quantile median na.omit
#' @importFrom grDevices png dev.off colorRampPalette
#' @importFrom graphics par plot points lines abline text mtext legend axis box barplot hist
#' @importFrom dplyr mutate filter select arrange group_by summarise ungroup n lag lead left_join row_number ntile case_when %>%
#' @importFrom tidyr drop_na
NULL

# =============================================================================
# CONFIGURATION
# =============================================================================

#' Calculate visual angle from physical measurements
#'
#' @description
#' Converts physical screen dimensions and viewing distance to visual angle
#' in degrees. Uses the formula: angle = 2 * atan(size / (2 * distance)) * (180/pi)
#'
#' @param size Physical size (width or height) in any unit
#' @param distance Viewing distance in the SAME unit as size
#'
#' @return Visual angle in degrees
#' @export
calculate_visual_angle <- function(size, distance) {
  angle_rad <- 2 * atan(size / (2 * distance))
  angle_deg <- angle_rad * (180 / pi)
  return(angle_deg)
}


#' Create uncertainty parameters from physical measurements
#'
#' @description
#' Creates analysis parameters by calculating visual angles from physical
#' screen dimensions and viewing distance. All measurements must be in the
#' same unit (cm, inches, etc.).
#'
#' @param screen_width Physical screen width
#' @param screen_height Physical screen height
#' @param viewing_distance Distance from eyes to screen (same units as above)
#' @param sampling_rate_hz Eye tracker sampling rate in Hz
#' @param units Character label for the units used (for documentation only)
#'
#' @return List of analysis parameters with computed visual angles
#' @export
#'
#' @examples
#' # Example: 24" monitor (53cm x 30cm) at 60cm viewing distance
#' params <- create_uncertainty_params(
#'   screen_width = 53,
#'   screen_height = 30,
#'   viewing_distance = 60,
#'   units = "cm"
#' )
create_uncertainty_params <- function(screen_width,
                                      screen_height,
                                      viewing_distance,
                                      sampling_rate_hz = 120,
                                      units = "cm") {

  # Calculate visual angles
  screen_width_deg <- calculate_visual_angle(screen_width, viewing_distance)
  screen_height_deg <- calculate_visual_angle(screen_height, viewing_distance)

  message(sprintf("Display configuration (%s):", units))
  message(sprintf("  Screen: %.1f x %.1f %s", screen_width, screen_height, units))
  message(sprintf("  Viewing distance: %.1f %s", viewing_distance, units))
  message(sprintf("  Visual angle: %.1f° x %.1f°", screen_width_deg, screen_height_deg))
  message(sprintf("  Sampling rate: %d Hz", sampling_rate_hz))

  list(
    # Physical measurements (for reference)
    screen_width_physical = screen_width,
    screen_height_physical = screen_height,
    viewing_distance = viewing_distance,
    units = units,

    # Computed visual angles
    screen_width_deg = screen_width_deg,
    screen_height_deg = screen_height_deg,

    # Eye tracker
    sampling_rate_hz = sampling_rate_hz,

    # Fixation detection
    velocity_threshold_deg_s = 30,
    min_fixation_duration_ms = 100,

    # Timing uncertainty (from your simulation)
    sigma_delta_sec = 0.020,      # 20ms timing uncertainty

    # Drift analysis
    drift_warning_threshold = 0.01,  # 1% drift triggers warning

    # Independence test
    independence_threshold = 0.3,    # |rho| < 0.3 for independence

    # Velocity coupling test
    velocity_coupling_r2_threshold = 0.7
  )
}


#' Default parameters for uncertainty analysis
#'
#' @description
#' Returns default parameters assuming a typical desktop setup.
#' For accurate analysis, use create_uncertainty_params() with your
#' actual physical measurements instead.
#'
#' @return List of analysis parameters
#' @export
default_uncertainty_params <- function() {
  message("Using default display parameters (24\" monitor at 60cm)")
  message("For accurate analysis, use create_uncertainty_params() with your measurements")

  # Typical 24" monitor: ~53cm x 30cm at 60cm viewing distance
  create_uncertainty_params(
    screen_width = 53,
    screen_height = 30,
    viewing_distance = 60,
    sampling_rate_hz = 120,
    units = "cm"
  )
}

# =============================================================================
# PHASE 2: EXTRACTION FUNCTIONS
# =============================================================================

#' Calculate instantaneous gaze velocity
#'
#' @description
#' Computes velocity from consecutive gaze samples. Velocity is calculated
#' in both normalized screen coordinates and degrees of visual angle.
#'
#' @param gaze_data Data frame with gaze_x, gaze_y, and time_sec columns
#' @param params Parameters from default_uncertainty_params()
#'
#' @return Data frame with added velocity columns
#' @export
calculate_gaze_velocity <- function(gaze_data, params = default_uncertainty_params()) {

  gaze_data <- gaze_data %>%
    arrange(time_sec) %>%
    mutate(
      # Time differences
      dt = c(NA, diff(time_sec)),

      # Position differences (normalized 0-1 coordinates)
      dx = c(NA, diff(gaze_x)),
      dy = c(NA, diff(gaze_y)),

      # Displacement magnitude (normalized)
      displacement = sqrt(dx^2 + dy^2),

      # Velocity in normalized units per second
      velocity_norm = displacement / dt,

      # Convert to degrees of visual angle
      dx_deg = dx * params$screen_width_deg,
      dy_deg = dy * params$screen_height_deg,
      displacement_deg = sqrt(dx_deg^2 + dy_deg^2),
      velocity_deg_s = displacement_deg / dt
    )

  # Remove infinite velocities from zero dt

  gaze_data <- gaze_data %>%
    mutate(
      velocity_deg_s = ifelse(is.infinite(velocity_deg_s), NA, velocity_deg_s),
      velocity_norm = ifelse(is.infinite(velocity_norm), NA, velocity_norm)
    )

  message(sprintf("Velocity calculation complete:"))
  message(sprintf("  - Valid samples: %d of %d",
                  sum(!is.na(gaze_data$velocity_deg_s)), nrow(gaze_data)))
  message(sprintf("  - Mean velocity: %.1f °/s",
                  mean(gaze_data$velocity_deg_s, na.rm = TRUE)))
  message(sprintf("  - Max velocity: %.1f °/s",
                  max(gaze_data$velocity_deg_s, na.rm = TRUE)))

  return(gaze_data)
}


#' Detect fixations based on velocity threshold
#'
#' @description
#' Identifies fixation periods where gaze velocity falls below threshold.
#' Also groups consecutive fixation samples into fixation events.
#'
#' @param gaze_data Data frame with velocity_deg_s column
#' @param params Parameters from default_uncertainty_params()
#'
#' @return Data frame with fixation labels
#' @export
detect_fixations <- function(gaze_data, params = default_uncertainty_params()) {

  if (!"velocity_deg_s" %in% names(gaze_data)) {
    gaze_data <- calculate_gaze_velocity(gaze_data, params)
  }

  gaze_data <- gaze_data %>%
    mutate(
      # Binary fixation flag
      is_fixation = !is.na(velocity_deg_s) &
        velocity_deg_s < params$velocity_threshold_deg_s,

      # Create fixation groups (consecutive fixation samples)
      fixation_change = is_fixation != lag(is_fixation, default = FALSE),
      fixation_group = cumsum(fixation_change)
    )

  # Label fixation groups with duration filter
  min_samples <- params$min_fixation_duration_ms / 1000 * params$sampling_rate_hz

  fixation_stats <- gaze_data %>%
    filter(is_fixation) %>%
    group_by(fixation_group) %>%
    summarise(
      n_samples = n(),
      duration_ms = (max(time_sec) - min(time_sec)) * 1000,
      .groups = "drop"
    ) %>%
    mutate(valid_fixation = n_samples >= min_samples)

  # Merge back
  gaze_data <- gaze_data %>%
    left_join(
      fixation_stats %>% select(fixation_group, valid_fixation),
      by = "fixation_group"
    ) %>%
    mutate(
      valid_fixation = replace_na(valid_fixation, FALSE),
      is_fixation = is_fixation & valid_fixation
    )

  n_fixations <- sum(fixation_stats$valid_fixation)
  n_fixation_samples <- sum(gaze_data$is_fixation)
  pct_fixation <- n_fixation_samples / nrow(gaze_data) * 100

  message(sprintf("\nFixation detection complete:"))
  message(sprintf("  - Fixations detected: %d", n_fixations))
  message(sprintf("  - Fixation samples: %d (%.1f%%)", n_fixation_samples, pct_fixation))
  message(sprintf("  - Saccade samples: %d (%.1f%%)",
                  nrow(gaze_data) - n_fixation_samples, 100 - pct_fixation))

  return(gaze_data)
}


#' Estimate tracker noise from fixation jitter
#'
#' @description
#' During stable fixations, position variance reflects hardware noise (σ_tracker).
#' This function estimates tracker noise by analyzing within-fixation variance.
#'
#' @param gaze_data Data frame with fixation labels
#' @param params Parameters from default_uncertainty_params()
#'
#' @return List with tracker noise estimates
#' @export
estimate_tracker_noise <- function(gaze_data, params = default_uncertainty_params()) {

  if (!"is_fixation" %in% names(gaze_data)) {
    gaze_data <- detect_fixations(gaze_data, params)
  }

  # Get fixation samples only
  fixation_data <- gaze_data %>%
    filter(is_fixation)

  if (nrow(fixation_data) < 10) {
    warning("Too few fixation samples for reliable noise estimation")
    return(list(
      sigma_tracker_norm = NA,
      sigma_tracker_deg = NA,
      sigma_x_norm = NA,
      sigma_y_norm = NA,
      n_fixation_samples = nrow(fixation_data),
      method = "insufficient_data"
    ))
  }

  # Method 1: Overall variance during fixations
  # (This overestimates if fixations are at different locations)

  # Method 2: Within-fixation variance (better)
  within_fixation_var <- fixation_data %>%
    group_by(fixation_group) %>%
    filter(n() >= 5) %>%  # Need enough samples per fixation
    summarise(
      var_x = var(gaze_x, na.rm = TRUE),
      var_y = var(gaze_y, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )

  # Weighted average variance (weight by samples)
  if (nrow(within_fixation_var) > 0) {
    total_n <- sum(within_fixation_var$n)
    weighted_var_x <- sum(within_fixation_var$var_x * within_fixation_var$n) / total_n
    weighted_var_y <- sum(within_fixation_var$var_y * within_fixation_var$n) / total_n

    sigma_x_norm <- sqrt(weighted_var_x)
    sigma_y_norm <- sqrt(weighted_var_y)
    sigma_tracker_norm <- sqrt(sigma_x_norm^2 + sigma_y_norm^2)

    # Convert to degrees
    sigma_x_deg <- sigma_x_norm * params$screen_width_deg
    sigma_y_deg <- sigma_y_norm * params$screen_height_deg
    sigma_tracker_deg <- sqrt(sigma_x_deg^2 + sigma_y_deg^2)

    method <- "within_fixation_variance"
  } else {
    # Fallback to overall variance
    sigma_x_norm <- sd(fixation_data$gaze_x, na.rm = TRUE)
    sigma_y_norm <- sd(fixation_data$gaze_y, na.rm = TRUE)
    sigma_tracker_norm <- sqrt(sigma_x_norm^2 + sigma_y_norm^2)

    sigma_x_deg <- sigma_x_norm * params$screen_width_deg
    sigma_y_deg <- sigma_y_norm * params$screen_height_deg
    sigma_tracker_deg <- sqrt(sigma_x_deg^2 + sigma_y_deg^2)

    method <- "overall_fixation_variance"
  }

  result <- list(
    sigma_tracker_norm = sigma_tracker_norm,
    sigma_tracker_deg = sigma_tracker_deg,
    sigma_x_norm = sigma_x_norm,
    sigma_y_norm = sigma_y_norm,
    sigma_x_deg = sigma_x_deg,
    sigma_y_deg = sigma_y_deg,
    n_fixation_samples = nrow(fixation_data),
    n_fixations_analyzed = nrow(within_fixation_var),
    method = method
  )

  message(sprintf("\nTracker noise estimation (%s):", method))
  message(sprintf("  - σ_tracker: %.4f (normalized), %.3f° (visual angle)",
                  sigma_tracker_norm, sigma_tracker_deg))
  message(sprintf("  - σ_x: %.4f norm (%.3f°)", sigma_x_norm, sigma_x_deg))
  message(sprintf("  - σ_y: %.4f norm (%.3f°)", sigma_y_norm, sigma_y_deg))
  message(sprintf("  - Based on %d samples from %d fixations",
                  nrow(fixation_data), nrow(within_fixation_var)))

  return(result)
}


#' Analyze velocity distribution
#'
#' @description
#' Characterizes the velocity distribution to understand viewing behavior:
#' proportion of fixations vs saccades, typical velocities, etc.
#'
#' @param gaze_data Data frame with velocity_deg_s column
#' @param params Parameters from default_uncertainty_params()
#'
#' @return List with velocity distribution statistics
#' @export
analyze_velocity_distribution <- function(gaze_data, params = default_uncertainty_params()) {

  if (!"velocity_deg_s" %in% names(gaze_data)) {
    gaze_data <- calculate_gaze_velocity(gaze_data, params)
  }

  velocities <- gaze_data$velocity_deg_s[!is.na(gaze_data$velocity_deg_s)]

  # Classify
  threshold <- params$velocity_threshold_deg_s
  fixation_velocities <- velocities[velocities < threshold]
  saccade_velocities <- velocities[velocities >= threshold]

  result <- list(
    # Overall statistics
    n_samples = length(velocities),
    mean_velocity = mean(velocities),
    median_velocity = median(velocities),
    sd_velocity = sd(velocities),
    max_velocity = max(velocities),

    # Fixation statistics
    n_fixation = length(fixation_velocities),
    pct_fixation = length(fixation_velocities) / length(velocities) * 100,
    mean_fixation_velocity = mean(fixation_velocities),

    # Saccade statistics
    n_saccade = length(saccade_velocities),
    pct_saccade = length(saccade_velocities) / length(velocities) * 100,
    mean_saccade_velocity = ifelse(length(saccade_velocities) > 0,
                                   mean(saccade_velocities), NA),
    max_saccade_velocity = ifelse(length(saccade_velocities) > 0,
                                  max(saccade_velocities), NA),

    # Raw data for histograms
    velocities = velocities,
    fixation_velocities = fixation_velocities,
    saccade_velocities = saccade_velocities,

    # Quantiles
    quantiles = quantile(velocities, probs = c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99))
  )

  message(sprintf("\nVelocity distribution analysis:"))
  message(sprintf("  - Total samples: %d", result$n_samples))
  message(sprintf("  - Fixation (<%d°/s): %d (%.1f%%), mean %.1f°/s",
                  threshold, result$n_fixation, result$pct_fixation,
                  result$mean_fixation_velocity))
  message(sprintf("  - Saccade (≥%d°/s): %d (%.1f%%), mean %.1f°/s, max %.1f°/s",
                  threshold, result$n_saccade, result$pct_saccade,
                  result$mean_saccade_velocity, result$max_saccade_velocity))

  return(result)
}


#' Analyze slice transitions for boundary uncertainty
#'
#' @description
#' Identifies slice transition boundaries and flags gaze points
#' that fall within the temporal uncertainty window.
#'
#' @param integrated_data Data frame from integrate_all_gaze_points()
#' @param z_axis Z-axis events data frame
#' @param params Parameters from default_uncertainty_params()
#'
#' @return List with transition analysis results
#' @export
analyze_slice_transitions <- function(integrated_data, z_axis,
                                      params = default_uncertainty_params()) {

  sigma_delta <- params$sigma_delta_sec


  # Normalize z_axis timestamps if needed
  z_clean <- z_axis %>%
    arrange(client_timestamp) %>%
    mutate(
      time_sec = (client_timestamp - min(client_timestamp)) / 1000
    )

  # Identify transitions (where slice_index changes)
  transitions <- z_clean %>%
    mutate(
      prev_index = lag(index, default = first(index)),
      is_transition = index != prev_index
    ) %>%
    filter(is_transition)

  if (nrow(transitions) == 0) {
    message("No slice transitions found in z-axis data")
    return(list(
      n_transitions = 0,
      n_ambiguous_gaze = 0,
      pct_ambiguous = 0,
      transition_times = numeric(0)
    ))
  }

  transition_times <- transitions$time_sec

  # For each gaze point, compute distance to nearest transition
  integrated_data <- integrated_data %>%
    mutate(
      dist_to_transition = sapply(time_aligned, function(t) {
        if (length(transition_times) == 0) return(Inf)
        min(abs(t - transition_times))
      }),

      # Classify by confidence level
      assignment_confidence = case_when(
        dist_to_transition > 3 * sigma_delta ~ "High (>99.7%)",
        dist_to_transition > 2 * sigma_delta ~ "Good (95-99.7%)",
        dist_to_transition > 1 * sigma_delta ~ "Moderate (68-95%)",
        TRUE ~ "Ambiguous (<68%)"
      ),

      # Binary flag for ambiguous
      is_ambiguous = dist_to_transition <= 3 * sigma_delta
    )

  # Summarize
  confidence_summary <- integrated_data %>%
    count(assignment_confidence) %>%
    mutate(pct = n / sum(n) * 100) %>%
    arrange(desc(pct))

  n_ambiguous <- sum(integrated_data$is_ambiguous)
  pct_ambiguous <- n_ambiguous / nrow(integrated_data) * 100

  result <- list(
    n_transitions = nrow(transitions),
    transition_times = transition_times,
    n_total_gaze = nrow(integrated_data),
    n_ambiguous_gaze = n_ambiguous,
    pct_ambiguous = pct_ambiguous,
    confidence_summary = confidence_summary,
    sigma_delta_used = sigma_delta,
    integrated_data_with_confidence = integrated_data
  )

  message(sprintf("\nSlice transition analysis:"))
  message(sprintf("  - Transitions detected: %d", nrow(transitions)))
  message(sprintf("  - σ_δ used: %.0f ms", sigma_delta * 1000))
  message(sprintf("  - Ambiguous gaze (within ±3σ): %d (%.2f%%)",
                  n_ambiguous, pct_ambiguous))
  message("\n  Confidence breakdown:")
  for (i in 1:nrow(confidence_summary)) {
    message(sprintf("    %s: %d (%.1f%%)",
                    confidence_summary$assignment_confidence[i],
                    confidence_summary$n[i],
                    confidence_summary$pct[i]))
  }

  return(result)
}


# =============================================================================
# PHASE 3: PREDICTION TEST FUNCTIONS
# =============================================================================

#' Test P1: Velocity-Spatial Coupling
#'
#' @description
#' Tests whether displacement variance scales with velocity, as predicted
#' by the theoretical framework.
#'
#' @param gaze_data Data frame with velocity and displacement columns
#' @param params Parameters from default_uncertainty_params()
#' @param n_bins Number of velocity bins
#'
#' @return List with P1 test results
#' @export
test_p1_velocity_coupling <- function(gaze_data,
                                      params = default_uncertainty_params(),
                                      n_bins = 10) {

  if (!"velocity_deg_s" %in% names(gaze_data)) {
    gaze_data <- calculate_gaze_velocity(gaze_data, params)
  }

  # Filter valid data
  valid_data <- gaze_data %>%
    filter(!is.na(velocity_deg_s), !is.na(displacement))

  # Create velocity bins
  velocity_range <- range(valid_data$velocity_deg_s)
  breaks <- seq(velocity_range[1], velocity_range[2], length.out = n_bins + 1)

  valid_data <- valid_data %>%
    mutate(velocity_bin = cut(velocity_deg_s, breaks = breaks, include.lowest = TRUE))

  # Compute statistics per bin
  binned_stats <- valid_data %>%
    group_by(velocity_bin) %>%
    summarise(
      mean_velocity = mean(velocity_deg_s, na.rm = TRUE),
      var_displacement = var(displacement, na.rm = TRUE),
      sd_displacement = sd(displacement, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    filter(n >= 10)  # Require minimum samples per bin

  if (nrow(binned_stats) < 3) {
    warning("Too few valid bins for P1 analysis")
    return(list(
      pass = FALSE,
      r_squared = NA,
      slope = NA,
      message = "Insufficient data for analysis"
    ))
  }

  # Fit linear model: sd_displacement ~ mean_velocity
  fit <- lm(sd_displacement ~ mean_velocity, data = binned_stats)
  r_squared <- summary(fit)$r.squared
  slope <- coef(fit)[2]

  # Determine pass/fail
  pass <- r_squared >= params$velocity_coupling_r2_threshold && slope > 0

  result <- list(
    pass = pass,
    r_squared = r_squared,
    slope = slope,
    intercept = coef(fit)[1],
    binned_data = binned_stats,
    model = fit,
    threshold = params$velocity_coupling_r2_threshold
  )

  message(sprintf("\n=== P1: Velocity-Spatial Coupling ==="))
  message(sprintf("  R² = %.3f (threshold: %.2f)", r_squared,
                  params$velocity_coupling_r2_threshold))
  message(sprintf("  Slope = %.6f", slope))
  message(sprintf("  Result: %s", ifelse(pass, "✓ PASS", "✗ FAIL")))

  return(result)
}


#' Test P2: Error Independence
#'
#' @description
#' Tests whether different error proxies are uncorrelated, supporting
#' the assumption that variances add directly.
#'
#' @param gaze_data Data frame with gaze data
#' @param params Parameters from default_uncertainty_params()
#'
#' @return List with P2 test results
#' @export
test_p2_independence <- function(gaze_data, params = default_uncertainty_params()) {

  if (!"velocity_deg_s" %in% names(gaze_data)) {
    gaze_data <- calculate_gaze_velocity(gaze_data, params)
  }

  # Compute error proxies
  gaze_data <- gaze_data %>%
    mutate(
      # Proxy 1: High-frequency jitter (deviation from local mean)
      local_mean_x = zoo::rollmean(gaze_x, k = 5, fill = NA, align = "center"),
      local_mean_y = zoo::rollmean(gaze_y, k = 5, fill = NA, align = "center"),
      jitter_x = gaze_x - local_mean_x,
      jitter_y = gaze_y - local_mean_y,
      jitter_magnitude = sqrt(jitter_x^2 + jitter_y^2),

      # Proxy 2: Displacement (velocity-related)
      # Already have displacement from velocity calculation

      # Proxy 3: Cumulative drift from session start
      drift_x = gaze_x - first(gaze_x[!is.na(gaze_x)]),
      drift_y = gaze_y - first(gaze_y[!is.na(gaze_y)]),
      drift_magnitude = sqrt(drift_x^2 + drift_y^2)
    )

  # Filter to valid data
  valid_data <- gaze_data %>%
    filter(!is.na(jitter_magnitude), !is.na(displacement), !is.na(drift_magnitude))

  if (nrow(valid_data) < 50) {
    warning("Too few samples for independence analysis")
    return(list(pass = FALSE, message = "Insufficient data"))
  }

  # Compute correlations
  cor_jitter_displacement <- cor(valid_data$jitter_magnitude,
                                 valid_data$displacement,
                                 use = "complete.obs")
  cor_jitter_drift <- cor(valid_data$jitter_magnitude,
                          valid_data$drift_magnitude,
                          use = "complete.obs")
  cor_displacement_drift <- cor(valid_data$displacement,
                                valid_data$drift_magnitude,
                                use = "complete.obs")

  # Create correlation matrix
  cor_matrix <- matrix(c(
    1, cor_jitter_displacement, cor_jitter_drift,
    cor_jitter_displacement, 1, cor_displacement_drift,
    cor_jitter_drift, cor_displacement_drift, 1
  ), nrow = 3, byrow = TRUE)
  rownames(cor_matrix) <- colnames(cor_matrix) <- c("Jitter", "Displacement", "Drift")

  max_cor <- max(abs(c(cor_jitter_displacement, cor_jitter_drift, cor_displacement_drift)))
  pass <- max_cor < params$independence_threshold

  result <- list(
    pass = pass,
    cor_jitter_displacement = cor_jitter_displacement,
    cor_jitter_drift = cor_jitter_drift,
    cor_displacement_drift = cor_displacement_drift,
    max_correlation = max_cor,
    correlation_matrix = cor_matrix,
    threshold = params$independence_threshold,
    proxy_data = valid_data %>%
      select(time_sec, jitter_magnitude, displacement, drift_magnitude)
  )

  message(sprintf("\n=== P2: Error Independence ==="))
  message(sprintf("  ρ(jitter, displacement) = %.3f", cor_jitter_displacement))
  message(sprintf("  ρ(jitter, drift) = %.3f", cor_jitter_drift))
  message(sprintf("  ρ(displacement, drift) = %.3f", cor_displacement_drift))
  message(sprintf("  max|ρ| = %.3f (threshold: %.2f)", max_cor,
                  params$independence_threshold))
  message(sprintf("  Result: %s", ifelse(pass, "✓ PASS", "✗ FAIL")))

  return(result)
}


#' Test P3: Temporal Uncertainty Accumulation
#'
#' @description
#' Tests whether timing uncertainty grows over the session duration,
#' as predicted by clock drift effects.
#'
#' @param integrated_data Data frame from integrate_all_gaze_points()
#' @param n_segments Number of session segments to analyze
#'
#' @return List with P3 test results
#' @export
test_p3_temporal_accumulation <- function(integrated_data, n_segments = 4) {

  # Compute timing residuals within each z-axis interval
  residual_data <- integrated_data %>%
    group_by(z_event_id) %>%
    mutate(
      interval_midpoint = mean(time_aligned, na.rm = TRUE),
      time_residual = time_aligned - interval_midpoint
    ) %>%
    ungroup()

  # Divide session into segments
  session_duration <- max(residual_data$time_sec, na.rm = TRUE) -
    min(residual_data$time_sec, na.rm = TRUE)

  residual_data <- residual_data %>%
    mutate(
      session_segment = ntile(time_sec, n_segments)
    )

  # Compute residual variance per segment
  variance_by_segment <- residual_data %>%
    group_by(session_segment) %>%
    summarise(
      residual_var = var(time_residual, na.rm = TRUE),
      residual_sd = sd(time_residual, na.rm = TRUE),
      mean_time = mean(time_sec, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )

  # Check if variance grows with time
  first_segment_var <- variance_by_segment$residual_var[1]
  last_segment_var <- variance_by_segment$residual_var[n_segments]
  variance_ratio <- last_segment_var / first_segment_var

  # Fit trend
  if (nrow(variance_by_segment) >= 3) {
    trend_fit <- lm(residual_sd ~ mean_time, data = variance_by_segment)
    trend_slope <- coef(trend_fit)[2]
    trend_positive <- trend_slope > 0
  } else {
    trend_slope <- NA
    trend_positive <- NA
  }

  grows_with_time <- variance_ratio > 1

  result <- list(
    grows_with_time = grows_with_time,
    variance_ratio = variance_ratio,
    first_segment_sd = sqrt(first_segment_var),
    last_segment_sd = sqrt(last_segment_var),
    trend_slope = trend_slope,
    variance_by_segment = variance_by_segment,
    session_duration = session_duration
  )

  message(sprintf("\n=== P3: Temporal Uncertainty Accumulation ==="))
  message(sprintf("  Session duration: %.1f seconds", session_duration))
  message(sprintf("  First segment σ: %.4f s", sqrt(first_segment_var)))
  message(sprintf("  Last segment σ: %.4f s", sqrt(last_segment_var)))
  message(sprintf("  Variance ratio (last/first): %.2f", variance_ratio))
  message(sprintf("  Grows with time: %s",
                  ifelse(grows_with_time, "✓ YES", "✗ NO")))

  return(result)
}


#' Test P4: Probabilistic Slice Assignment
#'
#' @description
#' Quantifies the proportion of gaze points with uncertain slice assignment
#' near transition boundaries.
#'
#' @param integrated_data Data frame from integrate_all_gaze_points()
#' @param z_axis Z-axis events data frame
#' @param params Parameters from default_uncertainty_params()
#'
#' @return List with P4 test results
#' @export
test_p4_slice_assignment <- function(integrated_data, z_axis,
                                     params = default_uncertainty_params()) {

  # Use the transition analysis function
  transition_results <- analyze_slice_transitions(integrated_data, z_axis, params)

  # P4 is informational - report the breakdown
  result <- list(
    n_transitions = transition_results$n_transitions,
    n_total_gaze = transition_results$n_total_gaze,
    n_ambiguous = transition_results$n_ambiguous_gaze,
    pct_ambiguous = transition_results$pct_ambiguous,
    confidence_breakdown = transition_results$confidence_summary,
    sigma_delta = transition_results$sigma_delta_used
  )

  message(sprintf("\n=== P4: Probabilistic Slice Assignment ==="))
  message(sprintf("  Slice transitions: %d", result$n_transitions))
  message(sprintf("  Gaze in ambiguous zone (±3σ): %.2f%%", result$pct_ambiguous))
  message("  (This is informational - no pass/fail criterion)")

  return(result)
}


#' Test P5: Calibration Drift
#'
#' @description
#' Estimates calibration drift rate by comparing mean gaze position
#' across session segments.
#'
#' @param gaze_data Data frame with gaze data
#' @param params Parameters from default_uncertainty_params()
#' @param n_segments Number of session segments
#'
#' @return List with P5 test results
#' @export
test_p5_calibration_drift <- function(gaze_data,
                                      params = default_uncertainty_params(),
                                      n_segments = 4) {

  # Divide session into segments
  gaze_data <- gaze_data %>%
    mutate(session_segment = ntile(time_sec, n_segments))

  # Compute mean position per segment
  segment_means <- gaze_data %>%
    group_by(session_segment) %>%
    summarise(
      mean_x = mean(gaze_x, na.rm = TRUE),
      mean_y = mean(gaze_y, na.rm = TRUE),
      mid_time = mean(time_sec, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )

  # Compute drift from first to last segment
  first_seg <- segment_means[1, ]
  last_seg <- segment_means[n_segments, ]

  drift_x <- last_seg$mean_x - first_seg$mean_x
  drift_y <- last_seg$mean_y - first_seg$mean_y
  total_drift <- sqrt(drift_x^2 + drift_y^2)

  # Compute drift rate
  session_duration <- max(gaze_data$time_sec) - min(gaze_data$time_sec)
  drift_rate_per_sec <- total_drift / session_duration
  drift_rate_per_min <- drift_rate_per_sec * 60

  # Convert to degrees
  drift_x_deg <- drift_x * params$screen_width_deg
  drift_y_deg <- drift_y * params$screen_height_deg
  total_drift_deg <- sqrt(drift_x_deg^2 + drift_y_deg^2)
  drift_rate_deg_per_min <- total_drift_deg / session_duration * 60

  # Fit linear trend
  trend_x <- lm(mean_x ~ mid_time, data = segment_means)
  trend_y <- lm(mean_y ~ mid_time, data = segment_means)

  result <- list(
    drift_x_norm = drift_x,
    drift_y_norm = drift_y,
    total_drift_norm = total_drift,
    drift_x_deg = drift_x_deg,
    drift_y_deg = drift_y_deg,
    total_drift_deg = total_drift_deg,
    drift_rate_norm_per_min = drift_rate_per_min,
    drift_rate_deg_per_min = drift_rate_deg_per_min,
    drift_rate_pct_per_min = drift_rate_per_min * 100,
    session_duration = session_duration,
    segment_means = segment_means,
    trend_x_slope = coef(trend_x)[2],
    trend_y_slope = coef(trend_y)[2]
  )

  message(sprintf("\n=== P5: Calibration Drift ==="))
  message(sprintf("  Session duration: %.1f seconds", session_duration))
  message(sprintf("  Total drift: %.4f (norm), %.3f° (visual angle)",
                  total_drift, total_drift_deg))
  message(sprintf("  Drift rate: %.4f%%/min (%.4f°/min)",
                  drift_rate_per_min * 100, drift_rate_deg_per_min))

  if (drift_rate_per_min > params$drift_warning_threshold) {
    message(sprintf("  ⚠ WARNING: Drift rate exceeds %.1f%%/min threshold",
                    params$drift_warning_threshold * 100))
  }

  return(result)
}


# =============================================================================
# PHASE 4: COMPLETE ANALYSIS WORKFLOW
# =============================================================================

#' Evaluate uncertainty in gaze-slice integration
#'
#' @description
#' Executes all extraction and prediction test functions on an empirical recording.
#'
#' @param gaze_data Raw gaze data frame
#' @param z_axis Z-axis events data frame
#' @param integrated_data Output from integrate_all_gaze_points()
#' @param params Optional parameters (uses defaults if NULL)
#'
#' @return List containing all analysis results
#' @export
evaluate_uncertainty <- function(gaze_data, z_axis, integrated_data,
                                 params = NULL) {

  if (is.null(params)) {
    params <- default_uncertainty_params()
  }

  message("=" %s+% strrep("=", 60))
  message("UNCERTAINTY EVALUATION")
  message("=" %s+% strrep("=", 60))
  message(sprintf("\nGaze samples: %d", nrow(gaze_data)))
  message(sprintf("Z-axis events: %d", nrow(z_axis)))
  message(sprintf("Integrated samples: %d", nrow(integrated_data)))

  # Phase 2: Extractions
  message("\n" %s+% strrep("-", 40))
  message("PHASE 2: Parameter Extraction")
  message(strrep("-", 40))

  # Calculate velocity and detect fixations
  gaze_data <- calculate_gaze_velocity(gaze_data, params)
  gaze_data <- detect_fixations(gaze_data, params)

  # Estimate tracker noise
  tracker_noise <- estimate_tracker_noise(gaze_data, params)

  # Velocity distribution
  velocity_dist <- analyze_velocity_distribution(gaze_data, params)

  # Slice transitions
  transitions <- analyze_slice_transitions(integrated_data, z_axis, params)

  # Phase 3: Prediction Tests
  message("\n" %s+% strrep("-", 40))
  message("PHASE 3: Prediction Tests")
  message(strrep("-", 40))

  p1_result <- test_p1_velocity_coupling(gaze_data, params)
  p2_result <- test_p2_independence(gaze_data, params)
  p3_result <- test_p3_temporal_accumulation(integrated_data)
  p4_result <- test_p4_slice_assignment(integrated_data, z_axis, params)
  p5_result <- test_p5_calibration_drift(gaze_data, params)

  # Compile results
  results <- list(
    params = params,

    # Extracted parameters
    tracker_noise = tracker_noise,
    velocity_distribution = velocity_dist,
    transitions = transitions,

    # Prediction test results
    p1_velocity_coupling = p1_result,
    p2_independence = p2_result,
    p3_temporal_accumulation = p3_result,
    p4_slice_assignment = p4_result,
    p5_calibration_drift = p5_result,

    # Processed data
    processed_gaze = gaze_data,

    # Summary
    summary = list(
      p1_pass = p1_result$pass,
      p2_pass = p2_result$pass,
      p3_grows = p3_result$grows_with_time,
      p4_pct_ambiguous = p4_result$pct_ambiguous,
      p5_drift_rate = p5_result$drift_rate_pct_per_min
    )
  )

  # Print summary
  message("\n" %s+% strrep("=", 60))
  message("SUMMARY")
  message(strrep("=", 60))
  message(sprintf("\nP1 Velocity Coupling:     %s (R²=%.3f)",
                  ifelse(p1_result$pass, "✓ PASS", "✗ FAIL"), p1_result$r_squared))
  message(sprintf("P2 Error Independence:    %s (max|ρ|=%.3f)",
                  ifelse(p2_result$pass, "✓ PASS", "✗ FAIL"), p2_result$max_correlation))
  message(sprintf("P3 Temporal Accumulation: %s (ratio=%.2f)",
                  ifelse(p3_result$grows_with_time, "✓ Grows", "— Stable"),
                  p3_result$variance_ratio))
  message(sprintf("P4 Slice Assignment:      %.2f%% ambiguous", p4_result$pct_ambiguous))
  message(sprintf("P5 Calibration Drift:     %.4f%%/min", p5_result$drift_rate_pct_per_min))

  message(sprintf("\nTracker noise (σ_tracker): %.4f norm (%.3f°)",
                  tracker_noise$sigma_tracker_norm, tracker_noise$sigma_tracker_deg))

  return(results)
}


# =============================================================================
# PHASE 5: VISUALIZATION
# =============================================================================

#' Visualize uncertainty evaluation results
#'
#' @description
#' Creates a 6-panel summary figure showing uncertainty analysis results.
#'
#' @param results Output from evaluate_uncertainty()
#' @param output_file Path for output PNG (NULL for screen display)
#' @param width Image width in pixels
#' @param height Image height in pixels
#'
#' @export
visualize_uncertainty <- function(results, output_file = NULL,
                                  width = 1200, height = 900) {

  if (!is.null(output_file)) {
    png(output_file, width = width, height = height)
  }

  par(mfrow = c(2, 3), mar = c(4, 4, 3, 2), oma = c(0, 0, 3, 0))

  # Panel 1: Velocity Distribution
  hist(results$velocity_distribution$velocities,
       breaks = 50, col = "lightblue", border = "white",
       main = "Velocity Distribution",
       xlab = "Velocity (°/s)", ylab = "Frequency")
  abline(v = results$params$velocity_threshold_deg_s, col = "red", lwd = 2, lty = 2)
  legend("topright", legend = c("Fixation threshold"),
         col = "red", lty = 2, lwd = 2, bty = "n", cex = 0.8)

  # Panel 2: P1 - Velocity vs Displacement
  p1 <- results$p1_velocity_coupling
  if (!is.null(p1$binned_data) && nrow(p1$binned_data) > 0) {
    plot(p1$binned_data$mean_velocity, p1$binned_data$sd_displacement,
         pch = 19, col = "steelblue", cex = 1.5,
         main = sprintf("P1: Velocity Coupling (R²=%.3f)", p1$r_squared),
         xlab = "Mean Velocity (°/s)", ylab = "Displacement SD")
    if (!is.na(p1$slope)) {
      abline(p1$model, col = "red", lwd = 2)
    }
    status_col <- ifelse(p1$pass, "darkgreen", "red")
    legend("topleft", legend = ifelse(p1$pass, "✓ PASS", "✗ FAIL"),
           text.col = status_col, bty = "n", cex = 1.2)
  } else {
    plot.new()
    text(0.5, 0.5, "Insufficient data\nfor P1 analysis", cex = 1.2)
  }

  # Panel 3: P2 - Independence Scatter
  p2 <- results$p2_independence
  if (!is.null(p2$proxy_data)) {
    plot(p2$proxy_data$jitter_magnitude, p2$proxy_data$displacement,
         pch = ".", col = rgb(0, 0, 0, 0.1),
         main = sprintf("P2: Independence (ρ=%.3f)", p2$cor_jitter_displacement),
         xlab = "Jitter Magnitude", ylab = "Displacement")
    status_col <- ifelse(p2$pass, "darkgreen", "red")
    legend("topleft", legend = ifelse(p2$pass, "✓ PASS", "✗ FAIL"),
           text.col = status_col, bty = "n", cex = 1.2)
  } else {
    plot.new()
    text(0.5, 0.5, "Insufficient data\nfor P2 analysis", cex = 1.2)
  }

  # Panel 4: P4 - Assignment Confidence
  p4 <- results$p4_slice_assignment
  if (!is.null(p4$confidence_breakdown)) {
    conf_colors <- c("High (>99.7%)" = "darkgreen",
                     "Good (95-99.7%)" = "lightgreen",
                     "Moderate (68-95%)" = "orange",
                     "Ambiguous (<68%)" = "red")
    bar_colors <- conf_colors[p4$confidence_breakdown$assignment_confidence]
    barplot(p4$confidence_breakdown$pct,
            names.arg = gsub(" ", "\n", p4$confidence_breakdown$assignment_confidence),
            col = bar_colors,
            main = "P4: Slice Assignment Confidence",
            ylab = "% of Gaze Points", las = 1, cex.names = 0.7)
  }

  # Panel 5: P5 - Drift Over Session
  p5 <- results$p5_calibration_drift
  if (!is.null(p5$segment_means)) {
    plot(p5$segment_means$mid_time, p5$segment_means$mean_x,
         type = "b", pch = 19, col = "steelblue", lwd = 2,
         main = sprintf("P5: Calibration Drift (%.4f%%/min)", p5$drift_rate_pct_per_min),
         xlab = "Session Time (s)", ylab = "Mean Gaze X (normalized)")
    abline(lm(mean_x ~ mid_time, data = p5$segment_means), col = "red", lty = 2)
  }

  # Panel 6: Summary
  plot.new()
  text(0.5, 0.95, "EMPIRICAL VALIDATION", font = 2, cex = 1.3)
  text(0.5, 0.85, "SUMMARY", font = 2, cex = 1.1)

  y_pos <- 0.70
  line_height <- 0.12

  # P1
  p1_text <- sprintf("P1 Velocity Coupling: %s (R²=%.3f)",
                     ifelse(results$summary$p1_pass, "✓", "✗"),
                     results$p1_velocity_coupling$r_squared)
  text(0.05, y_pos, p1_text, adj = 0,
       col = ifelse(results$summary$p1_pass, "darkgreen", "red"))

  # P2
  y_pos <- y_pos - line_height
  p2_text <- sprintf("P2 Independence: %s (max|ρ|=%.3f)",
                     ifelse(results$summary$p2_pass, "✓", "✗"),
                     results$p2_independence$max_correlation)
  text(0.05, y_pos, p2_text, adj = 0,
       col = ifelse(results$summary$p2_pass, "darkgreen", "red"))

  # P3
  y_pos <- y_pos - line_height
  p3_text <- sprintf("P3 Temporal: %s (ratio=%.2f)",
                     ifelse(results$summary$p3_grows, "↑ Grows", "— Stable"),
                     results$p3_temporal_accumulation$variance_ratio)
  text(0.05, y_pos, p3_text, adj = 0)

  # P4
  y_pos <- y_pos - line_height
  p4_text <- sprintf("P4 Ambiguous: %.2f%%", results$summary$p4_pct_ambiguous)
  text(0.05, y_pos, p4_text, adj = 0)

  # P5
  y_pos <- y_pos - line_height
  p5_text <- sprintf("P5 Drift: %.4f%%/min", results$summary$p5_drift_rate)
  text(0.05, y_pos, p5_text, adj = 0)

  # Tracker noise
  y_pos <- y_pos - line_height * 1.5
  text(0.05, y_pos, sprintf("σ_tracker: %.4f (%.3f°)",
                            results$tracker_noise$sigma_tracker_norm,
                            results$tracker_noise$sigma_tracker_deg), adj = 0)

  # Overall title
  mtext("Uncertainty Evaluation", outer = TRUE, cex = 1.5, line = 1)

  if (!is.null(output_file)) {
    dev.off()
    message(sprintf("\nVisualization saved to: %s", output_file))
  }
}


# =============================================================================
# HELPER: String concatenation operator
# =============================================================================
`%s+%` <- function(a, b) paste0(a, b)


# =============================================================================
# USAGE EXAMPLE
# =============================================================================

# Example usage (uncomment and modify paths):
#
# library(gazeneuro)
# library(tidyverse)
#
# # Create params with your physical measurements
# params <- create_uncertainty_params(
#   screen_width = 53,
#   screen_height = 30,
#   viewing_distance = 60,
#   sampling_rate_hz = 120,
#   units = "cm"
# )
#
# # Load data for case_12 sagittal
# nifti_data <- preload_nifti_data("path/to/case_12/sagittal.nii.gz")
# gaze_data <- read_csv("path/to/case_12/gaze_data.csv")
# z_axis <- read_csv("path/to/case_12/sagittal_slice_events.csv")
#
# # Integrate (using your existing function)
# integrated <- integrate_all_gaze_points(gaze_data, z_axis)
#
# # Run uncertainty evaluation
# results <- evaluate_uncertainty(gaze_data, z_axis, integrated, params)
#
# # Generate visualization
# visualize_uncertainty(results, "case_12_uncertainty_evaluation.png")
