#' Simulation Framework for gazeneuro Validation
#'
#' @description
#' Provides tools for generating synthetic gaze and z-axis data with
#' controllable uncertainty sources to validate the theoretical framework.
#'
#' The four uncertainty sources modeled:
#' - σ²_tracker: Eye tracker hardware noise
#' - σ²_drift: Calibration drift over time
#' - σ²_induced: Velocity-induced temporal-spatial coupling
#' - σ²_transform: Coordinate transformation error
#'
#' @name simulation
NULL

# =============================================================================
# SIMULATION PARAMETERS
# =============================================================================

#' Default simulation parameters
#'
#' @export
default_sim_params <- function() {

  list(
    # Timing parameters
    sampling_rate_hz = 120,          # Eye tracker sampling rate
    session_duration_sec = 60,       # Total session duration


    # Uncertainty parameters (in appropriate units)
    sigma_tracker = 0.005,           # Tracker noise (normalized coords, ~0.5% of screen)
    sigma_drift_rate = 0.001,        # Drift rate per second (normalized)
    sigma_latency = 0.020,           # Latency uncertainty (seconds, 20ms)
    sigma_transform = 1e-14,         # Transform error (voxels, machine precision)

    # Clock drift
    alpha_clock = 1e-5,              # Clock drift rate (dimensionless, ~10ppm)

    # Gaze behavior parameters
    fixation_duration_mean = 0.300,  # Mean fixation duration (seconds)
    fixation_duration_sd = 0.150,    # SD of fixation duration
    saccade_velocity_mean = 300,     # Mean saccade velocity (deg/sec)
    saccade_velocity_sd = 100,       # SD of saccade velocity
    saccade_duration = 0.030,        # Typical saccade duration (seconds)

    # Slice viewing parameters
    n_slices = 25,                   # Number of slices in volume
    slice_dwell_mean = 2.0,          # Mean time on each slice (seconds)
    slice_dwell_sd = 1.0,            # SD of slice dwell time

    # Image dimensions
    dims = c(512, 512, 25),
    pixel_dims = c(0.4296875, 0.4296875, 5.0)
  )
}

# =============================================================================
# GAZE STREAM GENERATION
# =============================================================================

#' Generate synthetic gaze stream
#'
#' @description
#' Creates a realistic gaze stream with fixations and saccades,
#' including controllable noise sources.
#'
#' @param params List of simulation parameters (from default_sim_params())
#' @param add_tracker_noise Add tracker hardware noise
#' @param add_drift Add calibration drift
#' @param seed Random seed for reproducibility
#'
#' @return Data frame with columns: device_time_stamp, gaze_point_on_display_area_x,
#'         gaze_point_on_display_area_y, true_x, true_y, velocity, is_fixation
#' @export
simulate_gaze_stream <- function(params = default_sim_params(),
                                 add_tracker_noise = TRUE,
                                 add_drift = TRUE,
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_samples <- params$sampling_rate_hz * params$session_duration_sec
  dt <- 1 / params$sampling_rate_hz

  # Initialize storage
  time_sec <- seq(0, params$session_duration_sec - dt, by = dt)
  true_x <- numeric(n_samples)
  true_y <- numeric(n_samples)
  velocity <- numeric(n_samples)
  is_fixation <- logical(n_samples)

  # Generate alternating fixations and saccades
  current_x <- runif(1, 0.2, 0.8)
  current_y <- runif(1, 0.2, 0.8)
  t <- 0
  i <- 1

  while (i <= n_samples) {
    # Generate fixation
    fix_duration <- max(0.1, rnorm(1, params$fixation_duration_mean, params$fixation_duration_sd))
    fix_samples <- min(round(fix_duration * params$sampling_rate_hz), n_samples - i + 1)

    if (fix_samples > 0) {
      idx <- i:(i + fix_samples - 1)
      true_x[idx] <- current_x
      true_y[idx] <- current_y
      velocity[idx] <- 0
      is_fixation[idx] <- TRUE
      i <- i + fix_samples
    }

    if (i > n_samples) break

    # Generate saccade to new location
    target_x <- runif(1, 0.1, 0.9)
    target_y <- runif(1, 0.1, 0.9)
    saccade_amp <- sqrt((target_x - current_x)^2 + (target_y - current_y)^2)

    sac_velocity <- max(100, rnorm(1, params$saccade_velocity_mean, params$saccade_velocity_sd))
    sac_samples <- max(2, round(params$saccade_duration * params$sampling_rate_hz))
    sac_samples <- min(sac_samples, n_samples - i + 1)

    if (sac_samples > 0) {
      idx <- i:(i + sac_samples - 1)
      # Linear interpolation during saccade
      t_frac <- seq(0, 1, length.out = sac_samples)
      true_x[idx] <- current_x + t_frac * (target_x - current_x)
      true_y[idx] <- current_y + t_frac * (target_y - current_y)
      velocity[idx] <- sac_velocity
      is_fixation[idx] <- FALSE

      current_x <- target_x
      current_y <- target_y
      i <- i + sac_samples
    }
  }

  # Apply tracker noise
  observed_x <- true_x
  observed_y <- true_y

  if (add_tracker_noise) {
    observed_x <- observed_x + rnorm(n_samples, 0, params$sigma_tracker)
    observed_y <- observed_y + rnorm(n_samples, 0, params$sigma_tracker)
  }

  # Apply calibration drift
  if (add_drift) {
    drift_x <- cumsum(rnorm(n_samples, 0, params$sigma_drift_rate / sqrt(params$sampling_rate_hz)))
    drift_y <- cumsum(rnorm(n_samples, 0, params$sigma_drift_rate / sqrt(params$sampling_rate_hz)))
    observed_x <- observed_x + drift_x
    observed_y <- observed_y + drift_y
  }

  # Clip to valid range
  observed_x <- pmax(0, pmin(1, observed_x))
  observed_y <- pmax(0, pmin(1, observed_y))

  # Create timestamps in microseconds (like real eye tracker)
  device_time_stamp <- time_sec * 1e6

  data.frame(
    device_time_stamp = device_time_stamp,
    gaze_point_on_display_area_x = observed_x,
    gaze_point_on_display_area_y = observed_y,
    true_x = true_x,
    true_y = true_y,
    velocity = velocity,
    is_fixation = is_fixation,
    time_sec = time_sec
  )
}

# =============================================================================
# Z-AXIS EVENT GENERATION
# =============================================================================

#' Generate synthetic z-axis (slice viewing) events
#'
#' @description
#' Creates a stream of slice navigation events simulating
#' a radiologist scrolling through slices.
#'
#' @param params List of simulation parameters
#' @param gaze_data Gaze data frame (to align durations)
#' @param latency_ms True latency between systems (milliseconds)
#' @param seed Random seed for reproducibility
#'
#' @return Data frame with columns: client_timestamp, plane, index, image_id,
#'         true_slice, event_time_sec
#' @export
simulate_z_axis_events <- function(params = default_sim_params(),
                                   gaze_data = NULL,
                                   latency_ms = 50,
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Determine session duration
  if (!is.null(gaze_data)) {
    duration <- max(gaze_data$time_sec)
  } else {
    duration <- params$session_duration_sec
  }

  # Generate slice viewing sequence
  events <- list()
  t <- 0
  event_id <- 1
  current_slice <- ceiling(params$n_slices / 2)  # Start at middle slice

  while (t < duration) {
    # Record current slice
    events[[event_id]] <- data.frame(
      event_time_sec = t,
      slice_num = current_slice,
      stringsAsFactors = FALSE
    )

    # Dwell time on this slice
    dwell <- max(0.5, rnorm(1, params$slice_dwell_mean, params$slice_dwell_sd))
    t <- t + dwell
    event_id <- event_id + 1

    # Navigate to next slice (random walk with boundaries)
    step <- sample(c(-2, -1, 1, 2), 1, prob = c(0.1, 0.4, 0.4, 0.1))
    current_slice <- max(1, min(params$n_slices, current_slice + step))
  }

  events_df <- do.call(rbind, events)

  # Convert to expected format
  data.frame(
    client_timestamp = (events_df$event_time_sec + latency_ms / 1000) * 1000,  # ms
    plane = "AXIAL",
    index = (events_df$slice_num - 1) / (params$n_slices - 1),  # Normalized 0-1
    image_id = "simulated_image",
    true_slice = events_df$slice_num,
    event_time_sec = events_df$event_time_sec,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# UNCERTAINTY SOURCE INJECTION
# =============================================================================

#' Inject velocity-induced temporal-spatial error (Prediction 1)
#'
#' @description
#' Applies the theoretical prediction: σ_induced ≈ ||v_gaze|| · σ_δ
#' where temporal uncertainty creates larger spatial errors during fast movements.
#'
#' @param gaze_data Gaze data frame with velocity column
#' @param sigma_latency Temporal uncertainty (seconds)
#' @param dt Sampling interval (seconds)
#'
#' @return Gaze data with induced spatial error
#' @export
inject_velocity_induced_error <- function(gaze_data, sigma_latency = 0.020, dt = 1/120) {
  n <- nrow(gaze_data)

  # Temporal offset drawn from latency uncertainty
  time_offset <- rnorm(n, 0, sigma_latency)

  # Approximate position shift using velocity
  # Convert velocity from deg/s to normalized units/s (assume ~60 deg = full screen)
  v_normalized <- gaze_data$velocity / 60

  # Direction of movement (use finite differences where available)
  dx <- c(0, diff(gaze_data$true_x))
  dy <- c(0, diff(gaze_data$true_y))
  mag <- sqrt(dx^2 + dy^2)
  mag[mag == 0] <- 1  # Avoid division by zero

  # Induced spatial error
  spatial_error_x <- v_normalized * time_offset * (dx / mag)
  spatial_error_y <- v_normalized * time_offset * (dy / mag)

  # Apply error
  gaze_data$gaze_point_on_display_area_x <- gaze_data$gaze_point_on_display_area_x + spatial_error_x
  gaze_data$gaze_point_on_display_area_y <- gaze_data$gaze_point_on_display_area_y + spatial_error_y

  # Store for analysis
  gaze_data$induced_error_x <- spatial_error_x
  gaze_data$induced_error_y <- spatial_error_y
  gaze_data$induced_error_magnitude <- sqrt(spatial_error_x^2 + spatial_error_y^2)

  gaze_data
}

#' Inject clock drift (Prediction 3)
#'
#' @description
#' Applies temporal uncertainty accumulation: σ_δ(t) ≈ √(σ²_δ,0 + (α_clock·t)²)
#'
#' @param z_axis_data Z-axis data frame
#' @param sigma_delta_0 Initial temporal uncertainty (seconds)
#' @param alpha_clock Clock drift rate (dimensionless)
#'
#' @return Z-axis data with accumulated drift
#' @export
inject_clock_drift <- function(z_axis_data, sigma_delta_0 = 0.001, alpha_clock = 1e-5) {
  t <- z_axis_data$event_time_sec

  # Accumulated uncertainty
  sigma_t <- sqrt(sigma_delta_0^2 + (alpha_clock * t)^2)

  # Apply drift to timestamps
  drift <- rnorm(nrow(z_axis_data), 0, sigma_t)
  z_axis_data$client_timestamp <- z_axis_data$client_timestamp + drift * 1000  # Convert to ms

  z_axis_data$accumulated_sigma <- sigma_t
  z_axis_data$applied_drift <- drift

  z_axis_data
}

# =============================================================================
# GROUND TRUTH AND VALIDATION
# =============================================================================

#' Calculate expected total uncertainty (Prediction 2)
#'
#' @description
#' Computes σ²_total ≈ σ²_tracker + σ²_drift + σ²_induced + σ²_transform
#' assuming weak correlations between error sources.
#'
#' @param params Simulation parameters
#' @param t Time since session start (seconds)
#' @param velocity Gaze velocity (deg/s)
#'
#' @return List with individual and total uncertainty estimates
#' @export
calculate_expected_uncertainty <- function(params = default_sim_params(),
                                           t = 30,
                                           velocity = 0) {
  # Tracker uncertainty (constant)
  sigma2_tracker <- params$sigma_tracker^2

  # Drift uncertainty (time-dependent)
  sigma2_drift <- (params$sigma_drift_rate * t)^2

  # Velocity-induced uncertainty
  v_normalized <- velocity / 60  # Convert to normalized units/s
  sigma2_induced <- (v_normalized * params$sigma_latency)^2

  # Transform uncertainty (machine precision)
  sigma2_transform <- params$sigma_transform^2

  # Total (assuming independence)
  sigma2_total <- sigma2_tracker + sigma2_drift + sigma2_induced + sigma2_transform

  list(
    sigma2_tracker = sigma2_tracker,
    sigma2_drift = sigma2_drift,
    sigma2_induced = sigma2_induced,
    sigma2_transform = sigma2_transform,
    sigma2_total = sigma2_total,
    sigma_total = sqrt(sigma2_total)
  )
}

#' Generate boundary condition data (Prediction 4)
#'
#' @description
#' Creates gaze points specifically at slice transition boundaries
#' for testing probabilistic slice assignment.
#'
#' @param params Simulation parameters
#' @param n_boundary_points Number of points per boundary
#' @param sigma_delta Temporal uncertainty for boundary definition
#'
#' @return Data frame with gaze points at boundaries and expected slice probabilities
#' @export
generate_boundary_test_data <- function(params = default_sim_params(),
                                        n_boundary_points = 100,
                                        sigma_delta = 0.020) {
  # Generate z-axis events
  z_events <- simulate_z_axis_events(params, seed = 42)

  # Get transition times
  transitions <- z_events$event_time_sec[-1]  # Skip first event
  n_transitions <- length(transitions)

  if (n_transitions == 0) {
    warning("No transitions found in z-axis data")
    return(NULL)
  }

  # Generate gaze points around each transition
  results <- list()

  for (i in seq_along(transitions)) {
    t_transition <- transitions[i]

    # Points within ±3σ of transition
    offsets <- seq(-3 * sigma_delta, 3 * sigma_delta, length.out = n_boundary_points)
    times <- t_transition + offsets

    # Calculate probability of correct slice assignment
    # P(correct) = P(gaze before transition | slice before) + P(gaze after | slice after)
    # Using cumulative normal distribution
    p_before_transition <- pnorm(0, mean = offsets, sd = sigma_delta)

    # Slice before and after
    slice_before <- z_events$true_slice[i]
    slice_after <- z_events$true_slice[i + 1]

    results[[i]] <- data.frame(
      transition_id = i,
      time_sec = times,
      offset_from_transition = offsets,
      offset_sigma_units = offsets / sigma_delta,
      p_slice_before = p_before_transition,
      p_slice_after = 1 - p_before_transition,
      slice_before = slice_before,
      slice_after = slice_after,
      is_ambiguous = abs(offsets) <= sigma_delta,  # Within ±1σ
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}

#' Generate calibration drift test data (Prediction 5)
#'
#' @description
#' Creates long-session data for measuring calibration drift rate.
#'
#' @param session_duration Total session duration (seconds)
#' @param n_recalibrations Number of recalibration points
#' @param true_drift_rate Known drift rate for validation
#'
#' @return Data frame with calibration checkpoints
#' @export
generate_drift_test_data <- function(session_duration = 600,  # 10 minutes
                                     n_recalibrations = 6,
                                     true_drift_rate = 0.001) {
  # Calibration times
  calib_times <- seq(0, session_duration, length.out = n_recalibrations)

  # True drift at each checkpoint
  true_drift_x <- cumsum(c(0, rep(true_drift_rate * diff(calib_times)[1] / sqrt(2), n_recalibrations - 1)))
  true_drift_y <- cumsum(c(0, rep(true_drift_rate * diff(calib_times)[1] / sqrt(2), n_recalibrations - 1)))

  # Add measurement noise
  measured_drift_x <- true_drift_x + rnorm(n_recalibrations, 0, 0.002)
  measured_drift_y <- true_drift_y + rnorm(n_recalibrations, 0, 0.002)

  data.frame(
    calibration_time_sec = calib_times,
    true_drift_x = true_drift_x,
    true_drift_y = true_drift_y,
    measured_drift_x = measured_drift_x,
    measured_drift_y = measured_drift_y,
    total_drift = sqrt(measured_drift_x^2 + measured_drift_y^2)
  )
}

# =============================================================================
# VALIDATION METRICS
# =============================================================================

#' Validate velocity-spatial coupling (Prediction 1)
#'
#' @description
#' Tests whether σ_induced ≈ ||v_gaze|| · σ_δ holds in simulated data.
#'
#' @param gaze_data Gaze data with velocity and induced error
#' @param sigma_latency Known latency uncertainty
#'
#' @return List with validation results including r², slope, and expected slope
#' @export
validate_velocity_coupling <- function(gaze_data, sigma_latency = 0.020) {
  # Filter out zero velocity points
  moving <- gaze_data$velocity > 0

  if (sum(moving) < 10) {
    return(list(
      valid = FALSE,
      message = "Insufficient moving samples"
    ))
  }

  # Convert velocity to normalized units
  v_normalized <- gaze_data$velocity[moving] / 60

  # Get induced errors
  errors <- gaze_data$induced_error_magnitude[moving]

  # Fit linear model
  fit <- lm(errors ~ v_normalized)

  # Expected slope is sigma_latency
  expected_slope <- sigma_latency
  observed_slope <- coef(fit)[2]

  list(
    valid = TRUE,
    r_squared = summary(fit)$r.squared,
    observed_slope = observed_slope,
    expected_slope = expected_slope,
    slope_ratio = observed_slope / expected_slope,
    slope_within_tolerance = abs(observed_slope - expected_slope) / expected_slope < 0.2,
    n_samples = sum(moving),
    fit = fit
  )
}

#' Validate error independence (Prediction 2)
#'
#' @description
#' Tests whether error sources are weakly correlated (|ρ| < 0.3).
#'
#' @param gaze_data Gaze data with error components
#'
#' @return List with correlation matrix and independence assessment
#' @export
validate_error_independence <- function(gaze_data) {
  # Compute individual error components
  tracker_error_x <- gaze_data$gaze_point_on_display_area_x - gaze_data$true_x
  tracker_error_y <- gaze_data$gaze_point_on_display_area_y - gaze_data$true_y

  induced_error <- gaze_data$induced_error_magnitude

  # Correlation matrix
  errors <- data.frame(
    tracker = sqrt(tracker_error_x^2 + tracker_error_y^2),
    induced = induced_error
  )

  # Remove NA values
  errors <- na.omit(errors)

  cor_matrix <- cor(errors)

  # Check weak correlation criterion
  off_diagonal <- abs(cor_matrix[1, 2])

  list(
    correlation_matrix = cor_matrix,
    max_off_diagonal = off_diagonal,
    independence_satisfied = off_diagonal < 0.3,
    n_samples = nrow(errors)
  )
}

#' Validate temporal accumulation (Prediction 3)
#'
#' @description
#' Tests whether σ_δ(t) ≈ √(σ²_δ,0 + (α_clock·t)²).
#'
#' @param z_axis_data Z-axis data with accumulated uncertainty
#' @param sigma_delta_0 Initial uncertainty
#' @param alpha_clock Clock drift rate
#'
#' @return List with fit statistics and parameter estimates
#' @export
validate_temporal_accumulation <- function(z_axis_data,
                                           sigma_delta_0 = 0.001,
                                           alpha_clock = 1e-5) {
  t <- z_axis_data$event_time_sec
  sigma_observed <- z_axis_data$accumulated_sigma

  # Expected values
  sigma_expected <- sqrt(sigma_delta_0^2 + (alpha_clock * t)^2)

  # Fit the model to observed
  # Using squared form for linear fitting
  sigma2_observed <- sigma_observed^2
  t2 <- t^2

  fit <- lm(sigma2_observed ~ t2)

  # Estimated parameters
  sigma2_0_est <- coef(fit)[1]
  alpha2_est <- coef(fit)[2]

  list(
    r_squared = summary(fit)$r.squared,
    sigma_delta_0_estimated = sqrt(abs(sigma2_0_est)),
    alpha_clock_estimated = sqrt(abs(alpha2_est)),
    sigma_delta_0_true = sigma_delta_0,
    alpha_clock_true = alpha_clock,
    parameter_recovery_good = abs(sqrt(abs(alpha2_est)) - alpha_clock) / alpha_clock < 0.5,
    mean_residual = mean(abs(sigma_observed - sigma_expected))
  )
}

#' Validate boundary slice assignment (Prediction 4)
#'
#' @description
#' Tests probabilistic slice assignment at boundaries.
#'
#' @param boundary_data Data from generate_boundary_test_data()
#'
#' @return List with ambiguous region statistics
#' @export
validate_boundary_assignment <- function(boundary_data) {
  # Calculate fraction of points in ambiguous region
  frac_ambiguous <- mean(boundary_data$is_ambiguous)

  # Expected: ~68% within ±1σ (by construction)
  # But we're looking at ±3σ total range, so ±1σ is 1/3 of range
  expected_ambiguous <- 2 * pnorm(1) - 1  # ~0.68

  # Check probability calibration
  # For points at exact boundary (offset = 0), p should be 0.5
  near_boundary <- abs(boundary_data$offset_sigma_units) < 0.1
  p_at_boundary <- mean(boundary_data$p_slice_before[near_boundary])

  list(
    fraction_ambiguous = frac_ambiguous,
    expected_fraction_within_1sigma = expected_ambiguous,
    mean_p_at_boundary = p_at_boundary,
    p_at_boundary_calibrated = abs(p_at_boundary - 0.5) < 0.1,
    n_transitions = length(unique(boundary_data$transition_id)),
    n_points = nrow(boundary_data)
  )
}

#' Validate calibration drift (Prediction 5)
#'
#' @description
#' Estimates drift rate from calibration checkpoint data.
#'
#' @param drift_data Data from generate_drift_test_data()
#'
#' @return List with estimated drift rate and recalibration recommendations
#' @export
validate_calibration_drift <- function(drift_data) {
  t <- drift_data$calibration_time_sec
  drift <- drift_data$total_drift

  # Fit linear model to estimate drift rate
  fit <- lm(drift ~ t)

  estimated_rate <- coef(fit)[2]

  # Recalibration interval recommendation

  # If acceptable drift is 0.01 (1% of screen), time = 0.01 / rate
  acceptable_drift <- 0.01
  recommended_interval <- acceptable_drift / estimated_rate

  list(
    estimated_drift_rate = estimated_rate,
    drift_rate_per_minute = estimated_rate * 60,
    r_squared = summary(fit)$r.squared,
    recommended_recalibration_interval_sec = recommended_interval,
    recommended_recalibration_interval_min = recommended_interval / 60
  )
}

# =============================================================================
# COMPLETE SIMULATION WORKFLOW
# =============================================================================

#' Run complete simulation with all uncertainty sources
#'
#' @description
#' Generates a complete simulated dataset for validating the gazeneuro framework.
#'
#' @param params Simulation parameters
#' @param latency_ms True latency between systems
#' @param seed Random seed
#'
#' @return List containing gaze_data, z_axis_data, ground_truth, and validation_results
#' @export
run_complete_simulation <- function(params = default_sim_params(),
                                    latency_ms = 50,
                                    seed = 123) {
  set.seed(seed)

  message("Generating gaze stream...")
  gaze_data <- simulate_gaze_stream(params, seed = seed)

  message("Generating z-axis events...")
  z_axis_data <- simulate_z_axis_events(params, gaze_data, latency_ms, seed = seed + 1)

  message("Injecting velocity-induced errors...")
  gaze_data <- inject_velocity_induced_error(gaze_data, params$sigma_latency)

  message("Injecting clock drift...")
  z_axis_data <- inject_clock_drift(z_axis_data, alpha_clock = params$alpha_clock)

  message("Generating boundary test data...")
  boundary_data <- generate_boundary_test_data(params, sigma_delta = params$sigma_latency)

  message("Generating drift test data...")
  drift_data <- generate_drift_test_data(
    session_duration = params$session_duration_sec * 10,
    true_drift_rate = params$sigma_drift_rate
  )

  # Store ground truth
  ground_truth <- list(
    true_latency_ms = latency_ms,
    params = params
  )

  # Run validations
  message("\nRunning validations...")

  validation_results <- list(
    P1_velocity_coupling = validate_velocity_coupling(gaze_data, params$sigma_latency),
    P2_independence = validate_error_independence(gaze_data),
    P3_temporal_accumulation = validate_temporal_accumulation(z_axis_data,
                                                              alpha_clock = params$alpha_clock),
    P4_boundary_assignment = validate_boundary_assignment(boundary_data),
    P5_calibration_drift = validate_calibration_drift(drift_data)
  )

  # Print summary
  message("\n=== VALIDATION SUMMARY ===")
  message("P1 (Velocity-Spatial Coupling): ",
          ifelse(validation_results$P1_velocity_coupling$slope_within_tolerance, "✓ PASS", "✗ FAIL"),
          sprintf(" (R² = %.3f)", validation_results$P1_velocity_coupling$r_squared))
  message("P2 (Error Independence): ",
          ifelse(validation_results$P2_independence$independence_satisfied, "✓ PASS", "✗ FAIL"),
          sprintf(" (max |ρ| = %.3f)", validation_results$P2_independence$max_off_diagonal))
  message("P3 (Temporal Accumulation): ",
          ifelse(validation_results$P3_temporal_accumulation$parameter_recovery_good, "✓ PASS", "✗ FAIL"),
          sprintf(" (R² = %.3f)", validation_results$P3_temporal_accumulation$r_squared))
  message("P4 (Boundary Assignment): ",
          ifelse(validation_results$P4_boundary_assignment$p_at_boundary_calibrated, "✓ PASS", "✗ FAIL"),
          sprintf(" (P at boundary = %.3f)", validation_results$P4_boundary_assignment$mean_p_at_boundary))
  message("P5 (Calibration Drift): ",
          sprintf("Rate = %.4f/min, Recalibrate every %.1f min",
                  validation_results$P5_calibration_drift$drift_rate_per_minute,
                  validation_results$P5_calibration_drift$recommended_recalibration_interval_min))

  list(
    gaze_data = gaze_data,
    z_axis_data = z_axis_data,
    boundary_data = boundary_data,
    drift_data = drift_data,
    ground_truth = ground_truth,
    validation_results = validation_results
  )
}
