# Tests for simulation.R
# Tests the simulation framework for gazeneuro validation

test_that("default_sim_params returns valid parameter list", {
  params <- default_sim_params()

  # Check required parameters exist
  expect_true("sampling_rate_hz" %in% names(params))
  expect_true("session_duration_sec" %in% names(params))
  expect_true("sigma_tracker" %in% names(params))
  expect_true("sigma_drift_rate" %in% names(params))
  expect_true("sigma_latency" %in% names(params))
  expect_true("alpha_clock" %in% names(params))

  # Check reasonable values

  expect_gte(params$sampling_rate_hz, 10)
  expect_lte(params$sampling_rate_hz, 2000)
  expect_gt(params$session_duration_sec, 0)
  expect_gt(params$sigma_tracker, 0)
  expect_lt(params$sigma_tracker, 0.1)  # Less than 10% of screen
  expect_gt(params$n_slices, 0)
})

# =============================================================================
# GAZE STREAM TESTS
# =============================================================================

test_that("simulate_gaze_stream generates correct number of samples", {
  params <- default_sim_params()
  params$sampling_rate_hz <- 120
  params$session_duration_sec <- 10

  gaze_data <- simulate_gaze_stream(params, seed = 123)

  expected_samples <- params$sampling_rate_hz * params$session_duration_sec
  expect_equal(nrow(gaze_data), expected_samples)
})

test_that("simulate_gaze_stream produces valid coordinate ranges", {
  params <- default_sim_params()
  gaze_data <- simulate_gaze_stream(params, seed = 123)

  # Observed coordinates should be in [0, 1]
  expect_true(all(gaze_data$gaze_point_on_display_area_x >= 0))
  expect_true(all(gaze_data$gaze_point_on_display_area_x <= 1))
  expect_true(all(gaze_data$gaze_point_on_display_area_y >= 0))
  expect_true(all(gaze_data$gaze_point_on_display_area_y <= 1))
})

test_that("simulate_gaze_stream contains both fixations and saccades", {
  params <- default_sim_params()
  gaze_data <- simulate_gaze_stream(params, seed = 123)

  # Should have both fixations and saccades
  expect_true(any(gaze_data$is_fixation == TRUE))
  expect_true(any(gaze_data$is_fixation == FALSE))

  # Fixations should have velocity near 0
  fixation_velocities <- gaze_data$velocity[gaze_data$is_fixation]
  expect_true(all(fixation_velocities == 0))

  # Saccades should have positive velocity
  saccade_velocities <- gaze_data$velocity[!gaze_data$is_fixation]
  expect_true(all(saccade_velocities > 0))
})

test_that("simulate_gaze_stream respects noise parameters", {
  params <- default_sim_params()

  # Generate without noise
  gaze_no_noise <- simulate_gaze_stream(params, add_tracker_noise = FALSE,
                                        add_drift = FALSE, seed = 123)

  # Generate with noise
  gaze_with_noise <- simulate_gaze_stream(params, add_tracker_noise = TRUE,
                                          add_drift = FALSE, seed = 123)

  # Without noise, observed should equal true (approximately, due to clipping)
  # For fixations only (no movement)
  fix_idx <- gaze_no_noise$is_fixation &
    gaze_no_noise$true_x > 0.1 & gaze_no_noise$true_x < 0.9
  expect_equal(
    gaze_no_noise$gaze_point_on_display_area_x[fix_idx][1:10],
    gaze_no_noise$true_x[fix_idx][1:10],
    tolerance = 1e-10
  )

  # With noise, there should be differences
  diff_with_noise <- abs(gaze_with_noise$gaze_point_on_display_area_x - gaze_with_noise$true_x)
  expect_gt(mean(diff_with_noise), 0)
})

test_that("simulate_gaze_stream is reproducible with seed", {
  params <- default_sim_params()

  gaze1 <- simulate_gaze_stream(params, seed = 42)
  gaze2 <- simulate_gaze_stream(params, seed = 42)

  expect_equal(gaze1$gaze_point_on_display_area_x, gaze2$gaze_point_on_display_area_x)
  expect_equal(gaze1$gaze_point_on_display_area_y, gaze2$gaze_point_on_display_area_y)
})

# =============================================================================
# Z-AXIS EVENT TESTS
# =============================================================================

test_that("simulate_z_axis_events generates events", {
  params <- default_sim_params()
  z_axis <- simulate_z_axis_events(params, seed = 123)

  expect_gt(nrow(z_axis), 0)
  expect_true("client_timestamp" %in% names(z_axis))
  expect_true("plane" %in% names(z_axis))
  expect_true("index" %in% names(z_axis))
})

test_that("simulate_z_axis_events produces valid slice indices", {
  params <- default_sim_params()
  z_axis <- simulate_z_axis_events(params, seed = 123)

  # Normalized index should be in [0, 1]
  expect_true(all(z_axis$index >= 0))
  expect_true(all(z_axis$index <= 1))

  # True slice should be in valid range
  expect_true(all(z_axis$true_slice >= 1))
  expect_true(all(z_axis$true_slice <= params$n_slices))
})

test_that("simulate_z_axis_events applies latency correctly", {
  params <- default_sim_params()
  latency_ms <- 100

  z_axis <- simulate_z_axis_events(params, latency_ms = latency_ms, seed = 123)

  # Client timestamp should be offset from event time by latency
  expected_offset <- latency_ms / 1000
  actual_offset <- z_axis$client_timestamp / 1000 - z_axis$event_time_sec

  expect_equal(mean(actual_offset), expected_offset, tolerance = 0.001)
})

# =============================================================================
# UNCERTAINTY INJECTION TESTS
# =============================================================================

test_that("inject_velocity_induced_error creates velocity-dependent errors", {
  params <- default_sim_params()
  gaze_data <- simulate_gaze_stream(params, add_tracker_noise = FALSE,
                                    add_drift = FALSE, seed = 123)

  gaze_with_error <- inject_velocity_induced_error(gaze_data, sigma_latency = 0.020)

  # Induced error should exist
  expect_true("induced_error_magnitude" %in% names(gaze_with_error))

  # Higher velocity should correlate with larger errors
  moving <- gaze_with_error$velocity > 0
  if (sum(moving) > 10) {
    cor_test <- cor(gaze_with_error$velocity[moving],
                    gaze_with_error$induced_error_magnitude[moving])
    expect_gt(cor_test, 0)  # Positive correlation
  }
})

test_that("inject_clock_drift accumulates over time", {
  params <- default_sim_params()
  params$session_duration_sec <- 300  # Longer session
  z_axis <- simulate_z_axis_events(params, seed = 123)

  z_with_drift <- inject_clock_drift(z_axis, sigma_delta_0 = 0.001, alpha_clock = 1e-4)

  # Accumulated sigma should increase over time
  expect_true("accumulated_sigma" %in% names(z_with_drift))

  # Later events should have higher accumulated uncertainty
  early_sigma <- mean(z_with_drift$accumulated_sigma[1:5])
  late_sigma <- mean(tail(z_with_drift$accumulated_sigma, 5))
  expect_gt(late_sigma, early_sigma)
})

# =============================================================================
# PREDICTION VALIDATION TESTS
# =============================================================================

test_that("validate_velocity_coupling returns expected structure", {
  params <- default_sim_params()
  gaze_data <- simulate_gaze_stream(params, seed = 123)
  gaze_data <- inject_velocity_induced_error(gaze_data, sigma_latency = params$sigma_latency)

  result <- validate_velocity_coupling(gaze_data, sigma_latency = params$sigma_latency)

  expect_true("valid" %in% names(result))
  expect_true("r_squared" %in% names(result))
  expect_true("observed_slope" %in% names(result))
  expect_true("expected_slope" %in% names(result))
})

test_that("validate_velocity_coupling recovers slope approximately", {
  params <- default_sim_params()
  sigma_latency <- 0.030

  gaze_data <- simulate_gaze_stream(params, add_tracker_noise = FALSE,
                                    add_drift = FALSE, seed = 123)
  gaze_data <- inject_velocity_induced_error(gaze_data, sigma_latency = sigma_latency)

  result <- validate_velocity_coupling(gaze_data, sigma_latency = sigma_latency)

  # Slope should be close to sigma_latency (within factor of 2)
  expect_gt(result$observed_slope, sigma_latency / 3)
  expect_lt(result$observed_slope, sigma_latency * 3)
})

test_that("validate_error_independence detects weak correlation", {
  params <- default_sim_params()
  gaze_data <- simulate_gaze_stream(params, seed = 123)
  gaze_data <- inject_velocity_induced_error(gaze_data, sigma_latency = params$sigma_latency)

  result <- validate_error_independence(gaze_data)

  expect_true("correlation_matrix" %in% names(result))
  expect_true("independence_satisfied" %in% names(result))

  # With independent error sources, correlation should be low
  expect_lt(result$max_off_diagonal, 0.5)
})

test_that("validate_temporal_accumulation recovers parameters", {
  params <- default_sim_params()
  z_axis <- simulate_z_axis_events(params, seed = 123)

  sigma_delta_0 <- 0.002
  alpha_clock <- 5e-5

  z_with_drift <- inject_clock_drift(z_axis, sigma_delta_0 = sigma_delta_0,
                                     alpha_clock = alpha_clock)

  result <- validate_temporal_accumulation(z_with_drift,
                                           sigma_delta_0 = sigma_delta_0,
                                           alpha_clock = alpha_clock)

  expect_true("r_squared" %in% names(result))
  expect_gt(result$r_squared, 0.5)  # Should fit well
})

test_that("validate_boundary_assignment produces calibrated probabilities", {
  params <- default_sim_params()
  boundary_data <- generate_boundary_test_data(params, n_boundary_points = 50,
                                               sigma_delta = 0.020)

  if (!is.null(boundary_data)) {
    result <- validate_boundary_assignment(boundary_data)

    expect_true("mean_p_at_boundary" %in% names(result))
    # Probability at exact boundary should be close to 0.5
    expect_gt(result$mean_p_at_boundary, 0.3)
    expect_lt(result$mean_p_at_boundary, 0.7)
  }
})

test_that("validate_calibration_drift estimates rate", {
  drift_data <- generate_drift_test_data(session_duration = 600,
                                         n_recalibrations = 10,
                                         true_drift_rate = 0.002)

  result <- validate_calibration_drift(drift_data)

  expect_true("estimated_drift_rate" %in% names(result))
  expect_true("recommended_recalibration_interval_sec" %in% names(result))

  # Estimated rate should be in reasonable range
  expect_gt(result$estimated_drift_rate, 0)
  expect_lt(result$estimated_drift_rate, 0.1)
})

# =============================================================================
# EXPECTED UNCERTAINTY CALCULATION TESTS
# =============================================================================

test_that("calculate_expected_uncertainty returns correct structure", {
  result <- calculate_expected_uncertainty()

  expect_true("sigma2_tracker" %in% names(result))
  expect_true("sigma2_drift" %in% names(result))
  expect_true("sigma2_induced" %in% names(result))
  expect_true("sigma2_transform" %in% names(result))
  expect_true("sigma2_total" %in% names(result))
  expect_true("sigma_total" %in% names(result))
})

test_that("calculate_expected_uncertainty composes correctly", {
  params <- default_sim_params()
  result <- calculate_expected_uncertainty(params, t = 30, velocity = 0)

  # Total should be sum of components (for independent sources)
  expected_total <- result$sigma2_tracker + result$sigma2_drift +
    result$sigma2_induced + result$sigma2_transform

  expect_equal(result$sigma2_total, expected_total, tolerance = 1e-20)
})

test_that("calculate_expected_uncertainty increases with velocity", {
  params <- default_sim_params()

  result_static <- calculate_expected_uncertainty(params, t = 30, velocity = 0)
  result_moving <- calculate_expected_uncertainty(params, t = 30, velocity = 300)

  # Moving should have higher total uncertainty
  expect_gt(result_moving$sigma2_total, result_static$sigma2_total)

  # The difference should be in the induced component
  expect_gt(result_moving$sigma2_induced, result_static$sigma2_induced)
})

test_that("calculate_expected_uncertainty increases with time due to drift", {
  params <- default_sim_params()

  result_early <- calculate_expected_uncertainty(params, t = 10, velocity = 0)
  result_late <- calculate_expected_uncertainty(params, t = 100, velocity = 0)

  # Later should have higher drift uncertainty
  expect_gt(result_late$sigma2_drift, result_early$sigma2_drift)
})

# =============================================================================
# COMPLETE SIMULATION TESTS
# =============================================================================

test_that("run_complete_simulation produces all expected outputs", {
  params <- default_sim_params()
  params$session_duration_sec <- 30  # Shorter for testing

  result <- run_complete_simulation(params, latency_ms = 50, seed = 123)

  expect_true("gaze_data" %in% names(result))
  expect_true("z_axis_data" %in% names(result))
  expect_true("boundary_data" %in% names(result))
  expect_true("drift_data" %in% names(result))
  expect_true("ground_truth" %in% names(result))
  expect_true("validation_results" %in% names(result))

  # Check validation results structure
  expect_true("P1_velocity_coupling" %in% names(result$validation_results))
  expect_true("P2_independence" %in% names(result$validation_results))
  expect_true("P3_temporal_accumulation" %in% names(result$validation_results))
  expect_true("P4_boundary_assignment" %in% names(result$validation_results))
  expect_true("P5_calibration_drift" %in% names(result$validation_results))
})

test_that("run_complete_simulation is reproducible", {
  params <- default_sim_params()
  params$session_duration_sec <- 10

  result1 <- run_complete_simulation(params, seed = 42)
  result2 <- run_complete_simulation(params, seed = 42)

  expect_equal(result1$gaze_data$gaze_point_on_display_area_x,
               result2$gaze_data$gaze_point_on_display_area_x)
})

# =============================================================================
# INTEGRATION WITH GAZENEURO PACKAGE TESTS
# =============================================================================

test_that("simulated data format matches gazeneuro expectations", {
  params <- default_sim_params()
  gaze_data <- simulate_gaze_stream(params, seed = 123)

  # Check column names match what integrate_all_gaze_points expects
  expect_true("device_time_stamp" %in% names(gaze_data))
  expect_true("gaze_point_on_display_area_x" %in% names(gaze_data))
  expect_true("gaze_point_on_display_area_y" %in% names(gaze_data))
})

test_that("simulated z_axis format matches gazeneuro expectations", {
  params <- default_sim_params()
  z_axis <- simulate_z_axis_events(params, seed = 123)

  # Check column names match what integrate_all_gaze_points expects
  expect_true("client_timestamp" %in% names(z_axis))
  expect_true("plane" %in% names(z_axis))
  expect_true("index" %in% names(z_axis))
  expect_true("image_id" %in% names(z_axis))
})

# Note: Integration with actual gazeneuro functions would require
# the package to be loaded. These tests verify format compatibility.

# =============================================================================
# EDGE CASE TESTS
# =============================================================================

test_that("simulation handles very short sessions", {
  params <- default_sim_params()
  params$session_duration_sec <- 1  # Very short

  expect_no_error({
    gaze_data <- simulate_gaze_stream(params, seed = 123)
    z_axis <- simulate_z_axis_events(params, seed = 123)
  })

  expect_gt(nrow(gaze_data), 0)
})

test_that("simulation handles very low sampling rate", {
  params <- default_sim_params()
  params$sampling_rate_hz <- 10  # Very low
  params$session_duration_sec <- 10

  gaze_data <- simulate_gaze_stream(params, seed = 123)

  expect_equal(nrow(gaze_data), 100)  # 10 Hz × 10 sec
})

test_that("simulation handles very high sampling rate", {
  params <- default_sim_params()
  params$sampling_rate_hz <- 1000  # Very high
  params$session_duration_sec <- 5

  gaze_data <- simulate_gaze_stream(params, seed = 123)

  expect_equal(nrow(gaze_data), 5000)  # 1000 Hz × 5 sec
})

test_that("simulation handles extreme latency values", {
  params <- default_sim_params()

  # Zero latency
  z_zero <- simulate_z_axis_events(params, latency_ms = 0, seed = 123)
  expect_gt(nrow(z_zero), 0)

  # Large latency
  z_large <- simulate_z_axis_events(params, latency_ms = 1000, seed = 123)
  expect_gt(nrow(z_large), 0)

  # Timestamps should differ by ~1000ms
  expect_gt(mean(z_large$client_timestamp - z_zero$client_timestamp), 900)
})
