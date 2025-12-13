#' Visualization of Simulation Framework Predictions
#'
#' Generates diagnostic plots demonstrating the five theoretical predictions
#' from the gazeneuro uncertainty framework.

library(ggplot2)
library(gridExtra)
library(dplyr)

# Source simulation functions (when running standalone)
# source("simulation.R")

#' Create comprehensive visualization of all five predictions
#'
#' @param output_file Path to save the visualization (NULL for display only
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @export
visualize_uncertainty_predictions <- function(output_file = NULL,
                                              width = 14,
                                              height = 10) {

  # Set seed for reproducibility
  set.seed(42)

  # Get default parameters
  params <- default_sim_params()
  params$session_duration_sec <- 120  # 2 minutes for better visualization

  message("Generating simulation data...")

  # Generate base data
  gaze_data <- simulate_gaze_stream(params, seed = 42)
  z_axis_data <- simulate_z_axis_events(params, gaze_data, latency_ms = 50, seed = 43)

  # Inject errors
  gaze_data <- inject_velocity_induced_error(gaze_data, sigma_latency = params$sigma_latency)
  z_axis_data <- inject_clock_drift(z_axis_data, sigma_delta_0 = 0.001, alpha_clock = 1e-5)

  # Generate specialized test data
  boundary_data <- generate_boundary_test_data(params, n_boundary_points = 100,
                                               sigma_delta = params$sigma_latency)
  drift_data <- generate_drift_test_data(session_duration = 600,
                                         n_recalibrations = 12,
                                         true_drift_rate = 0.001)

  message("Creating visualizations...")

  # =========================================================================
  # P1: Velocity-Spatial Coupling
  # =========================================================================

  # Compute binned statistics
  moving <- gaze_data$v_normalized > 0.01
  gaze_moving <- gaze_data[moving, ]

  # Create velocity bins
  gaze_moving$v_bin <- cut(gaze_moving$v_normalized,
                           breaks = quantile(gaze_moving$v_normalized,
                                             probs = seq(0, 1, 0.1)),
                           include.lowest = TRUE)

  bin_stats <- gaze_moving %>%
    group_by(v_bin) %>%
    summarise(
      v_mean = mean(v_normalized),
      error_mean = mean(induced_error_magnitude),
      error_se = sd(induced_error_magnitude) / sqrt(n()),
      n = n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(v_bin))

  # Theoretical line
  v_range <- seq(0, max(bin_stats$v_mean, na.rm = TRUE), length.out = 100)
  theoretical_error <- v_range * params$sigma_latency * sqrt(2/pi)
  theory_df <- data.frame(v = v_range, error = theoretical_error)

  p1 <- ggplot() +
    geom_ribbon(data = bin_stats,
                aes(x = v_mean,
                    ymin = error_mean - 1.96 * error_se,
                    ymax = error_mean + 1.96 * error_se),
                fill = "steelblue", alpha = 0.3) +
    geom_point(data = bin_stats,
               aes(x = v_mean, y = error_mean, size = n),
               color = "steelblue") +
    geom_line(data = theory_df,
              aes(x = v, y = error),
              color = "red", linetype = "dashed", linewidth = 1) +
    labs(
      title = "P1: Velocity-Spatial Coupling",
      subtitle = expression(sigma[induced] %~~% "||v|| ·" ~ sigma[delta]),
      x = "Normalized Velocity (screen widths/sec)",
      y = "Mean Induced Error (normalized)",
      size = "n samples"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      legend.position = c(0.15, 0.85)
    ) +
    annotate("text", x = max(bin_stats$v_mean) * 0.7,
             y = max(bin_stats$error_mean) * 0.3,
             label = "— Theoretical\n• Observed (±95% CI)",
             hjust = 0, size = 3, color = "gray30")

  # =========================================================================
  # P2: Error Independence
  # =========================================================================

  # Get error components for moving samples
  gaze_moving$tracker_mag <- sqrt(gaze_moving$tracker_noise_x^2 +
                                    gaze_moving$tracker_noise_y^2)

  # Sample for plotting (too many points otherwise)
  set.seed(42)
  sample_idx <- sample(nrow(gaze_moving), min(500, nrow(gaze_moving)))
  plot_data <- gaze_moving[sample_idx, ]

  # Compute correlation
  cor_val <- cor(plot_data$tracker_noise_x, plot_data$induced_error_x)

  p2 <- ggplot(plot_data, aes(x = tracker_noise_x, y = induced_error_x)) +
    geom_point(alpha = 0.3, color = "steelblue", size = 1) +
    geom_smooth(method = "lm", color = "red", se = TRUE, linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
    labs(
      title = "P2: Error Independence",
      subtitle = expression("Tracker vs Induced Error (" * rho < 0.3 * " required)"),
      x = "Tracker Error (X component)",
      y = "Induced Error (X component)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, color = "gray40")
    ) +
    annotate("label", x = min(plot_data$tracker_noise_x) * 0.8,
             y = max(plot_data$induced_error_x) * 0.8,
             label = sprintf("ρ = %.3f\n%s", cor_val,
                             ifelse(abs(cor_val) < 0.3, "✓ Independent", "✗ Correlated")),
             size = 3.5, fill = ifelse(abs(cor_val) < 0.3, "#d4edda", "#f8d7da"))

  # =========================================================================
  # P3: Temporal Uncertainty Accumulation
  # =========================================================================

  # Extended time range for visualization
  t_range <- seq(0, max(z_axis_data$event_time_sec), length.out = 100)
  sigma_delta_0 <- 0.001
  alpha_clock <- 1e-5

  # Theoretical curve
  sigma_theoretical <- sqrt(sigma_delta_0^2 + (alpha_clock * t_range)^2)
  theory_p3 <- data.frame(t = t_range, sigma = sigma_theoretical * 1000)  # Convert to ms

  # Observed from z-axis data
  z_axis_data$sigma_ms <- z_axis_data$accumulated_sigma * 1000

  p3 <- ggplot() +
    geom_line(data = theory_p3, aes(x = t, y = sigma),
              color = "red", linetype = "dashed", linewidth = 1) +
    geom_point(data = z_axis_data,
               aes(x = event_time_sec, y = sigma_ms),
               color = "steelblue", alpha = 0.6, size = 2) +
    labs(
      title = "P3: Temporal Uncertainty Accumulation",
      subtitle = expression(sigma[delta](t) %~~% sqrt(sigma[0]^2 + (alpha * t)^2)),
      x = "Session Time (seconds)",
      y = expression(sigma[delta] ~ "(ms)")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, color = "gray40")
    ) +
    annotate("text", x = max(z_axis_data$event_time_sec) * 0.6,
             y = max(z_axis_data$sigma_ms) * 0.3,
             label = sprintf("α = %.0e (clock drift rate)\nσ₀ = %.1f ms (initial)",
                             alpha_clock, sigma_delta_0 * 1000),
             hjust = 0, size = 3, color = "gray30")

  # =========================================================================
  # P4: Probabilistic Slice Assignment
  # =========================================================================

  if (!is.null(boundary_data) && nrow(boundary_data) > 0) {
    # Take first few transitions for clarity
    boundary_subset <- boundary_data %>%
      filter(transition_id <= 3)

    p4 <- ggplot(boundary_subset, aes(x = offset_sigma_units, y = p_slice_before)) +
      geom_line(aes(group = transition_id, color = factor(transition_id)),
                linewidth = 1) +
      geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray50") +
      geom_vline(xintercept = 0, linetype = "solid", color = "black") +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
      geom_rect(aes(xmin = -1, xmax = 1, ymin = 0, ymax = 1),
                fill = "yellow", alpha = 0.1) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
      labs(
        title = "P4: Probabilistic Slice Assignment",
        subtitle = "Ambiguous zone within ±1σ of transition",
        x = expression("Offset from transition (" * sigma * " units)"),
        y = "P(assigned to slice before transition)",
        color = "Transition"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        legend.position = "none"
      ) +
      annotate("text", x = 0, y = 0.15, label = "Ambiguous\nZone",
               size = 3, color = "orange3", fontface = "bold")
  } else {
    p4 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No boundary data") +
      theme_void()
  }

  # =========================================================================
  # P5: Calibration Drift
  # =========================================================================

  # Fit drift model
  fit <- lm(total_drift ~ calibration_time_sec, data = drift_data)
  drift_rate <- coef(fit)[2]

  # Acceptable drift threshold
  acceptable_drift <- 0.01  # 1% of screen
  recalib_time <- acceptable_drift / drift_rate

  p5 <- ggplot(drift_data, aes(x = calibration_time_sec / 60, y = total_drift * 100)) +
    geom_point(color = "steelblue", size = 3) +
    geom_smooth(method = "lm", color = "red", se = TRUE, linewidth = 1) +
    geom_hline(yintercept = acceptable_drift * 100,
               linetype = "dashed", color = "orange", linewidth = 1) +
    geom_vline(xintercept = recalib_time / 60,
               linetype = "dotted", color = "orange", linewidth = 1) +
    labs(
      title = "P5: Calibration Drift",
      subtitle = "Determines recalibration schedule",
      x = "Session Time (minutes)",
      y = "Total Drift (% of screen)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, color = "gray40")
    ) +
    annotate("label", x = max(drift_data$calibration_time_sec) / 60 * 0.6,
             y = acceptable_drift * 100 * 1.3,
             label = sprintf("Drift rate: %.3f%%/min\nRecalibrate every %.1f min",
                             drift_rate * 60 * 100, recalib_time / 60),
             size = 3, fill = "#fff3cd")

  # =========================================================================
  # Summary Panel: Total Uncertainty Composition
  # =========================================================================

  # Calculate uncertainty components over time for different velocities
  times <- seq(0, 300, by = 10)
  velocities <- c(0, 100, 300)  # deg/s

  uncertainty_data <- expand.grid(t = times, v = velocities) %>%
    mutate(
      sigma2_tracker = params$sigma_tracker^2,
      sigma2_drift = (params$sigma_drift_rate * t)^2,
      sigma2_induced = ((v / 60) * params$sigma_latency)^2,
      sigma2_transform = params$sigma_transform^2,
      sigma_total = sqrt(sigma2_tracker + sigma2_drift + sigma2_induced + sigma2_transform),
      velocity_label = factor(paste0(v, "°/s"),
                              levels = c("0°/s", "100°/s", "300°/s"))
    )

  p6 <- ggplot(uncertainty_data, aes(x = t, y = sigma_total * 100,
                                     color = velocity_label,
                                     linetype = velocity_label)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = c("0°/s" = "forestgreen",
                                  "100°/s" = "steelblue",
                                  "300°/s" = "firebrick")) +
    labs(
      title = "Total Uncertainty: σ²_total = σ²_tracker + σ²_drift + σ²_induced + σ²_transform",
      subtitle = "Uncertainty grows with time (drift) and velocity (induced error)",
      x = "Session Time (seconds)",
      y = "Total Spatial Uncertainty (% of screen)",
      color = "Gaze Velocity",
      linetype = "Gaze Velocity"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      legend.position = c(0.15, 0.8)
    )

  # =========================================================================
  # Combine all plots
  # =========================================================================

  combined <- grid.arrange(
    p1, p2, p3, p4, p5, p6,
    ncol = 3, nrow = 2,
    top = grid::textGrob("Gazeneuro Simulation Framework: Uncertainty Prediction Validation",
                         gp = grid::gpar(fontsize = 16, fontface = "bold"))
  )

  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, combined, width = width, height = height, dpi = 150)
    message("Saved visualization to: ", output_file)
  }

  invisible(combined)
}

# Run if executed directly
if (!interactive()) {
  visualize_uncertainty_predictions("uncertainty_predictions.png")
}
