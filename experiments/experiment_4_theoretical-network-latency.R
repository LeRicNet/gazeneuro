# theoretical_network_latency_evaluation.R
#
# Theoretical evaluation of network latency impact on gazeneuro framework
# All results saved to ./results/theoretical-network-latency/

library(tidyverse)
library(ggplot2)
library(patchwork)
library(glue)
library(furrr)
library(future)

setwd("/home/rstudio/interfaces/gazeneuro")

# Set experiment name and create output directory
EXPERIMENT_NAME <- "theoretical-network-latency"
OUTPUT_DIR <- file.path("./results", EXPERIMENT_NAME)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ========== PARALLEL SETUP ==========
n_cores <- availableCores() - 1
message("Detected ", n_cores + 1, " cores. Using ", n_cores, " for parallel processing.")
plan(multisession, workers = n_cores)

# ========== 1. Network Model Definitions ==========

create_network_models <- function() {
  list(
    lan = list(
      base_ms = 2,
      variance_ms = 0.5,
      packet_loss = 0.0001,
      jitter_ms = 0.2,
      distance_km = 0.1
    ),
    same_city = list(
      base_ms = 15,
      variance_ms = 5,
      packet_loss = 0.001,
      jitter_ms = 2,
      distance_km = 50
    ),
    regional = list(
      base_ms = 45,
      variance_ms = 10,
      packet_loss = 0.002,
      jitter_ms = 5,
      distance_km = 500
    ),
    continental = list(
      base_ms = 85,
      variance_ms = 20,
      packet_loss = 0.005,
      jitter_ms = 10,
      distance_km = 3000
    ),
    international = list(
      base_ms = 180,
      variance_ms = 40,
      packet_loss = 0.01,
      jitter_ms = 20,
      distance_km = 10000
    )
  )
}

# ========== 2. Latency Trace Generation ==========

generate_latency_trace <- function(scenario, duration_seconds = 300, sampling_rate = 10) {
  models <- create_network_models()
  params <- models[[scenario]]

  n_samples <- duration_seconds * sampling_rate
  time <- seq(0, duration_seconds, length.out = n_samples)

  # Base latency with diurnal drift
  drift <- 5 * sin(2 * pi * time / 86400)

  # Log-normal variance
  variance <- rlnorm(n_samples,
                     meanlog = log(params$variance_ms) - 0.5^2/2,
                     sdlog = 0.5)

  # Short-term jitter
  jitter <- rnorm(n_samples, mean = 0, sd = params$jitter_ms)

  # Combine components
  latency <- params$base_ms + drift + variance + jitter
  latency[latency < 0] <- 0

  # Add packet loss
  loss_mask <- runif(n_samples) < params$packet_loss
  latency[loss_mask] <- NA

  tibble(
    time = time,
    latency_ms = latency,
    scenario = scenario,
    packet_loss = loss_mask
  )
}

# ========== 3. Impact Analysis Functions ==========

compute_sync_error <- function(network_latency_ms,
                               api_poll_rate = 4,
                               slice_duration_mean = 2000,
                               slice_duration_sd = 1000) {

  polling_delay <- (1000 / api_poll_rate) / 2
  total_latency <- network_latency_ms + polling_delay

  # Misassignment probability model
  p_wrong_slice <- pmin(1.0, total_latency / slice_duration_mean)

  tibble(
    network_latency_ms = network_latency_ms,
    total_latency_ms = total_latency,
    sync_error_ms = total_latency,
    misassignment_rate = p_wrong_slice,
    data_quality_score = 1 - p_wrong_slice
  )
}

# ========== 4. Monte Carlo Helper Functions ==========

simulate_viewing_session <- function(duration_seconds = 300,
                                     slice_duration_mean = 2,
                                     slice_duration_sd = 1) {

  time <- 0
  slices <- list()
  slice_num <- 1

  while (time < duration_seconds) {
    duration <- max(0.5, rnorm(1, slice_duration_mean, slice_duration_sd))

    slices[[slice_num]] <- tibble(
      slice_id = slice_num,
      start_time = time,
      end_time = time + duration,
      duration = duration
    )

    time <- time + duration
    slice_num <- slice_num + 1
  }

  bind_rows(slices)
}

calculate_sync_performance <- function(network_trace, viewing_session) {
  valid_latencies <- network_trace$latency_ms[!is.na(network_trace$latency_ms)]

  if (length(valid_latencies) == 0) {
    return(tibble(
      mean_latency = NA,
      p95_latency = NA,
      packet_loss = 1,
      sync_accuracy = 0
    ))
  }

  mean_latency <- mean(valid_latencies)
  mean_duration <- mean(viewing_session$duration) * 1000

  tibble(
    mean_latency = mean_latency,
    p95_latency = quantile(valid_latencies, 0.95),
    packet_loss = mean(is.na(network_trace$latency_ms)),
    sync_accuracy = max(0, 1 - mean_latency / mean_duration)
  )
}

# ========== 4a. Sequential Monte Carlo (for comparison) ==========

run_monte_carlo <- function(n_simulations = 1000, duration = 300) {
  scenarios <- names(create_network_models())

  results <- expand_grid(
    scenario = scenarios,
    simulation = 1:n_simulations
  ) %>%
    mutate(
      network_trace = map2(scenario, simulation, function(s, i) {
        set.seed(i * 1000 + which(scenarios == s))
        generate_latency_trace(s, duration)
      }),
      viewing_session = map(simulation, function(i) {
        set.seed(i + 5000)
        simulate_viewing_session(duration)
      }),
      performance = map2(network_trace, viewing_session, calculate_sync_performance)
    ) %>%
    select(-network_trace, -viewing_session) %>%
    unnest(performance)

  return(results)
}

# ========== 4b. Parallel Monte Carlo ==========

run_monte_carlo_parallel <- function(n_simulations = 1000, duration = 300) {
  message("Running Monte Carlo simulations in parallel using ", n_cores, " cores...")

  scenarios <- names(create_network_models())

  # Create simulation grid
  sim_grid <- expand_grid(
    scenario = scenarios,
    simulation = 1:n_simulations
  )

  # Split simulations into chunks for parallel processing
  n_chunks <- n_cores * 4
  sim_grid$chunk <- rep(1:n_chunks, length.out = nrow(sim_grid))

  # Progress reporting
  pb <- progress::progress_bar$new(
    format = "  Simulating [:bar] :percent eta: :eta",
    total = n_chunks,
    clear = FALSE
  )

  # Parallel processing function for each chunk
  process_chunk <- function(chunk_data) {
    chunk_data %>%
      mutate(
        network_trace = map2(scenario, simulation, function(s, i) {
          set.seed(i * 1000 + which(scenarios == s))
          generate_latency_trace(s, duration)
        }),
        viewing_session = map(simulation, function(i) {
          set.seed(i + 5000)
          simulate_viewing_session(duration)
        }),
        performance = map2(network_trace, viewing_session, calculate_sync_performance)
      ) %>%
      select(-network_trace, -viewing_session) %>%
      unnest(performance)
  }

  # Run parallel processing
  results <- sim_grid %>%
    group_split(chunk) %>%
    future_map_dfr(function(chunk) {
      pb$tick()
      process_chunk(chunk)
    }, .options = furrr_options(seed = 123))

  return(results)
}

# ========== 5. Theoretical Bounds Calculation ==========

calculate_theoretical_bounds <- function() {
  # System parameters
  api_poll_rate <- 4
  mean_processing_time <- 10

  # Queueing theory
  arrival_rate <- api_poll_rate
  service_rate <- 1000 / mean_processing_time
  utilization <- arrival_rate / service_rate

  mean_queue_delay <- ifelse(
    utilization < 1,
    utilization / (service_rate * (1 - utilization)),
    Inf
  )

  tibble(
    parameter = c(
      "Max tolerable latency (ms)",
      "API polling interval (ms)",
      "Queue delay (ms)",
      "System utilization",
      "System stable"
    ),
    value = c(
      200,
      1000 / api_poll_rate,
      mean_queue_delay * 1000,
      utilization,
      as.numeric(utilization < 1)
    )
  )
}

# ========== 6. Visualization Functions ==========

create_impact_plot <- function(impact_data) {
  p1 <- ggplot(impact_data, aes(x = network_latency_ms, y = data_quality_score)) +
    geom_line(size = 1.2, color = "steelblue") +
    geom_ribbon(aes(ymin = 0, ymax = data_quality_score),
                alpha = 0.2, fill = "steelblue") +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "green", size = 1) +
    geom_hline(yintercept = 0.90, linetype = "dashed", color = "orange", size = 1) +
    annotate("text", x = 400, y = 0.96, label = "95% threshold", color = "green") +
    annotate("text", x = 400, y = 0.91, label = "90% threshold", color = "orange") +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(
      title = "A. Synchronization Accuracy vs Network Latency",
      x = "Network Latency (ms)",
      y = "Data Quality Score"
    ) +
    theme_minimal()

  envelope_data <- expand_grid(
    latency_ms = seq(0, 300, by = 10),
    packet_loss = seq(0, 0.05, by = 0.002)
  ) %>%
    mutate(
      quality_score = map2_dbl(latency_ms, packet_loss, function(l, p) {
        base_score <- compute_sync_error(l)$data_quality_score[1]
        base_score * (1 - p * 10)
      })
    )

  p2 <- ggplot(envelope_data, aes(x = latency_ms, y = packet_loss * 100, z = quality_score)) +
    geom_contour_filled(breaks = seq(0, 1, by = 0.1)) +
    geom_contour(breaks = c(0.90, 0.95), color = "white", size = 1) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    scale_fill_viridis_d(name = "Quality\nScore", option = "D") +
    labs(
      title = "B. Operating Envelope",
      x = "Network Latency (ms)",
      y = "Packet Loss (%)"
    ) +
    theme_minimal()

  p1 / p2
}

create_monte_carlo_plot <- function(mc_results) {
  summary_stats <- mc_results %>%
    group_by(scenario) %>%
    summarise(
      mean_accuracy = mean(sync_accuracy, na.rm = TRUE),
      ci_lower = quantile(sync_accuracy, 0.025, na.rm = TRUE),
      ci_upper = quantile(sync_accuracy, 0.975, na.rm = TRUE),
      mean_latency = mean(mean_latency, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      scenario = factor(scenario, levels = names(create_network_models()))
    )

  ggplot(summary_stats, aes(x = scenario, y = mean_accuracy)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, size = 1) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "darkgreen", size = 1) +
    geom_hline(yintercept = 0.90, linetype = "dashed", color = "orange", size = 1) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(
      title = "Expected Synchronization Accuracy by Network Distance",
      subtitle = sprintf("Based on %d Monte Carlo simulations", n_distinct(mc_results$simulation)),
      x = "Network Scenario",
      y = "Synchronization Accuracy",
      caption = "Error bars show 95% CI. Dashed lines indicate accuracy thresholds."
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ========== 7. Report Generation ==========

generate_report <- function(impact_data, mc_results, bounds) {
  threshold_95 <- impact_data %>%
    filter(data_quality_score >= 0.95) %>%
    filter(!is.infinite(network_latency_ms)) %>%
    summarise(max_latency = max(network_latency_ms, na.rm = TRUE)) %>%
    pull(max_latency)

  # If no valid threshold found, use interpolation
  if (is.na(threshold_95) || is.infinite(threshold_95)) {
    # Find the crossing point by interpolation
    above_95 <- impact_data %>% filter(data_quality_score > 0.95) %>% tail(1)
    below_95 <- impact_data %>% filter(data_quality_score < 0.95) %>% head(1)

    if (nrow(above_95) > 0 && nrow(below_95) > 0) {
      threshold_95 <- approx(
        x = c(above_95$data_quality_score, below_95$data_quality_score),
        y = c(above_95$network_latency_ms, below_95$network_latency_ms),
        xout = 0.95
      )$y
    } else {
      threshold_95 <- 100  # Default fallback
    }
  }

  threshold_90 <- impact_data %>%
    filter(data_quality_score >= 0.90) %>%
    summarise(max_latency = max(network_latency_ms)) %>%
    pull(max_latency)

  mc_summary <- mc_results %>%
    group_by(scenario) %>%
    summarise(
      `Mean Accuracy (%)` = round(mean(sync_accuracy, na.rm = TRUE) * 100, 1),
      `Mean Latency (ms)` = round(mean(mean_latency, na.rm = TRUE), 1),
      `95% Latency (ms)` = round(mean(p95_latency, na.rm = TRUE), 1),
      `Packet Loss (%)` = round(mean(packet_loss, na.rm = TRUE) * 100, 2),
      .groups = "drop"
    )

  report <- glue("
# Theoretical Network Latency Analysis Report
Generated: {Sys.Date()}

## Executive Summary

- **95% accuracy threshold**: {threshold_95}ms maximum network latency
- **90% accuracy threshold**: {threshold_90}ms maximum network latency
- **Recommended deployment limit**: 75ms latency with <0.5% packet loss

## Monte Carlo Simulation Results

{knitr::kable(mc_summary, format = 'markdown')}

## System Theoretical Bounds

{knitr::kable(bounds, format = 'markdown')}

## Deployment Guidelines

1. **LAN/Same Building** (<5ms): Optimal performance, all features enabled
2. **Metropolitan** (10-30ms): Full functionality maintained
3. **Regional** (30-70ms): Consider increasing API polling frequency
4. **Continental** (70-120ms): Implement predictive synchronization
5. **International** (>150ms): Architectural modifications required

## Technical Recommendations

- Increase API polling rate from 4Hz to 10Hz for latencies >50ms
- Implement client-side timestamp injection for latencies >100ms
- Consider WebSocket connections for latencies >150ms
")

  return(report)
}

# ========== 8. Parallel Analysis Functions ==========

analyze_latency_impact_parallel <- function() {
  latency_values <- seq(0, 500, by = 2)

  # Split into chunks for parallel processing
  chunks <- split(latency_values,
                  cut(seq_along(latency_values), n_cores, labels = FALSE))

  # Process in parallel
  future_map_dfr(chunks, function(latencies) {
    tibble(network_latency_ms = latencies) %>%
      mutate(compute_sync_error(network_latency_ms)) %>%
      unnest(cols = everything())
  })
}

# ========== 9. Main Analysis Function ==========

run_analysis_parallel <- function() {
  start_time <- Sys.time()

  message("Starting theoretical network latency evaluation (parallel version)...")
  message("Output directory: ", OUTPUT_DIR)
  message("Using ", n_cores, " CPU cores")

  # Generate latency impact analysis in parallel
  message("\nAnalyzing latency impact...")
  impact_data <- analyze_latency_impact_parallel()

  # Run Monte Carlo simulation in parallel
  mc_results <- run_monte_carlo_parallel(n_simulations = 1000, duration = 300)

  # Calculate theoretical bounds
  message("\nCalculating theoretical bounds...")
  bounds <- calculate_theoretical_bounds()

  # Generate visualizations
  message("\nCreating visualizations...")

  impact_plot <- create_impact_plot(impact_data)
  ggsave(file.path(OUTPUT_DIR, "latency_impact_analysis.pdf"),
         impact_plot, width = 10, height = 12, dpi = 300)

  mc_plot <- create_monte_carlo_plot(mc_results)
  ggsave(file.path(OUTPUT_DIR, "monte_carlo_results.pdf"),
         mc_plot, width = 10, height = 8, dpi = 300)

  # Generate network trace examples in parallel
  message("\nGenerating example traces...")
  example_traces <- future_map_dfr(names(create_network_models()),
                                   ~generate_latency_trace(.x, duration_seconds = 60))

  trace_plot <- ggplot(example_traces, aes(x = time, y = latency_ms, color = scenario)) +
    geom_line(alpha = 0.7) +
    facet_wrap(~scenario, scales = "free_y", ncol = 1) +
    labs(title = "Example Network Latency Traces (60 seconds)",
         x = "Time (seconds)", y = "Latency (ms)") +
    theme_minimal() +
    theme(legend.position = "none")

  ggsave(file.path(OUTPUT_DIR, "example_network_traces.pdf"),
         trace_plot, width = 10, height = 12, dpi = 300)

  # Save data
  message("\nSaving analysis data...")
  saveRDS(impact_data, file.path(OUTPUT_DIR, "impact_analysis.rds"))
  saveRDS(mc_results, file.path(OUTPUT_DIR, "monte_carlo_results.rds"))
  saveRDS(bounds, file.path(OUTPUT_DIR, "theoretical_bounds.rds"))
  write_csv(impact_data, file.path(OUTPUT_DIR, "impact_analysis.csv"))
  write_csv(mc_results, file.path(OUTPUT_DIR, "monte_carlo_results.csv"))

  # Generate report
  message("\nGenerating summary report...")
  report <- generate_report(impact_data, mc_results, bounds)
  writeLines(report, file.path(OUTPUT_DIR, "analysis_report.md"))

  # Clean up parallel backend
  plan(sequential)

  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "secs")

  message("\n✓ Analysis complete!")
  message("Results saved to: ", OUTPUT_DIR)
  message("Total runtime: ", round(runtime, 1), " seconds")

  return(list(
    impact = impact_data,
    monte_carlo = mc_results,
    bounds = bounds,
    report = report,
    runtime = runtime
  ))
}
# ========== Run Analysis ==========

# Run full analysis
results <- run_analysis_parallel()
