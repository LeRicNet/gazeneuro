#!/usr/bin/env Rscript
# ==============================================================================
# SRCV Framework Validation Study
# ==============================================================================
# This script implements a comprehensive validation study for the Synchronized
# Residual Cross-Validation (SRCV) latency estimation framework
#
# To run different experiments, change the 'experiment_number' variable
# All outputs are saved under: results/experiment{N}/
#
# Author: [Your Name]
# Date: 2025-01-10
# Version: 1.1
# ==============================================================================

# Load required libraries ------------------------------------------------------
library(tidyverse)
library(parallel)
library(lme4)
library(mgcv)
library(pROC)
library(plotly)
library(viridis)
library(ggpubr)
library(knitr)
library(kableExtra)
library(signal)
library(pracma)

# Set global parameters --------------------------------------------------------
set.seed(42)  # For reproducibility
n_cores <- detectCores() - 1  # Leave one core free

# TEST MODE: Set to TRUE for a quick test run, FALSE for full experiment
TEST_MODE <- FALSE  # Change to TRUE for testing

if (TEST_MODE) {
  n_replicates <- 2  # Just 2 replicates for testing
  message("*** RUNNING IN TEST MODE - Results will be limited ***")
} else {
  n_replicates <- 50  # Full 50 replicates
}

save_interval <- 1000  # Save results every N iterations

# Define experiment number (easy to change for different experiments)
experiment_number <- 5  # Change this for different experiment runs

# Define base output directory
base_dir <- sprintf("results/experiment%d", experiment_number)

# Create output directories
dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "figures"), showWarnings = FALSE)
dir.create(file.path(base_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(base_dir, "diagnostics"), showWarnings = FALSE)

# Define experimental parameters -----------------------------------------------
# In test mode, use smaller parameter grid
if (TEST_MODE) {
  params <- expand.grid(
    true_latency = c(0, 50, 200),  # Reduced
    snr = c(0, 10, 20),  # Reduced
    n_events = c(5, 20, 50),  # Reduced
    sampling_freq = c(60, 120),  # Reduced
    replicate = 1:n_replicates
  )
} else {
  params <- expand.grid(
    true_latency = c(0, 20, 50, 100, 200, 500),  # ms
    snr = c(0, 3, 6, 10, 15, 20, 30),  # dB
    n_events = c(3, 5, 10, 20, 50, 100),
    sampling_freq = c(20, 60, 120, 240),  # Hz
    replicate = 1:n_replicates
  )
}

message(sprintf("Total simulations: %d", nrow(params)))
message(sprintf("Output directory: %s", base_dir))
if (TEST_MODE) {
  message("*** TEST MODE: Using reduced parameter grid ***")
}

# ==============================================================================
# DATA GENERATION FUNCTIONS
# ==============================================================================

#' Generate viewing events using Poisson process
#'
#' @param duration Total experiment duration in seconds
#' @param lambda Event rate (events per second)
#' @param n_events Number of events to generate
#' @return Data frame with event start times and durations
generate_viewing_events <- function(duration = 200, lambda = 0.5, n_events = NULL) {
  if (!is.null(n_events)) {
    # Generate events uniformly distributed over duration
    # Or use exponential but truncate to duration
    start_times <- sort(runif(n_events, min = 0, max = duration * 0.8))  # Keep some buffer
  } else {
    # Generate based on Poisson process
    n_expected <- duration * lambda
    inter_arrivals <- rexp(n_expected * 2, rate = lambda)
    start_times <- cumsum(inter_arrivals)
    start_times <- start_times[start_times < duration]
  }

  # Generate event durations: Normal(2000ms, 500ms), truncated [500, 5000]
  durations <- pmax(500, pmin(5000, rnorm(length(start_times), 2000, 500)))

  # Make sure events don't exceed duration
  valid_events <- (start_times * 1000 + durations) < (duration * 1000)

  data.frame(
    event_id = seq_along(start_times)[valid_events],
    start_time = start_times[valid_events] * 1000,  # Convert to ms
    duration = durations[valid_events],
    end_time = start_times[valid_events] * 1000 + durations[valid_events]
  )
}

#' Generate correlated gaze response
#'
#' @param events Data frame with event timings
#' @param snr_db Target signal-to-noise ratio in dB
#' @param fs Sampling frequency in Hz
#' @param latency True latency to apply (ms)
#' @return List with time series and ground truth
generate_gaze_response <- function(events, snr_db, fs, latency) {
  # Create time vector
  duration <- max(events$end_time) + 5000  # Add buffer
  dt <- 1000 / fs  # Sampling period in ms
  time <- seq(0, duration, by = dt)
  n_samples <- length(time)

  # Initialize base signal (step function at event transitions)
  base_signal <- rep(0, n_samples)

  for (i in 1:nrow(events)) {
    # Find samples during this event
    event_start_idx <- which.min(abs(time - (events$start_time[i] + latency)))
    event_end_idx <- which.min(abs(time - (events$end_time[i] + latency)))

    if (event_start_idx <= n_samples && event_end_idx <= n_samples) {
      base_signal[event_start_idx:event_end_idx] <- 1

      # Add exponential decay response (physiological component)
      tau <- 150  # ms decay constant
      decay_duration <- min(500, event_end_idx - event_start_idx)
      if (decay_duration > 0) {
        decay_time <- seq(0, decay_duration * dt, by = dt)
        decay_response <- exp(-decay_time / tau)
        decay_idx <- event_start_idx:(event_start_idx + length(decay_response) - 1)
        decay_idx <- decay_idx[decay_idx <= n_samples]
        base_signal[decay_idx] <- base_signal[decay_idx] + decay_response[1:length(decay_idx)]
      }
    }
  }

  # Normalize signal
  base_signal <- (base_signal - mean(base_signal)) / sd(base_signal)

  # Add noise to achieve target SNR
  signal_power <- mean(base_signal^2)
  snr_linear <- 10^(snr_db / 10)
  noise_power <- signal_power / snr_linear
  noise <- rnorm(n_samples, 0, sqrt(noise_power))

  # Combine signal and noise
  gaze_signal <- base_signal + noise

  # Add temporal jitter (±5% of sampling period)
  jitter <- runif(n_samples, -0.05, 0.05) * dt

  list(
    time = time + jitter,
    gaze = gaze_signal,
    base_signal = base_signal,
    events = events,
    true_latency = latency,
    snr_db = snr_db,
    fs = fs
  )
}

# ==============================================================================
# SRCV LATENCY ESTIMATION FUNCTIONS
# ==============================================================================

#' Compute cross-correlation for latency estimation
#'
#' @param gaze_signal Gaze time series
#' @param event_signal Event indicator time series
#' @param max_lag Maximum lag to consider (samples)
#' @return List with correlation values and lag indices
compute_cross_correlation <- function(gaze_signal, event_signal, max_lag = NULL) {
  n <- length(gaze_signal)
  if (is.null(max_lag)) {
    max_lag <- round(n / 4)
  }

  # Normalize signals
  gaze_norm <- (gaze_signal - mean(gaze_signal)) / sd(gaze_signal)
  event_norm <- (event_signal - mean(event_signal)) / sd(event_signal)

  # Compute cross-correlation
  lags <- -max_lag:max_lag
  correlations <- numeric(length(lags))

  for (i in seq_along(lags)) {
    lag <- lags[i]
    if (lag < 0) {
      # Event leads gaze
      idx_event <- 1:(n + lag)
      idx_gaze <- (-lag + 1):n
    } else {
      # Gaze leads event
      idx_event <- (lag + 1):n
      idx_gaze <- 1:(n - lag)
    }
    correlations[i] <- cor(gaze_norm[idx_gaze], event_norm[idx_event])
  }

  optimal_lag = -1 * lags[which.max(abs(correlations))]

  list(
    lags = lags,
    correlations = correlations,
    max_correlation = max(abs(correlations)),
    optimal_lag = optimal_lag
  )
}

#' Estimate latency using SRCV method
#'
#' @param gaze_data List from generate_gaze_response
#' @return List with estimation results
estimate_latency_srcv <- function(gaze_data) {
  # Create event indicator signal
  time <- gaze_data$time
  events <- gaze_data$events
  event_signal <- rep(0, length(time))

  for (i in 1:nrow(events)) {
    # event_idx <- which(time >= events$start_time[i] & time <= events$end_time[i])
    event_idx <- which(time > events$start_time[i] & time < events$end_time[i])
    event_signal[event_idx] <- 1
  }

  # Compute cross-correlation
  max_lag <- round(1000 / (1000 / gaze_data$fs))  # 1 second max lag
  xcorr <- compute_cross_correlation(gaze_data$gaze, event_signal, max_lag)

  # Estimate latency
  estimated_latency <- xcorr$optimal_lag * (1000 / gaze_data$fs)

  # Compute prominence ratio
  correlations_abs <- abs(xcorr$correlations)
  peak_idx <- which.max(correlations_abs)
  peak_value <- correlations_abs[peak_idx]
  other_values <- correlations_abs[-peak_idx]
  prominence_ratio <- peak_value / median(other_values)

  # Estimate uncertainty (simplified version of Eq. 20)
  n_effective <- nrow(events)
  correlation_peak <- xcorr$max_correlation
  sigma_delta <- (1000 / gaze_data$fs) / (sqrt(n_effective) * correlation_peak)

  list(
    estimated_latency = estimated_latency,
    true_latency = gaze_data$true_latency,
    error = abs(estimated_latency - gaze_data$true_latency),
    prominence_ratio = prominence_ratio,
    correlation_peak = correlation_peak,
    sigma_delta = sigma_delta,
    n_events = nrow(events),
    snr = gaze_data$snr_db,
    fs = gaze_data$fs
  )
}

#' Bootstrap uncertainty estimation
#'
#' @param gaze_data Original gaze data
#' @param B Number of bootstrap samples
#' @return Bootstrap standard deviation
bootstrap_uncertainty <- function(gaze_data, B = 100) {
  bootstrap_estimates <- numeric(B)
  n <- length(gaze_data$gaze)

  for (b in 1:B) {
    # Resample with replacement
    idx <- sample(1:n, n, replace = TRUE)
    gaze_data_boot <- gaze_data
    gaze_data_boot$gaze <- gaze_data$gaze[idx]

    # Estimate latency on bootstrap sample
    result <- estimate_latency_srcv(gaze_data_boot)
    bootstrap_estimates[b] <- result$estimated_latency
  }

  sd(bootstrap_estimates)
}

# ==============================================================================
# MAIN SIMULATION FUNCTION
# ==============================================================================

#' Run single simulation
#'
#' @param params_row Single row from parameter grid
#' @return Data frame with results
run_simulation <- function(params_row) {
  tryCatch({
    # Generate viewing events
    events <- generate_viewing_events(
      duration = 200,
      n_events = params_row$n_events
    )

    # Generate gaze response
    gaze_data <- generate_gaze_response(
      events = events,
      snr_db = params_row$snr,
      fs = params_row$sampling_freq,
      latency = params_row$true_latency
    )

    # Estimate latency
    result <- estimate_latency_srcv(gaze_data)

    # Bootstrap uncertainty (only for subset to save time)
    if (params_row$replicate <= 5) {
      result$sigma_delta_boot <- bootstrap_uncertainty(gaze_data, B = 100)
    } else {
      result$sigma_delta_boot <- NA
    }

    # Create output data frame - explicitly convert to data.frame
    output_df <- data.frame(
      true_latency = params_row$true_latency,
      snr = params_row$snr,
      n_events = params_row$n_events,
      sampling_freq = params_row$sampling_freq,
      replicate = params_row$replicate,
      estimated_latency = result$estimated_latency,
      error = result$error,
      prominence_ratio = result$prominence_ratio,
      correlation_peak = result$correlation_peak,
      sigma_delta = result$sigma_delta,
      sigma_delta_boot = result$sigma_delta_boot,
      stringsAsFactors = FALSE
    )

    return(output_df)

  }, error = function(e) {
    # Return NA results on error - also as data.frame
    output_df <- data.frame(
      true_latency = params_row$true_latency,
      snr = params_row$snr,
      n_events = params_row$n_events,
      sampling_freq = params_row$sampling_freq,
      replicate = params_row$replicate,
      estimated_latency = NA,
      error = NA,
      prominence_ratio = NA,
      correlation_peak = NA,
      sigma_delta = NA,
      sigma_delta_boot = NA,
      stringsAsFactors = FALSE
    )

    return(output_df)
  })
}

# ==============================================================================
# RUN SIMULATIONS
# ==============================================================================

# Initialize timing variables
start_time <- Sys.time()
end_time <- start_time  # Default in case we load existing results

# Check if we should load existing results or run new simulations
existing_results_file <- file.path(base_dir, "simulation_results.rds")
if (file.exists(existing_results_file) && !TEST_MODE) {
  response <- readline(prompt = "Existing results found. Load them instead of re-running? (y/n): ")
  if (tolower(response) == "y") {
    message("Loading existing results...")
    results <- readRDS(existing_results_file)
    message(sprintf("Loaded %d simulation results", nrow(results)))
    end_time <- Sys.time()  # Set end time for loaded results
  } else {
    # Run new simulations
    message("Starting new simulations...")
    start_time <- Sys.time()

    # Run simulations in parallel
    if (n_cores > 1) {
      cl <- makeCluster(n_cores)
      clusterEvalQ(cl, {
        library(tidyverse)
        library(signal)
      })
      clusterExport(cl, c("generate_viewing_events", "generate_gaze_response",
                          "compute_cross_correlation", "estimate_latency_srcv",
                          "bootstrap_uncertainty", "run_simulation"))

      # Split params into list for parallel processing
      params_list <- split(params, seq(nrow(params)))
      results_list <- parLapply(cl, params_list, function(x) run_simulation(x))
      stopCluster(cl)
    } else {
      params_list <- split(params, seq(nrow(params)))
      results_list <- lapply(params_list, function(x) run_simulation(x))
    }

    # Combine results with error checking
    message("Combining results...")
    results <- bind_rows(results_list)

    # Check if results were properly created
    if (is.null(results) || nrow(results) == 0) {
      stop("No valid results were generated. Check the simulation functions.")
    }

    message(sprintf("Successfully processed %d/%d simulations",
                    sum(!is.na(results$error)), nrow(results)))

    # Save intermediate results
    saveRDS(results, file.path(base_dir, "simulation_results.rds"))
    write.csv(results, file.path(base_dir, "simulation_results.csv"), row.names = FALSE)

    end_time <- Sys.time()
    message(sprintf("Simulations completed in %.2f minutes",
                    as.numeric(end_time - start_time, units = "mins")))
  }
} else {
  # No existing results or in test mode - run simulations
  message("Starting simulations...")
  start_time <- Sys.time()

  # Run simulations in parallel
  if (n_cores > 1) {
    cl <- makeCluster(n_cores)
    clusterEvalQ(cl, {
      library(tidyverse)
      library(signal)
    })
    clusterExport(cl, c("generate_viewing_events", "generate_gaze_response",
                        "compute_cross_correlation", "estimate_latency_srcv",
                        "bootstrap_uncertainty", "run_simulation"))

    # Split params into list for parallel processing
    params_list <- split(params, seq(nrow(params)))
    results_list <- parLapply(cl, params_list, function(x) run_simulation(x))
    stopCluster(cl)
  } else {
    params_list <- split(params, seq(nrow(params)))
    results_list <- lapply(params_list, function(x) run_simulation(x))
  }

  # Combine results with error checking
  message("Combining results...")
  results <- bind_rows(results_list)

  # Check if results were properly created
  if (is.null(results) || nrow(results) == 0) {
    stop("No valid results were generated. Check the simulation functions.")
  }

  message(sprintf("Successfully processed %d/%d simulations",
                  sum(!is.na(results$error)), nrow(results)))

  # Save intermediate results
  saveRDS(results, file.path(base_dir, "simulation_results.rds"))
  write.csv(results, file.path(base_dir, "simulation_results.csv"), row.names = FALSE)

  end_time <- Sys.time()
  message(sprintf("Simulations completed in %.2f minutes",
                  as.numeric(end_time - start_time, units = "mins")))
}

# ==============================================================================
# STATISTICAL ANALYSIS
# ==============================================================================

message("Running statistical analyses...")

# Ensure results is a proper data frame
if (!is.data.frame(results)) {
  stop("Results object is not a data frame. Check simulation output.")
}

# Check for required columns
required_cols <- c("error", "n_events", "snr", "sampling_freq", "true_latency",
                   "prominence_ratio", "sigma_delta")
missing_cols <- setdiff(required_cols, names(results))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns in results:", paste(missing_cols, collapse = ", ")))
}

# Clean results - use dplyr explicitly to avoid conflicts
results_clean <- dplyr::filter(results, !is.na(error))

if (nrow(results_clean) == 0) {
  stop("No valid results after filtering. All simulations may have failed.")
}

results_clean <- results_clean %>%
  dplyr::mutate(
    n_snr_product = n_events * snr,
    log_error = log(error + 1),  # Add 1 to handle zeros
    log_n = log(n_events),
    log_snr = log(snr + 1),  # Add 1 to handle zeros
    log_n_snr = log(n_snr_product + 1),
    log_fs = log(sampling_freq),
    success = error < 50  # Binary success indicator
  )

message(sprintf("Analyzing %d valid simulations", nrow(results_clean)))

# Model 1: Power function relationship (H1) -----------------------------------
model1 <- lm(log_error ~ log_n + log_snr + log_n_snr + log_fs,
             data = results_clean)

summary_model1 <- summary(model1)
confint_model1 <- confint(model1)

# Test H1: coefficient of log(N×SNR) < -0.5
h1_test <- coef(model1)["log_n_snr"] < -0.5
h1_pvalue <- pt(coef(model1)["log_n_snr"] / summary_model1$coefficients["log_n_snr", "Std. Error"],
                df = model1$df.residual)

# Calculate success probability for N×SNR > 60
results_clean$success_nsnr60 <- results_clean$n_snr_product > 60
success_rate_nsnr60 <- mean(results_clean$success[results_clean$success_nsnr60])

# Model 2: Quality metric validation (H2) -------------------------------------
# Check if we have enough data for the model
if (length(unique(results_clean$success)) > 1 && length(unique(results_clean$true_latency)) > 1) {
  model2 <- glmer(success ~ prominence_ratio + (1|true_latency),
                  data = results_clean,
                  family = binomial,
                  control = glmerControl(optimizer = "bobyqa"))

  # ROC analysis
  roc_prominence <- roc(results_clean$success, results_clean$prominence_ratio)
  optimal_threshold <- coords(roc_prominence, "best", ret = "threshold")
  auc_value <- auc(roc_prominence)

  # Test H2: AUC > 0.90
  h2_test <- auc_value > 0.90
} else {
  message("Warning: Not enough variation in data for Model 2. Using simplified analysis.")
  model2 <- glm(success ~ prominence_ratio, data = results_clean, family = binomial)
  roc_prominence <- roc(results_clean$success, results_clean$prominence_ratio)
  optimal_threshold <- list(threshold = median(results_clean$prominence_ratio, na.rm = TRUE))
  auc_value <- 0.5  # Default for insufficient data
  h2_test <- FALSE
}

# Model 3: Uncertainty calibration (H3) ---------------------------------------
uncertainty_data <- results_clean %>%
  as_tibble() %>%
  dplyr::filter(!is.na(sigma_delta_boot)) %>%
  mutate(
    log_sigma_boot = log(sigma_delta_boot + 1),
    log_sigma_analytical = log(sigma_delta + 1)
  )

# Check if we have enough data for uncertainty calibration
if (nrow(uncertainty_data) > 2) {
  model3 <- lm(log_sigma_boot ~ log_sigma_analytical, data = uncertainty_data)
  r_squared_model3 <- summary(model3)$r.squared

  if ("log_sigma_analytical" %in% names(coef(model3))) {
    slope_model3 <- coef(model3)["log_sigma_analytical"]
  } else {
    slope_model3 <- NA
  }

  # Test H3: slope ∈ [0.95, 1.05] and R² > 0.90
  h3_test <- !is.na(slope_model3) && (slope_model3 >= 0.95 & slope_model3 <= 1.05) & (r_squared_model3 > 0.90)
} else {
  message("Warning: Not enough bootstrap data for Model 3. Skipping uncertainty calibration.")
  # Create a dummy model for compatibility
  model3 <- lm(log_sigma_boot ~ 1, data = data.frame(log_sigma_boot = 0))
  r_squared_model3 <- 0
  slope_model3 <- NA
  h3_test <- FALSE
}

# Model 4: Minimum requirements analysis (H4) ---------------------------------
# Check if we have enough data points for GAM
if (nrow(results_clean) > 30) {
  model4 <- gam(success ~ s(n_events, snr, sampling_freq, k = 3),
                data = results_clean,
                family = binomial)
} else {
  # Use simpler model for small datasets
  model4 <- glm(success ~ n_events * snr * sampling_freq,
                data = results_clean,
                family = binomial)
}

# Test H4: P(success | N=3, SNR=20, fs>60) > 0.90
test_data_h4 <- data.frame(
  n_events = 3,
  snr = 20,
  sampling_freq = c(120, 240)
)

pred_h4 <- predict(model4, newdata = test_data_h4, type = "response")
h4_test <- all(pred_h4 > 0.90)

# ==============================================================================
# GENERATE OUTPUTS
# ==============================================================================

message("Generating figures and tables...")

# Figure 1: 3D surface plot of error as function of N and SNR ----------------
fig1_data <- results_clean %>%
  group_by(n_events, snr, sampling_freq) %>%
  summarise(
    mean_error = mean(error),
    se_error = sd(error) / sqrt(n()),
    .groups = "drop"
  )

fig1 <- plot_ly(
  data = fig1_data,
  x = ~n_events,
  y = ~snr,
  z = ~mean_error,
  color = ~factor(sampling_freq),
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 5)
) %>%
  layout(
    title = "Latency Estimation Error as Function of N and SNR",
    scene = list(
      xaxis = list(title = "Number of Events (N)"),
      yaxis = list(title = "SNR (dB)"),
      zaxis = list(title = "Mean Error (ms)")
    )
  )

htmlwidgets::saveWidget(fig1, file.path(base_dir, "figures", "figure1_3d_surface.html"))

# Figure 2: ROC curves for prominence ratio -----------------------------------
fig2_data <- list()
for (snr_level in unique(results_clean$snr)) {
  subset_data <- results_clean %>% dplyr::filter(snr == snr_level)
  roc_obj <- roc(subset_data$success, subset_data$prominence_ratio)
  fig2_data[[as.character(snr_level)]] <- data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities,
    snr = snr_level
  )
}

fig2 <- ggplot(bind_rows(fig2_data), aes(x = fpr, y = tpr, color = factor(snr))) +
  geom_line(size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5) +
  labs(
    title = "ROC Curves for Prominence Ratio at Different SNR Levels",
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = "SNR (dB)"
  ) +
  theme_minimal() +
  scale_color_viridis_d()

ggsave(file.path(base_dir, "figures", "figure2_roc_curves.png"), fig2, width = 10, height = 8, dpi = 300)

# Figure 3: Scatter plot of analytical vs bootstrap uncertainty ---------------
if (nrow(uncertainty_data) > 2 && !is.na(slope_model3)) {
  fig3 <- ggplot(uncertainty_data, aes(x = sigma_delta, y = sigma_delta_boot)) +
    geom_point(alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    geom_smooth(method = "lm", se = TRUE) +
    labs(
      title = "Analytical vs Bootstrap Uncertainty Estimates",
      x = "Analytical σδ (ms)",
      y = "Bootstrap σδ (ms)",
      subtitle = sprintf("R² = %.3f, Slope = %.3f", r_squared_model3, slope_model3)
    ) +
    theme_minimal() +
    coord_equal()
} else {
  # Create a placeholder plot if not enough data
  fig3 <- ggplot(data.frame(x = 1, y = 1), aes(x = x, y = y)) +
    geom_point(size = 0) +
    labs(
      title = "Analytical vs Bootstrap Uncertainty Estimates",
      subtitle = "Insufficient data for comparison",
      x = "Analytical σδ (ms)",
      y = "Bootstrap σδ (ms)"
    ) +
    theme_minimal()
}

ggsave(file.path(base_dir, "figures", "figure3_uncertainty_calibration.png"), fig3, width = 8, height = 8, dpi = 300)

# Figure 4: Heatmap of success probability ------------------------------------
# Create prediction grid
pred_grid <- expand.grid(
  n_events = seq(3, 100, length.out = 20),
  snr = seq(0, 30, length.out = 20),
  sampling_freq = c(20, 60, 120, 240)
)

pred_grid$success_prob <- predict(model4, newdata = pred_grid, type = "response")

fig4 <- ggplot(pred_grid, aes(x = n_events, y = snr, fill = success_prob)) +
  geom_tile() +
  facet_wrap(~sampling_freq, labeller = label_both) +
  scale_fill_viridis(name = "P(Success)", limits = c(0, 1)) +
  labs(
    title = "Success Probability Across Parameter Space",
    x = "Number of Events",
    y = "SNR (dB)"
  ) +
  theme_minimal() +
  geom_contour(aes(z = success_prob), color = "white", alpha = 0.5,
               breaks = c(0.5, 0.9, 0.95))

ggsave(file.path(base_dir, "figures", "figure4_success_heatmap.png"), fig4, width = 12, height = 10, dpi = 300)

# Table 1: Regression coefficients --------------------------------------------
# Prepare coefficient data depending on model types
if (class(model2)[1] == "glmerMod") {
  model2_coefs <- fixef(model2)
  model2_ci <- suppressWarnings(confint(model2, method = "Wald", parm = "beta_"))
} else {
  model2_coefs <- coef(model2)
  model2_ci <- suppressWarnings(confint(model2))
}

# Ensure CI dimensions match coefficients
if (nrow(model2_ci) != length(model2_coefs)) {
  model2_ci <- matrix(rep(c(NA, NA), length(model2_coefs)), ncol = 2, byrow = TRUE)
  rownames(model2_ci) <- names(model2_coefs)
}

# Handle model3 coefficients
if ("log_sigma_analytical" %in% names(coef(model3))) {
  model3_coefs <- coef(model3)
  model3_ci <- suppressWarnings(confint(model3))
} else {
  # Dummy values for insufficient data case
  model3_coefs <- c("(Intercept)" = 0, "log_sigma_analytical" = NA)
  model3_ci <- matrix(c(0, 0, NA, NA), ncol = 2)
  rownames(model3_ci) <- names(model3_coefs)
}

table1 <- data.frame(
  Model = c(rep("Model 1: Power Function", length(coef(model1))),
            rep("Model 2: Quality Metric", length(model2_coefs)),
            rep("Model 3: Uncertainty", length(model3_coefs))),
  Coefficient = c(names(coef(model1)),
                  names(model2_coefs),
                  names(model3_coefs)),
  Estimate = c(coef(model1),
               model2_coefs,
               model3_coefs),
  CI_Lower = c(confint_model1[,1],
               model2_ci[,1],
               model3_ci[,1]),
  CI_Upper = c(confint_model1[,2],
               model2_ci[,2],
               model3_ci[,2]),
  stringsAsFactors = FALSE
) %>%
  mutate(across(where(is.numeric), ~round(., 3)))

write.csv(table1, file.path(base_dir, "tables", "table1_regression_coefficients.csv"), row.names = FALSE)

kable(table1, format = "html", caption = "Table 1: Regression Coefficients with 95% CIs") %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  save_kable(file.path(base_dir, "tables", "table1_regression_coefficients.html"))

# Table 2: Confusion matrix ---------------------------------------------------
results_clean$predicted_success <- results_clean$prominence_ratio > optimal_threshold$threshold

confusion_matrix <- table(
  Actual = results_clean$success,
  Predicted = results_clean$predicted_success
)

table2 <- as.data.frame.matrix(confusion_matrix) %>%
  mutate(
    Actual = c("Failure", "Success"),
    .before = 1
  )

# Add performance metrics
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
sensitivity <- confusion_matrix[2,2] / sum(confusion_matrix[2,])
specificity <- confusion_matrix[1,1] / sum(confusion_matrix[1,])

metrics_row <- data.frame(
  Actual = "Performance",
  `FALSE` = sprintf("Acc: %.3f", accuracy),
  `TRUE` = sprintf("Sens: %.3f, Spec: %.3f", sensitivity, specificity)
)
names(metrics_row) <- names(table2)

table2 <- bind_cols(table2, metrics_row)

write.csv(table2, file.path(base_dir, "tables", "table2_confusion_matrix.csv"), row.names = FALSE)

kable(table2, format = "html",
      caption = sprintf("Table 2: Confusion Matrix (Threshold = %.2f)", optimal_threshold$threshold)) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  save_kable(file.path(base_dir, "tables", "table2_confusion_matrix.html"))

# ==============================================================================
# GENERATE DIAGNOSTIC PLOTS FOR FAILED SYNCHRONIZATIONS
# ==============================================================================

# Identify failed synchronizations
failed_sync <- results_clean %>%
  dplyr::filter(error > 100) %>%
  slice_sample(n = 10)

if (nrow(failed_sync) > 0) {
  message("Generating diagnostic plots for failed synchronizations...")

  for (i in 1:nrow(failed_sync)) {
    # Recreate the failed case
    set.seed(failed_sync$replicate[i])

    events <- generate_viewing_events(
      duration = 200,
      n_events = failed_sync$n_events[i]
    )

    gaze_data <- generate_gaze_response(
      events = events,
      snr_db = failed_sync$snr[i],
      fs = failed_sync$sampling_freq[i],
      latency = failed_sync$true_latency[i]
    )

    # Create diagnostic plot
    p1 <- ggplot(data.frame(time = gaze_data$time, signal = gaze_data$gaze)) +
      geom_line(aes(x = time/1000, y = signal), alpha = 0.7) +
      labs(title = sprintf("Failed Sync: Error = %.1f ms", failed_sync$error[i]),
           subtitle = sprintf("True δ = %d ms, SNR = %d dB, N = %d, fs = %d Hz",
                              failed_sync$true_latency[i], failed_sync$snr[i],
                              failed_sync$n_events[i], failed_sync$sampling_freq[i]),
           x = "Time (s)", y = "Gaze Signal") +
      theme_minimal()

    # Add event markers
    for (j in 1:nrow(events)) {
      p1 <- p1 +
        geom_vline(xintercept = events$start_time[j]/1000, color = "red", alpha = 0.3) +
        geom_vline(xintercept = events$end_time[j]/1000, color = "blue", alpha = 0.3)
    }

    ggsave(file.path(base_dir, "diagnostics", sprintf("failed_sync_%d.png", i)), p1, width = 12, height = 6, dpi = 150)
  }
}

# ==============================================================================
# SUMMARY REPORT
# ==============================================================================

# Calculate runtime properly
runtime_mins <- as.numeric(difftime(end_time, start_time, units = "mins"))
sims_per_sec <- ifelse(runtime_mins > 0,
                       nrow(params) / (runtime_mins * 60),
                       NA)

summary_report <- sprintf("
================================================================================
SRCV FRAMEWORK VALIDATION STUDY - SUMMARY REPORT
================================================================================
Date: %s
Mode: %s
Total Simulations: %d
Successful Runs: %d (%.1f%%)

HYPOTHESIS TESTING RESULTS:
--------------------------------------------------------------------------------
H1: Power function relationship
   - Coefficient of log(N×SNR): %.3f %s
   - Success rate when N×SNR > 60: %.1f%%
   - TEST RESULT: %s

H2: Prominence ratio validation
   - AUC: %.3f
   - Optimal threshold: %.2f
   - TEST RESULT: %s

H3: Uncertainty calibration
   - Slope: %s
   - R²: %.3f
   - TEST RESULT: %s

H4: Minimum requirements
   - P(success | N=3, SNR=20, fs>60): %.1f%%
   - TEST RESULT: %s

COMPUTATIONAL PERFORMANCE:
--------------------------------------------------------------------------------
- Total runtime: %.2f minutes
- Simulations per second: %s
- Cores used: %d

FILES GENERATED:
--------------------------------------------------------------------------------
- Results: %s/simulation_results.csv
- Figure 1: %s/figures/figure1_3d_surface.html
- Figure 2: %s/figures/figure2_roc_curves.png
- Figure 3: %s/figures/figure3_uncertainty_calibration.png
- Figure 4: %s/figures/figure4_success_heatmap.png
- Table 1: %s/tables/table1_regression_coefficients.html
- Table 2: %s/tables/table2_confusion_matrix.html
- Diagnostics: %s/diagnostics/failed_sync_*.png

================================================================================
",
                          Sys.Date(),
                          ifelse(TEST_MODE, "TEST", "FULL"),
                          nrow(params),
                          nrow(results_clean),
                          100 * nrow(results_clean) / nrow(params),
                          coef(model1)["log_n_snr"],
                          ifelse(coef(model1)["log_n_snr"] < -0.5, "✓ (< -0.5)", "✗ (≥ -0.5)"),
                          100 * success_rate_nsnr60,
                          ifelse(h1_test & success_rate_nsnr60 > 0.95, "PASS", "FAIL"),
                          auc_value,
                          optimal_threshold$threshold,
                          ifelse(h2_test, "PASS", "FAIL"),
                          ifelse(is.na(slope_model3), "N/A", sprintf("%.3f", slope_model3)),
                          r_squared_model3,
                          ifelse(h3_test, "PASS", "FAIL"),
                          100 * mean(pred_h4),
                          ifelse(h4_test, "PASS", "FAIL"),
                          runtime_mins,
                          ifelse(is.na(sims_per_sec), "N/A", sprintf("%.1f", sims_per_sec)),
                          n_cores,
                          base_dir, base_dir, base_dir, base_dir, base_dir, base_dir, base_dir, base_dir
)

cat(summary_report)
writeLines(summary_report, file.path(base_dir, "summary_report.txt"))

message("\nValidation study completed successfully!")
message(sprintf("Check the %s directory for all outputs.", base_dir))
