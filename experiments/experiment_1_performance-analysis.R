# Calculate Performance Metrics for Experiment 1
# Analyzes integration results to compute all metrics from Methods

library(tidyverse)
library(gazeneuro)

# Load processed data
output_dir <- "experiment1_data"
all_sessions <- readRDS(file.path(output_dir, "all_sessions.rds"))
integration_results <- readRDS(file.path(output_dir, "integration_results.rds"))
processing_log <- read.csv(file.path(output_dir, "processing_log.csv"))

# Function to calculate spatial accuracy
calculate_spatial_accuracy <- function(original_gaze, integrated, z_axis) {
  # Match integrated points back to original based on timestamps
  # This is approximate since integration may shift timestamps

  errors <- c()

  if (nrow(integrated) > 0) {
    # Sample 100 random points for efficiency
    n_sample <- min(100, nrow(integrated))
    sample_idx <- sample(1:nrow(integrated), n_sample)

    for (idx in sample_idx) {
      int_point <- integrated[idx, ]

      # Find closest original point by time
      time_diff <- abs(original_gaze$device_time_stamp - int_point$gaze_id)
      closest_idx <- which.min(time_diff)

      if (length(closest_idx) > 0 && time_diff[closest_idx] < 1000) {  # Within 1ms
        orig_point <- original_gaze[closest_idx, ]

        # Calculate Euclidean distance
        dist <- sqrt((int_point$gaze_x - orig_point$gaze_point_on_display_area_x)^2 +
                       (int_point$gaze_y - orig_point$gaze_point_on_display_area_y)^2)
        errors <- c(errors, dist)
      }
    }
  }

  if (length(errors) == 0) errors <- c(0)

  return(list(
    mean_error = mean(errors),
    sd_error = sd(errors),
    max_error = max(errors),
    accurate_pct = sum(errors < 0.01) / length(errors) * 100  # <1% of screen
  ))
}

cat("=== Calculating Performance Metrics for Experiment 1 ===\n\n")

# Initialize comprehensive metrics storage
all_metrics <- data.frame()

# Process each modality
modalities <- names(integration_results)

for (modality in modalities) {
  cat(sprintf("--- Analyzing %s modality ---\n", toupper(modality)))

  modality_results <- integration_results[[modality]]
  original_sessions <- all_sessions[[modality]]

  # Calculate metrics for each session
  for (i in 1:length(modality_results)) {
    result <- modality_results[[i]]
    original <- original_sessions[[i]]

    # 1. Latency estimation error (already calculated)
    latency_error <- abs(result$estimated_latency_ms - result$true_latency_ms)

    # 2. Assignment success rate (already calculated)
    assignment_rate <- result$matched_pct

    # 3. Processing speed
    processing_speed_per_1000 <- (result$processing_time_sec / result$n_gaze_original) * 1000

    # 4. Spatial accuracy
    # Compare assigned gaze positions with original (accounting for noise)
    spatial_errors <- calculate_spatial_accuracy(
      original$gaze_data,
      result$integrated,
      original$z_axis
    )

    # Compile metrics for this session
    session_metrics <- data.frame(
      modality = modality,
      session_id = i,

      # Latency metrics
      true_latency_ms = result$true_latency_ms,
      estimated_latency_ms = result$estimated_latency_ms,
      latency_error_ms = latency_error,
      latency_accurate = latency_error < 50,  # <50ms threshold from paper

      # Assignment metrics
      n_gaze_original = result$n_gaze_original,
      n_gaze_matched = result$n_gaze_matched,
      assignment_rate_pct = assignment_rate,

      # Processing metrics
      processing_time_sec = result$processing_time_sec,
      processing_speed_ms_per_1000 = processing_speed_per_1000,
      points_per_second = result$n_gaze_original / result$processing_time_sec,

      # Spatial accuracy metrics
      mean_spatial_error = spatial_errors$mean_error,
      sd_spatial_error = spatial_errors$sd_error,
      max_spatial_error = spatial_errors$max_error,
      spatial_accurate_pct = spatial_errors$accurate_pct
    )

    all_metrics <- rbind(all_metrics, session_metrics)
  }
}



# Save detailed metrics
write.csv(all_metrics, file.path(output_dir, "detailed_metrics.csv"), row.names = FALSE)
cat("\nSaved detailed metrics to:", file.path(output_dir, "detailed_metrics.csv"), "\n")

# Generate summary statistics by modality
cat("\n=== Performance Summary by Modality ===\n")

summary_metrics <- all_metrics %>%
  group_by(modality) %>%
  summarise(
    n_sessions = n(),

    # Latency metrics
    mean_latency_error = mean(latency_error_ms),
    sd_latency_error = sd(latency_error_ms),
    pct_latency_accurate = mean(latency_accurate) * 100,

    # Assignment metrics
    mean_assignment_rate = mean(assignment_rate_pct),
    sd_assignment_rate = sd(assignment_rate_pct),

    # Processing metrics
    mean_proc_speed_ms = mean(processing_speed_ms_per_1000),
    sd_proc_speed_ms = sd(processing_speed_ms_per_1000),
    mean_points_per_sec = mean(points_per_second),

    # Spatial metrics
    mean_spatial_error = mean(mean_spatial_error),
    sd_spatial_error = mean(sd_spatial_error),

    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 2)))

print(summary_metrics)

# Save summary
write.csv(summary_metrics, file.path(output_dir, "performance_summary.csv"), row.names = FALSE)

# Create publication-ready comparison table
cat("\n=== Creating Publication Table ===\n")

pub_table <- summary_metrics %>%
  select(
    Modality = modality,
    `Latency Error (ms)` = mean_latency_error,
    `Assignment Rate (%)` = mean_assignment_rate,
    `Processing (ms/1000pts)` = mean_proc_speed_ms,
    `Spatial Error` = mean_spatial_error
  ) %>%
  mutate(
    `Latency Error (ms)` = sprintf("%.1f ± %.1f",
                                   summary_metrics$mean_latency_error,
                                   summary_metrics$sd_latency_error),
    `Assignment Rate (%)` = sprintf("%.1f ± %.1f",
                                    summary_metrics$mean_assignment_rate,
                                    summary_metrics$sd_assignment_rate),
    `Processing (ms/1000pts)` = sprintf("%.2f", summary_metrics$mean_proc_speed_ms),
    `Spatial Error` = sprintf("%.4f", summary_metrics$mean_spatial_error)
  )

print(pub_table)
write.csv(pub_table, file.path(output_dir, "table1_modality_comparison.csv"), row.names = FALSE)

# Create final visualization
cat("\n=== Creating Results Visualization ===\n")

library(ggplot2)
library(patchwork)

# Prepare data for plotting
plot_data <- all_metrics %>%
  mutate(modality = factor(modality, levels = c("reference", "research", "webcam")))

# Create individual plots
p1 <- ggplot(plot_data, aes(x = modality, y = latency_error_ms, fill = modality)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("reference" = "#E74C3C", "research" = "#3498DB", "webcam" = "#2ECC71")) +
  labs(title = "A. Latency Estimation Error", y = "Error (ms)", x = "") +
  theme_minimal() +
  theme(legend.position = "none")

p2 <- ggplot(plot_data, aes(x = modality, y = assignment_rate_pct, fill = modality)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("reference" = "#E74C3C", "research" = "#3498DB", "webcam" = "#2ECC71")) +
  labs(title = "B. Gaze Assignment Success", y = "Success Rate (%)", x = "") +
  theme_minimal() +
  theme(legend.position = "none")

p3 <- ggplot(plot_data, aes(x = modality, y = processing_speed_ms_per_1000, fill = modality)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("reference" = "#E74C3C", "research" = "#3498DB", "webcam" = "#2ECC71")) +
  labs(title = "C. Processing Speed", y = "Time per 1000 points (ms)", x = "") +
  theme_minimal() +
  theme(legend.position = "none")

p4 <- ggplot(plot_data, aes(x = modality, y = mean_spatial_error, fill = modality)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("reference" = "#E74C3C", "research" = "#3498DB", "webcam" = "#2ECC71")) +
  labs(title = "D. Spatial Accuracy", y = "Mean Error (normalized)", x = "") +
  theme_minimal() +
  theme(legend.position = "none")

# Combine plots
combined_plot <- (p1 | p2) / (p3 | p4)

ggsave(file.path(output_dir, "figure1_performance_metrics.pdf"),
       combined_plot, width = 10, height = 8)
ggsave(file.path(output_dir, "figure1_performance_metrics.png"),
       combined_plot, width = 10, height = 8, dpi = 300)

cat("Saved performance figure to:", file.path(output_dir, "figure1_performance_metrics.pdf"), "\n")

# Statistical comparisons
cat("\n=== Statistical Comparisons ===\n")

# Compare research vs webcam on key metrics
research_data <- all_metrics[all_metrics$modality == "research", ]
webcam_data <- all_metrics[all_metrics$modality == "webcam", ]

# Latency error comparison
t_latency <- t.test(research_data$latency_error_ms, webcam_data$latency_error_ms)
cat(sprintf("Latency error - Research vs Webcam: t = %.2f, p = %.4f\n",
            t_latency$statistic, t_latency$p.value))

# Assignment rate comparison
t_assignment <- t.test(research_data$assignment_rate_pct, webcam_data$assignment_rate_pct)
cat(sprintf("Assignment rate - Research vs Webcam: t = %.2f, p = %.4f\n",
            t_assignment$statistic, t_assignment$p.value))

cat("\nAnalysis complete!\n")
