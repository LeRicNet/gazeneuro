devtools::load_all()
# Load required packages
library(gazeneuro)
library(tidyverse)

# Set working directory
setwd("~/interfaces/gazeneuro/data/test_data/tc_001")

# Create results directory
dir.create("results", showWarnings = FALSE)

# Initialize results storage (now per case and plane)
latency_results <- data.frame(
  case_id = character(),
  plane = character(),
  gaze_duration = numeric(),
  z_duration = numeric(),
  estimated_latency = numeric(),
  n_gaze_points = integer(),
  n_z_events = integer(),
  match_rate = numeric(),
  stringsAsFactors = FALSE
)

# Consider these possible planes, but detect availability per case
candidate_planes <- c("axial", "sagittal", "coronal", "transversal")

# Loop through 19 cases and only the planes that exist in each case
for (case_num in 2:19) {
  cat("\n=== Processing Case", case_num, "===\n")

  # Define paths common to the case
  case_dir <- sprintf("case_%02d", case_num)
  gaze_path <- file.path(case_dir, "gaze_data.csv")

  # Load gaze once per case
  gaze_data <- read_csv(gaze_path)

  # Detect which planes exist in this case
  available_planes <- candidate_planes[c(
    file.exists(file.path(case_dir, sprintf("%s.nii.gz", candidate_planes))) &
      file.exists(file.path(case_dir, sprintf("%s_slice_events.csv", candidate_planes)))
  )]

  if (length(available_planes) == 0) {
    warning(sprintf("No planes found for %s; skipping.", case_dir))
    next
  }

  for (pl in available_planes) {
    cat("  - Plane:", pl, "\n")

    # Define plane-specific paths
    nifti_path <- file.path(case_dir, sprintf("%s.nii.gz", pl))
    slice_path <- file.path(case_dir, sprintf("%s_slice_events.csv", pl))

    # Load data
    nifti_data <- preload_nifti_data(nifti_path)
    z_axis <- read_csv(slice_path)

    # Integrate
    integrated <- integrate_all_gaze_points(gaze_data, z_axis)

    # Calculate metrics
    gaze_duration <- (max(gaze_data$device_time_stamp) - min(gaze_data$device_time_stamp)) / 1e6
    z_duration <- (max(z_axis$client_timestamp) - min(z_axis$client_timestamp)) / 1000
    latency <- z_duration - gaze_duration
    match_rate <- nrow(integrated) / nrow(gaze_data)

    # Store results
    latency_results <- rbind(latency_results, data.frame(
      case_id = sprintf("case_%02d", case_num),
      plane = pl,
      gaze_duration = gaze_duration,
      z_duration = z_duration,
      estimated_latency = latency,
      n_gaze_points = nrow(gaze_data),
      n_z_events = nrow(z_axis),
      match_rate = match_rate
    ))

    # Ensure results dir exists
    dir.create("results", showWarnings = FALSE)

    # Quality check visualization
    png(sprintf("results/case_%02d_%s_quality_check.png", case_num, pl),
        width = 1200, height = 800)
    check_integration_quality(gaze_data, z_axis, integrated)
    dev.off()

    # Example visualization: pick a representative slice index if needed
    png(sprintf("results/case_%02d_%s_slice_12.png", case_num, pl),
        width = 800, height = 800)
    plot_slice_with_all_gaze(nifti_data, integrated, slice_num = 12)
    dev.off()
  }
}

# Save results (now includes plane column)
write_csv(latency_results, "results/latency_analysis.csv")

# Summary statistics
summary(latency_results)

# Latency consistency
cat("\nLatency Statistics:\n")
cat("Mean:", mean(latency_results$estimated_latency), "seconds\n")
cat("SD:", sd(latency_results$estimated_latency), "seconds\n")
cat("Range:", range(latency_results$estimated_latency), "seconds\n")

# Match rate analysis
cat("\nMatch Rate Statistics:\n")
cat("Mean:", mean(latency_results$match_rate) * 100, "%\n")
cat("Min:", min(latency_results$match_rate) * 100, "%\n")

# Visualize latency across cases
png("results/latency_distribution.png", width = 1000, height = 600)
par(mfrow = c(1, 2))
barplot(latency_results$estimated_latency,
        names.arg = latency_results$case_id,
        main = "Estimated Latency by Case",
        ylab = "Latency (seconds)",
        las = 2, cex.names = 0.8)
abline(h = mean(latency_results$estimated_latency), col = "red", lwd = 2)

barplot(latency_results$match_rate * 100,
        names.arg = latency_results$case_id,
        main = "Match Rate by Case",
        ylab = "Match Rate (%)",
        las = 2, cex.names = 0.8)
abline(h = 95, col = "green", lwd = 2, lty = 2)
dev.off()
