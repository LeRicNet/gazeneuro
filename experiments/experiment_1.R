# Generate Full Dataset for Experiment 1
# Creates 100 sessions × 3 modalities = 300 synthetic sessions

library(tidyverse)
library(gazeneuro)

# Source the generator functions (assuming they're loaded)
# source("synthetic_data_generator.R")

# Set up output directory
output_dir <- "experiment1_data"
dir.create(output_dir, showWarnings = FALSE)

# Initialize results storage
all_sessions <- list()
generation_log <- data.frame()

# Generate data for each modality
modalities <- c("research", "webcam", "reference")
n_sessions_per_modality <- 100
base_seed <- 42  # For reproducibility

cat("=== Generating Synthetic Data for Experiment 1 ===\n")
cat(sprintf("Total sessions to generate: %d\n\n",
            length(modalities) * n_sessions_per_modality))

for (modality in modalities) {
  cat(sprintf("\n--- Generating %s modality ---\n", toupper(modality)))

  # Generate sessions for this modality
  start_time <- Sys.time()
  sessions <- generate_modality_sessions(modality,
                                         n_sessions = n_sessions_per_modality,
                                         base_seed = base_seed)
  end_time <- Sys.time()

  # Store sessions
  all_sessions[[modality]] <- sessions

  # Calculate summary statistics
  validation_results <- map_df(sessions, validate_synthetic_data)

  # Log generation details
  generation_log <- rbind(generation_log, data.frame(
    modality = modality,
    n_sessions = n_sessions_per_modality,
    generation_time = as.numeric(difftime(end_time, start_time, units = "secs")),
    mean_samples = mean(validation_results$n_samples),
    sd_samples = sd(validation_results$n_samples),
    mean_rate_error = mean(validation_results$rate_error),
    total_gaze_points = sum(validation_results$n_samples)
  ))

  # Print summary
  cat(sprintf("Generated %d sessions in %.1f seconds\n",
              n_sessions_per_modality,
              as.numeric(difftime(end_time, start_time, units = "secs"))))
  cat(sprintf("Mean samples per session: %.0f (SD: %.0f)\n",
              mean(validation_results$n_samples),
              sd(validation_results$n_samples)))
  cat(sprintf("Mean sampling rate error: %.2f Hz\n",
              mean(validation_results$rate_error)))

  # Update base seed for next modality
  base_seed <- base_seed + 1000
}

# Save the generated data
cat("\n=== Saving Generated Data ===\n")

# Save as RDS for easy loading
saveRDS(all_sessions, file.path(output_dir, "all_sessions.rds"))
cat("Saved all sessions to:", file.path(output_dir, "all_sessions.rds"), "\n")

# Save generation log
write.csv(generation_log, file.path(output_dir, "generation_log.csv"), row.names = FALSE)
cat("Saved generation log to:", file.path(output_dir, "generation_log.csv"), "\n")

# Create summary statistics
cat("\n=== Generation Summary ===\n")
print(generation_log)

# Calculate total statistics
total_gaze_points <- sum(generation_log$total_gaze_points)
total_time <- sum(generation_log$generation_time)

cat(sprintf("\nTotal gaze points generated: %s\n",
            format(total_gaze_points, big.mark = ",")))
cat(sprintf("Total generation time: %.1f seconds\n", total_time))
cat(sprintf("Average time per session: %.2f seconds\n",
            total_time / (n_sessions_per_modality * length(modalities))))

# Create validation plots
cat("\n=== Creating Validation Plots ===\n")

# Combine validation data for all modalities
all_validation <- list()
for (modality in modalities) {
  val_data <- map_df(all_sessions[[modality]], validate_synthetic_data)
  val_data$modality <- modality
  all_validation[[modality]] <- val_data
}
validation_df <- bind_rows(all_validation)

# Create validation plots
pdf(file.path(output_dir, "validation_plots.pdf"), width = 10, height = 8)

# Plot 1: Sampling rates by modality
par(mfrow = c(2, 2))

boxplot(actual_rate ~ modality, data = validation_df,
        main = "Actual Sampling Rates by Modality",
        ylab = "Sampling Rate (Hz)",
        col = c("skyblue", "lightgreen", "lightcoral"))

# Plot 2: Number of samples distribution
boxplot(n_samples ~ modality, data = validation_df,
        main = "Number of Samples per Session",
        ylab = "Sample Count",
        col = c("skyblue", "lightgreen", "lightcoral"))

# Plot 3: Spatial spread
plot(validation_df$spatial_spread_x, validation_df$spatial_spread_y,
     col = factor(validation_df$modality),
     pch = 19, cex = 0.5,
     main = "Spatial Spread by Modality",
     xlab = "X Spread (SD)", ylab = "Y Spread (SD)")
legend("topright", legend = modalities,
       col = 1:3, pch = 19, cex = 0.8)

# Plot 4: Rate errors
hist(validation_df$rate_error, breaks = 30,
     main = "Sampling Rate Errors",
     xlab = "Rate Error (Hz)",
     col = "gray80")

dev.off()

cat("Saved validation plots to:", file.path(output_dir, "validation_plots.pdf"), "\n")

# Save a sample session for each modality as CSV for inspection
cat("\n=== Saving Sample Sessions ===\n")
for (modality in modalities) {
  sample_session <- all_sessions[[modality]][[1]]

  # Save gaze data
  write.csv(sample_session$gaze_data,
            file.path(output_dir, sprintf("sample_%s_gaze.csv", modality)),
            row.names = FALSE)

  # Save z-axis data
  write.csv(sample_session$z_axis,
            file.path(output_dir, sprintf("sample_%s_zaxis.csv", modality)),
            row.names = FALSE)
}

cat("\nData generation complete!\n")
cat("All files saved to:", output_dir, "\n")
