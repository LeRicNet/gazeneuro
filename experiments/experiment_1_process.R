# Process Experiment 1 Data with gazeneuro
# Runs integrate_all_gaze_points() on all 300 synthetic sessions

library(tidyverse)
library(gazeneuro)
library(tictoc)  # For timing

# Load the generated data
output_dir <- "experiment1_data"
all_sessions <- readRDS(file.path(output_dir, "all_sessions.rds"))

cat("=== Processing Synthetic Data with gazeneuro ===\n")
cat(sprintf("Total sessions to process: %d\n\n",
            sum(sapply(all_sessions, length))))

# Initialize results storage
integration_results <- list()
processing_log <- data.frame()

# Process each modality
modalities <- names(all_sessions)

for (modality in modalities) {
  cat(sprintf("\n--- Processing %s modality ---\n", toupper(modality)))

  sessions <- all_sessions[[modality]]
  n_sessions <- length(sessions)

  # Storage for this modality
  modality_results <- list()

  # Process each session
  tic()

  for (i in 1:n_sessions) {
    if (i %% 10 == 0) {
      cat(sprintf("Processing session %d/%d\n", i, n_sessions))
    }

    session <- sessions[[i]]

    # Measure processing time
    proc_start <- Sys.time()

    # Run gazeneuro integration
    # Capture messages to extract estimated latency
    msg_output <- capture.output({
      integrated <- integrate_all_gaze_points(
        gaze_data = session$gaze_data,
        z_axis = session$z_axis,
        latency_correction = NULL  # Let it auto-detect
      )
    }, type = "message")

    proc_end <- Sys.time()
    proc_time <- as.numeric(difftime(proc_end, proc_start, units = "secs"))

    # Extract estimated latency from messages
    latency_line <- grep("Estimated latency:", msg_output, value = TRUE)
    estimated_latency <- NA
    if (length(latency_line) > 0) {
      estimated_latency <- as.numeric(gsub(".*Estimated latency: ([0-9.-]+) sec.*", "\\1", latency_line[1])) * 1000
    }

    # Extract matched percentage
    matched_line <- grep("Gaze points matched:", msg_output, value = TRUE)
    matched_pct <- NA
    if (length(matched_line) > 0) {
      matched_pct <- as.numeric(gsub(".*\\(([0-9.]+)%\\).*", "\\1", matched_line[1]))
    }

    # Store results
    modality_results[[i]] <- list(
      session_id = i,
      modality = modality,
      integrated = integrated,
      true_latency_ms = session$metadata$latency_seconds,
      estimated_latency_ms = estimated_latency,
      n_gaze_original = nrow(session$gaze_data),
      n_gaze_matched = nrow(integrated),
      matched_pct = matched_pct,
      processing_time_sec = proc_time,
      n_slices = nrow(session$z_axis)
    )

    # Log this session
    processing_log <- rbind(processing_log, data.frame(
      modality = modality,
      session_id = i,
      true_latency_ms = session$metadata$latency_seconds,
      estimated_latency_ms = estimated_latency,
      latency_error_ms = abs(estimated_latency - session$metadata$latency_seconds),
      n_gaze_original = nrow(session$gaze_data),
      n_gaze_matched = nrow(integrated),
      matched_pct = matched_pct,
      processing_time_sec = proc_time,
      points_per_sec = nrow(session$gaze_data) / proc_time
    ))
  }

  elapsed <- toc(quiet = TRUE)

  # Store modality results
  integration_results[[modality]] <- modality_results

  # Print summary for this modality
  modality_log <- processing_log[processing_log$modality == modality, ]

  cat(sprintf("\nProcessed %d sessions in %.1f seconds\n",
              n_sessions, elapsed$toc - elapsed$tic))
  cat(sprintf("Mean latency error: %.2f ms (SD: %.2f)\n",
              mean(modality_log$latency_error_ms, na.rm = TRUE),
              sd(modality_log$latency_error_ms, na.rm = TRUE)))
  cat(sprintf("Mean matched percentage: %.1f%% (SD: %.1f)\n",
              mean(modality_log$matched_pct, na.rm = TRUE),
              sd(modality_log$matched_pct, na.rm = TRUE)))
  cat(sprintf("Mean processing speed: %.0f points/sec\n",
              mean(modality_log$points_per_sec, na.rm = TRUE)))
}

# Save integration results
cat("\n=== Saving Integration Results ===\n")

# Save full results
saveRDS(integration_results, file.path(output_dir, "integration_results.rds"))
cat("Saved integration results to:", file.path(output_dir, "integration_results.rds"), "\n")

# Save processing log
write.csv(processing_log, file.path(output_dir, "processing_log.csv"), row.names = FALSE)
cat("Saved processing log to:", file.path(output_dir, "processing_log.csv"), "\n")

# Create summary statistics by modality
summary_stats <- processing_log %>%
  group_by(modality) %>%
  summarise(
    n_sessions = n(),
    mean_latency_error = mean(latency_error_ms, na.rm = TRUE),
    sd_latency_error = sd(latency_error_ms, na.rm = TRUE),
    mean_matched_pct = mean(matched_pct, na.rm = TRUE),
    sd_matched_pct = sd(matched_pct, na.rm = TRUE),
    mean_proc_time = mean(processing_time_sec, na.rm = TRUE),
    total_proc_time = sum(processing_time_sec),
    mean_points_per_sec = mean(points_per_sec, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== Processing Summary by Modality ===\n")
print(summary_stats)

# Save summary
write.csv(summary_stats, file.path(output_dir, "modality_summary.csv"), row.names = FALSE)

# Create diagnostic plots
cat("\n=== Creating Diagnostic Plots ===\n")

pdf(file.path(output_dir, "processing_diagnostics.pdf"), width = 12, height = 10)

par(mfrow = c(2, 2))

# Plot 1: Latency estimation error by modality
boxplot(latency_error_ms ~ modality, data = processing_log,
        main = "Latency Estimation Error by Modality",
        ylab = "Error (ms)",
        col = c("skyblue", "lightgreen", "lightcoral"))
abline(h = 50, lty = 2, col = "red")  # 50ms threshold from paper

# Plot 2: Match percentage by modality
boxplot(matched_pct ~ modality, data = processing_log,
        main = "Gaze Point Match Percentage",
        ylab = "Matched (%)",
        col = c("skyblue", "lightgreen", "lightcoral"))

# Plot 3: Scatter plot of true vs estimated latency
plot(processing_log$true_latency_ms, processing_log$estimated_latency_ms,
     col = factor(processing_log$modality),
     pch = 19, cex = 0.6,
     main = "True vs Estimated Latency",
     xlab = "True Latency (ms)", ylab = "Estimated Latency (ms)")
abline(0, 1, lty = 2)
legend("topleft", legend = unique(processing_log$modality),
       col = 1:3, pch = 19, cex = 0.8)

# Plot 4: Processing speed by modality
boxplot(points_per_sec ~ modality, data = processing_log,
        main = "Processing Speed",
        ylab = "Points per Second",
        col = c("skyblue", "lightgreen", "lightcoral"))

dev.off()

cat("Saved diagnostic plots to:", file.path(output_dir, "processing_diagnostics.pdf"), "\n")

# Quick check: any sessions with poor matching?
poor_matches <- processing_log[processing_log$matched_pct < 90, ]
if (nrow(poor_matches) > 0) {
  cat(sprintf("\nWARNING: %d sessions had <90%% match rate\n", nrow(poor_matches)))
  print(table(poor_matches$modality))
} else {
  cat("\nAll sessions achieved >90% match rate\n")
}

cat("\nProcessing complete!\n")
