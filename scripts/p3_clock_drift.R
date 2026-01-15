library(tidyverse)
library(furrr)

plan(multisession, workers = 4)

process_case_p3 <- function(case_dir) {
  case_name <- basename(case_dir)
  gaze_file <- file.path(case_dir, "gaze_data.csv")

  if (!file.exists(gaze_file)) {
    return(tibble(case = case_name, error = "no gaze file"))
  }

  tryCatch({
    gaze <- read_csv(gaze_file, show_col_types = FALSE) %>%
      arrange(device_time_stamp) %>%
      mutate(
        device_sec = (device_time_stamp - min(device_time_stamp)) / 1e6,
        system_sec = (system_time_stamp - min(system_time_stamp)) / 1e6,
        drift_ms = (system_sec - device_sec) * 1000
      )

    # Relative to start
    gaze <- gaze %>% mutate(drift_ms = drift_ms - drift_ms[1])

    duration <- max(gaze$device_sec)
    total_drift <- gaze$drift_ms[nrow(gaze)]

    # Fit linear model
    fit <- lm(drift_ms ~ device_sec, data = gaze)
    drift_rate <- coef(fit)[2]  # ms/sec
    r_squared <- summary(fit)$r.squared

    tibble(
      case = case_name,
      n_samples = nrow(gaze),
      duration_sec = duration,
      total_drift_ms = total_drift,
      drift_rate_us_sec = drift_rate * 1000,  # Convert to µs/sec
      drift_rate_ppm = drift_rate * 1000,     # µs/sec = ppm
      r_squared = r_squared,
      error = NA_character_
    )

  }, error = function(e) {
    tibble(case = case_name, error = as.character(e$message))
  })
}

case_dirs <- list.dirs('data/test_data/tc_001', recursive = FALSE, full.names = TRUE)
case_dirs <- case_dirs[!str_detect(case_dirs, "results")]

cat("Processing", length(case_dirs), "cases...\n")
results_p3 <- future_map_dfr(case_dirs, process_case_p3, .progress = TRUE)

plan(sequential)

cat("\n\nResults:\n")
print(results_p3 %>% select(case, duration_sec, total_drift_ms, drift_rate_us_sec, r_squared))

cat("\n\nSummary:\n")
results_p3 %>%
  filter(is.na(error)) %>%
  summarise(
    n_cases = n(),
    mean_duration = mean(duration_sec),
    mean_drift_rate = mean(drift_rate_us_sec),
    sd_drift_rate = sd(drift_rate_us_sec),
    mean_r2 = mean(r_squared),
    min_r2 = min(r_squared)
  ) %>%
  print()


# Panel A: Drift rate distribution
p1 <- ggplot(results_p3, aes(x = drift_rate_us_sec)) +
  geom_histogram(bins = 10, fill = "steelblue", color = "white") +
  geom_vline(xintercept = mean(results_p3$drift_rate_us_sec), linetype = "dashed", color = "red") +
  labs(
    title = "A: Drift Rate Distribution",
    x = "Drift Rate (us/sec)",
    y = "Count"
  ) +
  theme_minimal()

# Panel B: Total drift vs duration (should be linear)
p2 <- ggplot(results_p3, aes(x = duration_sec, y = total_drift_ms)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
  labs(
    title = "B: Total Drift vs Session Duration",
    x = "Session Duration (sec)",
    y = "Total Drift (ms)"
  ) +
  theme_minimal()

# Panel C: R² values (should all be ~1)
p3 <- ggplot(results_p3, aes(x = case, y = r_squared)) +
  geom_col(fill = "darkorange") +
  geom_hline(yintercept = 0.99, linetype = "dashed", color = "red") +
  coord_cartesian(ylim = c(0.995, 1.001)) +
  labs(
    title = "C: Linearity (R-squared)",
    x = "Case",
    y = "R-squared"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

library(patchwork)
p_combined <- p1 + p2 + p3 +
  plot_annotation(
    title = "P3: Clock Drift Accumulation Across 19 Sessions",
    subtitle = sprintf("Mean drift rate = %.1f us/sec (SD = %.2f), Mean R-squared = %.4f",
                       mean(results_p3$drift_rate_us_sec), sd(results_p3$drift_rate_us_sec),
                       mean(results_p3$r_squared))
  )

print(p_combined)

ggsave("data/test_data/tc_001/results/p3_clock_drift_summary.png", p_combined, width = 12, height = 4)
cat("\nSaved to results/p3_clock_drift_summary.png\n")
