library(tidyverse)
library(furrr)

# Set up parallel processing
plan(multisession, workers = 4)

# Get complete cases
case_dirs <- list.dirs('data/test_data/tc_001', recursive = FALSE, full.names = TRUE)
case_dirs <- case_dirs[!str_detect(case_dirs, "results")]

# Function to process one case
process_case_p1 <- function(case_dir) {
  case_name <- basename(case_dir)
  gaze_file <- file.path(case_dir, "gaze_data.csv")

  if (!file.exists(gaze_file)) {
    return(tibble(case = case_name, error = "no gaze file"))
  }

  tryCatch({
    gaze <- read_csv(gaze_file, show_col_types = FALSE) %>%
      arrange(device_time_stamp) %>%
      filter(!is.na(gaze_point_on_display_area_x), !is.na(gaze_point_on_display_area_y)) %>%
      mutate(time_sec = (device_time_stamp - min(device_time_stamp)) / 1e6)

    if (nrow(gaze) < 100) {
      return(tibble(case = case_name, error = "insufficient data"))
    }

    k <- 5
    gaze <- gaze %>%
      mutate(
        x_smooth = zoo::rollmean(gaze_point_on_display_area_x, k = k, fill = NA, align = "center"),
        y_smooth = zoo::rollmean(gaze_point_on_display_area_y, k = k, fill = NA, align = "center"),
        x_resid = gaze_point_on_display_area_x - x_smooth,
        y_resid = gaze_point_on_display_area_y - y_smooth,
        resid_magnitude = sqrt(x_resid^2 + y_resid^2),
        dx = (lead(x_smooth) - lag(x_smooth)) / (lead(time_sec) - lag(time_sec)),
        dy = (lead(y_smooth) - lag(y_smooth)) / (lead(time_sec) - lag(time_sec)),
        velocity = sqrt(dx^2 + dy^2)
      ) %>%
      filter(!is.na(velocity), !is.na(resid_magnitude), velocity > 0, resid_magnitude > 0)

    # Compute binned statistics
    velocity_bins <- gaze %>%
      mutate(velocity_bin = cut(velocity, breaks = quantile(velocity, probs = seq(0, 1, 0.1)), include.lowest = TRUE)) %>%
      group_by(velocity_bin) %>%
      summarise(
        n = n(),
        mean_velocity = mean(velocity),
        mean_residual = mean(resid_magnitude),
        .groups = "drop"
      ) %>%
      filter(!is.na(velocity_bin))

    # Correlation on binned data
    r_binned <- cor(velocity_bins$mean_velocity, velocity_bins$mean_residual)

    # Log-log slope
    fit <- lm(log10(mean_residual) ~ log10(mean_velocity), data = velocity_bins)
    slope <- coef(fit)[2]
    r_squared <- summary(fit)$r.squared

    tibble(
      case = case_name,
      n_samples = nrow(gaze),
      r_binned = r_binned,
      log_log_slope = slope,
      log_log_r2 = r_squared,
      error = NA_character_
    )

  }, error = function(e) {
    tibble(case = case_name, error = as.character(e$message))
  })
}

# Run in parallel
cat("Processing", length(case_dirs), "cases...\n")
results_p1 <- future_map_dfr(case_dirs, process_case_p1, .progress = TRUE)

plan(sequential)

# Summary
cat("\n\nResults Summary:\n")
print(results_p1 %>% select(-error))

cat("\nAggregate statistics:\n")
results_p1 %>%
  filter(is.na(error)) %>%
  summarise(
    n_cases = n(),
    mean_r = mean(r_binned),
    sd_r = sd(r_binned),
    mean_slope = mean(log_log_slope),
    sd_slope = sd(log_log_slope),
    mean_r2 = mean(log_log_r2)
  ) %>%
  print()

# Panel A: Distribution of correlations
p1 <- ggplot(results_p1, aes(x = r_binned)) +
  geom_histogram(bins = 10, fill = "steelblue", color = "white") +
  geom_vline(xintercept = mean(results_p1$r_binned), linetype = "dashed", color = "red") +
  labs(
    title = "A: Correlation Distribution",
    x = "Pearson r (velocity vs residual)",
    y = "Count"
  ) +
  theme_minimal()

# Panel B: Distribution of slopes
p2 <- ggplot(results_p1, aes(x = log_log_slope)) +
  geom_histogram(bins = 10, fill = "darkorange", color = "white") +
  geom_vline(xintercept = mean(results_p1$log_log_slope), linetype = "dashed", color = "red") +
  labs(
    title = "B: Power Law Exponent",
    subtitle = expression(residual %prop% velocity^beta),
    x = expression(beta ~ "(log-log slope)"),
    y = "Count"
  ) +
  theme_minimal()

# Panel C: r vs slope scatter
p3 <- ggplot(results_p1, aes(x = log_log_slope, y = r_binned)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "gray50") +
  labs(
    title = "C: Slope vs Correlation",
    x = expression(beta),
    y = "r"
  ) +
  theme_minimal()

library(patchwork)
p_combined <- p1 + p2 + p3 +
  plot_annotation(
    title = "P1: Velocity-Spatial Coupling Across 19 Sessions",
    subtitle = sprintf("Mean r = %.3f (SD = %.3f), Mean β = %.3f (SD = %.3f)",
                       mean(results_p1$r_binned), sd(results_p1$r_binned),
                       mean(results_p1$log_log_slope), sd(results_p1$log_log_slope))
  )

print(p_combined)

ggsave("data/test_data/tc_001/results/p1_velocity_coupling_summary.png",
       p_combined, width = 12, height = 4)
cat("\nSaved to results/p1_velocity_coupling_summary.png\n")
