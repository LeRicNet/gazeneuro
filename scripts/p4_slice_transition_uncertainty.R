devtools::load_all()
library(tidyverse)
library(gazeneuro)
library(furrr)

plan(multisession, workers = 4)

process_case_p4 <- function(case_dir) {
  case_name <- basename(case_dir)
  gaze_file <- file.path(case_dir, "gaze_data.csv")
  axial_file <- file.path(case_dir, "axial_slice_events.csv")
  sagittal_file <- file.path(case_dir, "sagittal_slice_events.csv")

  if (!file.exists(gaze_file)) {
    return(tibble(case = case_name, error = "no gaze file"))
  }

  tryCatch({
    gaze <- read_csv(gaze_file, show_col_types = FALSE)

    # Load slice events
    axial <- if (file.exists(axial_file)) read_csv(axial_file, show_col_types = FALSE) else tibble()
    sagittal <- if (file.exists(sagittal_file)) read_csv(sagittal_file, show_col_types = FALSE) else tibble()

    slice_events <- bind_rows(axial, sagittal) %>% arrange(client_timestamp)

    if (nrow(slice_events) < 5) {
      return(tibble(case = case_name, error = "insufficient transitions"))
    }

    # Integrate
    integrated <- integrate_all_gaze_points(gaze, slice_events)
    integrated$time_sec <- integrated$time_sec.x

    # Get transition times
    transition_times <- slice_events %>%
      arrange(client_timestamp) %>%
      mutate(time_sec = (client_timestamp - min(client_timestamp)) / 1000) %>%
      pull(time_sec)

    # Distance to nearest transition
    dist_to_nearest <- function(t, transitions) {
      min(abs(t - transitions))
    }

    integrated <- integrated %>%
      mutate(dist_to_transition = map_dbl(time_sec, ~dist_to_nearest(.x, transition_times)))

    # Timing uncertainty (20ms conservative estimate)
    sigma_delta_sec <- 0.020

    integrated <- integrated %>%
      mutate(
        p_correct = pnorm(dist_to_transition / sigma_delta_sec),
        in_ambiguous_zone = dist_to_transition < 3 * sigma_delta_sec
      )

    tibble(
      case = case_name,
      n_gaze = nrow(integrated),
      n_transitions = length(transition_times),
      pct_ambiguous = mean(integrated$in_ambiguous_zone) * 100,
      mean_p_correct = mean(integrated$p_correct),
      median_dist_ms = median(integrated$dist_to_transition) * 1000,
      n_near_boundary = sum(integrated$dist_to_transition < 0.005),
      mean_p_near_boundary = if (sum(integrated$dist_to_transition < 0.005) > 0)
        mean(integrated$p_correct[integrated$dist_to_transition < 0.005]) else NA_real_,
      error = NA_character_
    )

  }, error = function(e) {
    tibble(case = case_name, error = as.character(e$message))
  })
}

case_dirs <- list.dirs('data/test_data/tc_001', recursive = FALSE, full.names = TRUE)
case_dirs <- case_dirs[!str_detect(case_dirs, "results")]

cat("Processing", length(case_dirs), "cases...\n")
results_p4 <- future_map_dfr(case_dirs, process_case_p4, .progress = TRUE)

plan(sequential)

cat("\n\nResults:\n")
print(results_p4 %>% dplyr::filter(is.na(error)) %>%
        dplyr::select(case, n_transitions, pct_ambiguous, mean_p_correct, median_dist_ms))

cat("\n\nSummary (cases with sufficient transitions):\n")
results_p4 %>%
  dplyr::filter(is.na(error)) %>%
  summarise(
    n_cases = n(),
    mean_pct_ambiguous = mean(pct_ambiguous),
    sd_pct_ambiguous = sd(pct_ambiguous),
    mean_p_correct = mean(mean_p_correct),
    mean_median_dist_ms = mean(median_dist_ms)
  ) %>%
  print()

# TODO: needs summary figure like other pX scripts in this folder. Currently this script fails on 'n_transitions doesn't exist'
