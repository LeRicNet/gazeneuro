devtools::load_all()
library(tidyverse)
library(gazeneuro)

setwd("~/interfaces/gazeneuro/data/test_data/tc_001")

# Check actual content of slice event files
failure_cases <- sprintf("case_%02d", 2:9)
planes <- c("axial", "sagittal")

for (case in failure_cases) {
  for (plane in planes) {
    file_path <- file.path(case, sprintf("%s_slice_events.csv", plane))

    if (file.exists(file_path)) {
      data <- read_csv(file_path)
      cat(sprintf("%s %s: %d rows\n", case, plane, nrow(data)))

      # Print full content if only 1-2 rows
      if (nrow(data) <= 2) {
        print(data)
      }
    }
  }
}

# Analyze gaze recording patterns for failure cases
for (case in failure_cases) {
  gaze_path <- file.path(case, "gaze_data.csv")
  gaze <- read_csv(gaze_path)

  # Convert timestamps to seconds
  gaze_time <- (gaze$device_time_stamp - min(gaze$device_time_stamp)) / 1e6

  cat(sprintf("\n%s:\n", case))
  cat(sprintf("  Duration: %.2f sec\n", max(gaze_time)))
  cat(sprintf("  N points: %d\n", nrow(gaze)))
  cat(sprintf("  Sampling rate: %.1f Hz\n", nrow(gaze) / max(gaze_time)))

  # Check for temporal gaps
  time_diffs <- diff(gaze_time)
  large_gaps <- which(time_diffs > 1.0)
  if (length(large_gaps) > 0) {
    cat(sprintf("  Large gaps (>1s): %d instances\n", length(large_gaps)))
  }
}

# Compare characteristics of success cases vs failure cases
success_cases <- c("case_10", "case_11", "case_12", "case_14", "case_15")
failure_cases <- sprintf("case_%02d", 2:9)

# Analyze viewing software logs if available
# Look for differences in:
# - Session initialization sequences
# - Software version or configuration
# - User interaction patterns
# - Hardware status indicators
