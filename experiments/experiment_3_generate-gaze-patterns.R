# 02_generate_gaze_patterns.R
# Generate synthetic gaze patterns targeting specific anatomical coordinates
# for coordinate transformation validation

library(qs)
library(tidyverse)
library(gazeneuro)

# Load synthetic NIfTI volumes from previous step
volumes <- qread("coordinate_validation/synthetic_nifti_volumes.qs")

#' Define anatomical target points in world coordinates
#' @param volume_info Information about the NIfTI volume
#' @return Data frame of target points with world coordinates
define_anatomical_targets <- function(volume_info) {
  tc <- volume_info$test_case
  affine <- volume_info$affine

  # Calculate volume center in voxel coordinates
  center_voxel <- (tc$dims - 1) / 2

  # Convert to world coordinates
  center_world <- affine %*% c(center_voxel, 1)
  center_world <- center_world[1:3]

  # Define target points relative to center
  targets <- data.frame(
    target_id = 1:9,
    description = c(
      "center",
      "anterior_20mm", "posterior_20mm",
      "left_20mm", "right_20mm",
      "superior_20mm", "inferior_20mm",
      "corner_anterior_left_superior",
      "corner_posterior_right_inferior"
    ),
    world_x = c(
      center_world[1],
      center_world[1] + 20, center_world[1] - 20,
      center_world[1] - 20, center_world[1] + 20,
      center_world[1], center_world[1],
      center_world[1] + 40, center_world[1] - 40
    ),
    world_y = c(
      center_world[2],
      center_world[2], center_world[2],
      center_world[2], center_world[2],
      center_world[2] + 20, center_world[2] - 20,
      center_world[2] + 40, center_world[2] - 40
    ),
    world_z = c(
      center_world[3],
      center_world[3], center_world[3],
      center_world[3], center_world[3],
      center_world[3], center_world[3],
      center_world[3] + 30, center_world[3] - 30
    )
  )

  return(targets)
}

#' Convert world coordinates to voxel and normalized screen coordinates
#' @param world_coords World coordinates (x, y, z)
#' @param affine Affine transformation matrix
#' @param dims Volume dimensions
#' @param plane Viewing plane (AXIAL, SAGITTAL, CORONAL)
#' @param slice_idx Slice index for the viewing plane
#' @return List with voxel and screen coordinates
world_to_screen <- function(world_coords, affine, dims, plane, slice_idx) {
  # Convert world to voxel
  inv_affine <- solve(affine)
  voxel <- inv_affine %*% c(world_coords, 1)
  voxel <- voxel[1:3]

  # Check if voxel is within bounds
  if (any(voxel < 0) || any(voxel >= dims)) {
    return(list(valid = FALSE))
  }

  # Extract 2D coordinates based on plane
  if (plane == "AXIAL") {
    # Check if target is on this slice
    if (abs(voxel[3] - slice_idx) > 0.5) {
      return(list(valid = FALSE))
    }
    screen_x <- voxel[1] / (dims[1] - 1)
    screen_y <- voxel[2] / (dims[2] - 1)
  } else if (plane == "SAGITTAL") {
    if (abs(voxel[1] - slice_idx) > 0.5) {
      return(list(valid = FALSE))
    }
    screen_x <- voxel[2] / (dims[2] - 1)
    screen_y <- voxel[3] / (dims[3] - 1)
  } else { # CORONAL
    if (abs(voxel[2] - slice_idx) > 0.5) {
      return(list(valid = FALSE))
    }
    screen_x <- voxel[1] / (dims[1] - 1)
    screen_y <- voxel[3] / (dims[3] - 1)
  }

  return(list(
    valid = TRUE,
    voxel = voxel,
    screen_x = screen_x,
    screen_y = screen_y
  ))
}

#' Generate fixation pattern around a target
#' @param target_x Target x coordinate (normalized)
#' @param target_y Target y coordinate (normalized)
#' @param duration Duration in seconds
#' @param sampling_rate Sampling rate in Hz
#' @param noise_sd Standard deviation of noise in normalized units
#' @return Data frame of gaze points
generate_fixation <- function(target_x, target_y, duration,
                              sampling_rate = 120, noise_sd = 0.002) {
  n_samples <- round(duration * sampling_rate)

  # Generate gaze points with Gaussian noise
  gaze <- data.frame(
    time_offset = seq(0, duration, length.out = n_samples),
    gaze_x = pmax(0, pmin(1, target_x + rnorm(n_samples, 0, noise_sd))),
    gaze_y = pmax(0, pmin(1, target_y + rnorm(n_samples, 0, noise_sd))),
    pattern_type = "fixation"
  )

  return(gaze)
}

#' Generate saccade between two points
#' @param start_x Starting x coordinate
#' @param start_y Starting y coordinate
#' @param end_x Ending x coordinate
#' @param end_y Ending y coordinate
#' @param duration Saccade duration in seconds
#' @param sampling_rate Sampling rate in Hz
#' @return Data frame of gaze points
generate_saccade <- function(start_x, start_y, end_x, end_y,
                             duration = 0.03, sampling_rate = 120) {
  n_samples <- round(duration * sampling_rate)

  # Sigmoid function for smooth transition
  t <- seq(0, 1, length.out = n_samples)
  sigmoid <- 1 / (1 + exp(-10 * (t - 0.5)))

  gaze <- data.frame(
    time_offset = seq(0, duration, length.out = n_samples),
    gaze_x = start_x + (end_x - start_x) * sigmoid,
    gaze_y = start_y + (end_y - start_y) * sigmoid,
    pattern_type = "saccade"
  )

  return(gaze)
}

#' Generate scanning pattern across slice
#' @param dims Slice dimensions
#' @param duration Total duration
#' @param sampling_rate Sampling rate in Hz
#' @return Data frame of gaze points
generate_scanning_pattern <- function(dims, duration = 5, sampling_rate = 120) {
  # Create a raster scan pattern
  n_rows <- 5
  n_cols <- 5

  x_points <- seq(0.1, 0.9, length.out = n_cols)
  y_points <- seq(0.1, 0.9, length.out = n_rows)

  # Create scan path
  scan_path <- expand.grid(x = x_points, y = y_points)
  scan_path <- scan_path[order(scan_path$y, scan_path$x), ]

  # Generate gaze pattern
  gaze_pattern <- data.frame()
  current_time <- 0
  dwell_time <- duration / nrow(scan_path)

  for (i in 1:nrow(scan_path)) {
    # Add fixation at this point
    fixation <- generate_fixation(
      scan_path$x[i], scan_path$y[i],
      dwell_time * 0.8, sampling_rate
    )
    fixation$time_offset <- fixation$time_offset + current_time
    gaze_pattern <- rbind(gaze_pattern, fixation)
    current_time <- max(fixation$time_offset)

    # Add saccade to next point (if not last)
    if (i < nrow(scan_path)) {
      saccade <- generate_saccade(
        scan_path$x[i], scan_path$y[i],
        scan_path$x[i + 1], scan_path$y[i + 1],
        duration = dwell_time * 0.2, sampling_rate
      )
      saccade$time_offset <- saccade$time_offset + current_time
      saccade$pattern_type <- "scan_saccade"
      gaze_pattern <- rbind(gaze_pattern, saccade)
      current_time <- max(saccade$time_offset)
    }
  }

  return(gaze_pattern)
}

# Generate gaze patterns for each test volume
all_gaze_patterns <- list()
pattern_id <- 1

for (volume_name in names(volumes)) {
  volume_info <- volumes[[volume_name]]
  tc <- volume_info$test_case

  message(sprintf("Generating gaze patterns for: %s", volume_name))

  # Define anatomical targets
  targets <- define_anatomical_targets(volume_info)

  # For each viewing plane
  for (plane in c("AXIAL", "SAGITTAL", "CORONAL")) {
    # Select appropriate slice indices
    if (plane == "AXIAL") {
      slice_indices <- c(
        round(tc$dims[3] * 0.3),  # Lower slice
        round(tc$dims[3] * 0.5),  # Middle slice
        round(tc$dims[3] * 0.7)   # Upper slice
      )
    } else if (plane == "SAGITTAL") {
      slice_indices <- c(
        round(tc$dims[1] * 0.3),
        round(tc$dims[1] * 0.5),
        round(tc$dims[1] * 0.7)
      )
    } else {  # CORONAL
      slice_indices <- c(
        round(tc$dims[2] * 0.3),
        round(tc$dims[2] * 0.5),
        round(tc$dims[2] * 0.7)
      )
    }

    for (slice_idx in slice_indices) {
      # Pattern 1: Target fixations
      target_gaze <- data.frame()
      current_time <- 0
      valid_targets <- 0

      for (i in 1:nrow(targets)) {
        target_screen <- world_to_screen(
          as.numeric(targets[i, c("world_x", "world_y", "world_z")]),
          volume_info$affine,
          tc$dims,
          plane,
          slice_idx
        )

        if (target_screen$valid) {
          # Generate 1-second fixation
          fixation <- generate_fixation(
            target_screen$screen_x,
            target_screen$screen_y,
            duration = 1.0
          )
          fixation$time_offset <- fixation$time_offset + current_time
          fixation$target_id <- targets$target_id[i]
          fixation$target_desc <- targets$description[i]
          target_gaze <- rbind(target_gaze, fixation)
          current_time <- max(fixation$time_offset) + 0.2  # 200ms gap
          valid_targets <- valid_targets + 1
        }
      }

      if (valid_targets > 0) {
        target_gaze$pattern_id <- pattern_id
        target_gaze$volume_name <- volume_name
        target_gaze$plane <- plane
        target_gaze$slice_idx <- slice_idx
        target_gaze$pattern_name <- "target_fixations"

        all_gaze_patterns[[pattern_id]] <- target_gaze
        pattern_id <- pattern_id + 1
      }

      # Pattern 2: Scanning pattern
      scan_gaze <- generate_scanning_pattern(tc$dims, duration = 5)
      scan_gaze$pattern_id <- pattern_id
      scan_gaze$volume_name <- volume_name
      scan_gaze$plane <- plane
      scan_gaze$slice_idx <- slice_idx
      scan_gaze$pattern_name <- "scanning"
      scan_gaze$target_id <- NA
      scan_gaze$target_desc <- NA

      all_gaze_patterns[[pattern_id]] <- scan_gaze
      pattern_id <- pattern_id + 1
    }
  }
}

# Combine all patterns
combined_patterns <- bind_rows(all_gaze_patterns)

# Add microsecond timestamps (simulating eye tracker output)
combined_patterns <- combined_patterns %>%
  group_by(pattern_id) %>%
  mutate(
    device_time_stamp = round((time_offset + runif(1, 0, 1000)) * 1e6),
    gaze_point_on_display_area_x = gaze_x,
    gaze_point_on_display_area_y = gaze_y
  ) %>%
  ungroup()

# Save results
qsave(combined_patterns, "coordinate_validation/synthetic_gaze_patterns.qs")

# Create summary
pattern_summary <- combined_patterns %>%
  group_by(pattern_id, volume_name, plane, slice_idx, pattern_name) %>%
  summarise(
    n_points = n(),
    duration_sec = max(time_offset),
    n_fixations = sum(pattern_type == "fixation"),
    n_targets = n_distinct(target_id[!is.na(target_id)]),
    .groups = "drop"
  )

write_csv(pattern_summary, "coordinate_validation/gaze_patterns_summary.csv")

# Print summary
cat("\nGenerated", n_distinct(combined_patterns$pattern_id), "gaze patterns\n")
cat("Total gaze points:", nrow(combined_patterns), "\n")
cat("Patterns per volume:", n_distinct(combined_patterns$pattern_id) / length(volumes), "\n")

cat("\nPattern distribution:\n")
print(table(combined_patterns$pattern_name))

cat("\nData saved to coordinate_validation/synthetic_gaze_patterns.qs\n")
