# 03_validate_coordinate_transformations.R
# Validate coordinate transformations through the gazeneuro pipeline
# Tests accuracy of screen → voxel → world transformations

library(qs)
library(tidyverse)
library(gazeneuro)
library(RNifti)

# Load test data
volumes <- qread("coordinate_validation/synthetic_nifti_volumes.qs")
gaze_patterns <- qread("coordinate_validation/synthetic_gaze_patterns.qs")

#' Perform round-trip transformation test
#' @param gaze_point Single gaze point with screen coordinates
#' @param nifti_data Preloaded NIfTI data (from gazeneuro)
#' @param plane Viewing plane
#' @param slice_idx Slice index
#' @return Data frame with transformation results and errors
test_round_trip <- function(gaze_point, nifti_data, plane, slice_idx) {
  nvimg <- nifti_data$nvimage
  dims <- nifti_data$dims

  # Step 1: Screen to voxel coordinates
  # Convert normalized screen coords to voxel based on plane
  if (plane == "AXIAL") {
    voxel_x <- gaze_point$gaze_x * (dims[1] - 1)
    voxel_y <- gaze_point$gaze_y * (dims[2] - 1)
    voxel_z <- slice_idx
  } else if (plane == "SAGITTAL") {
    voxel_x <- slice_idx
    voxel_y <- gaze_point$gaze_x * (dims[2] - 1)
    voxel_z <- gaze_point$gaze_y * (dims[3] - 1)
  } else {  # CORONAL
    voxel_x <- gaze_point$gaze_x * (dims[1] - 1)
    voxel_y <- slice_idx
    voxel_z <- gaze_point$gaze_y * (dims[3] - 1)
  }

  # Step 2: Voxel to world coordinates using NVImage
  voxel_coords <- c(voxel_x, voxel_y, voxel_z)
  frac_coords <- voxel_coords / nvimg$dimsRAS[2:4]
  world_coords <- nvimg$convertFrac2MM(frac_coords)[1:3]

  # Step 3: World back to voxel
  voxel_back <- nvimg$mm2vox(world_coords, frac = TRUE)

  # Step 4: Voxel back to screen
  if (plane == "AXIAL") {
    screen_x_back <- voxel_back[1] / (dims[1] - 1)
    screen_y_back <- voxel_back[2] / (dims[2] - 1)
    voxel_z_back <- voxel_back[3]
  } else if (plane == "SAGITTAL") {
    voxel_x_back <- voxel_back[1]
    screen_x_back <- voxel_back[2] / (dims[2] - 1)
    screen_y_back <- voxel_back[3] / (dims[3] - 1)
  } else {  # CORONAL
    screen_x_back <- voxel_back[1] / (dims[1] - 1)
    voxel_y_back <- voxel_back[2]
    screen_y_back <- voxel_back[3] / (dims[3] - 1)
  }

  # Calculate errors
  result <- data.frame(
    # Original values
    orig_screen_x = gaze_point$gaze_x,
    orig_screen_y = gaze_point$gaze_y,
    # Intermediate values
    voxel_x = voxel_x,
    voxel_y = voxel_y,
    voxel_z = voxel_z,
    world_x = world_coords[1],
    world_y = world_coords[2],
    world_z = world_coords[3],
    # Round-trip values
    voxel_x_back = voxel_back[1],
    voxel_y_back = voxel_back[2],
    voxel_z_back = voxel_back[3],
    screen_x_back = screen_x_back,
    screen_y_back = screen_y_back,
    # Errors
    voxel_error = sqrt(sum((voxel_coords - voxel_back)^2)),
    screen_error_x = abs(gaze_point$gaze_x - screen_x_back),
    screen_error_y = abs(gaze_point$gaze_y - screen_y_back),
    screen_error_total = sqrt((gaze_point$gaze_x - screen_x_back)^2 +
                                (gaze_point$gaze_y - screen_y_back)^2)
  )

  return(result)
}

#' Validate known target coordinates
#' @param target_patterns Gaze patterns with known targets
#' @param volume_info Volume information including affine
#' @return Data frame with target validation results
validate_targets <- function(target_patterns, volume_info) {
  # Extract unique targets from the fixation patterns
  targets <- target_patterns %>%
    filter(pattern_type == "fixation", !is.na(target_desc)) %>%
    group_by(target_id, target_desc) %>%
    summarise(
      mean_gaze_x = mean(gaze_x),
      mean_gaze_y = mean(gaze_y),
      plane = first(plane),
      slice_idx = first(slice_idx),
      .groups = "drop"
    )

  # For each target, reconstruct world coordinates and compare
  validation_results <- data.frame()

  for (i in 1:nrow(targets)) {
    target <- targets[i, ]

    # Expected world coordinates based on target description
    # Parse from the original target generation logic
    center_offset <- case_when(
      grepl("center", target$target_desc) ~ c(0, 0, 0),
      grepl("anterior", target$target_desc) ~ c(20, 0, 0),
      grepl("posterior", target$target_desc) ~ c(-20, 0, 0),
      grepl("left", target$target_desc) ~ c(-20, 0, 0),
      grepl("right", target$target_desc) ~ c(20, 0, 0),
      grepl("superior", target$target_desc) ~ c(0, 20, 0),
      grepl("inferior", target$target_desc) ~ c(0, -20, 0),
      grepl("corner_anterior_left_superior", target$target_desc) ~ c(40, 40, 30),
      grepl("corner_posterior_right_inferior", target$target_desc) ~ c(-40, -40, -30),
      TRUE ~ c(0, 0, 0)
    )

    # This is a placeholder - in real validation, we'd compute from known center
    # For now, we'll validate the transformation consistency
    result <- data.frame(
      target_id = target$target_id,
      target_desc = target$target_desc,
      plane = target$plane,
      slice_idx = target$slice_idx,
      gaze_x = target$mean_gaze_x,
      gaze_y = target$mean_gaze_y,
      expected_offset_x = center_offset[1],
      expected_offset_y = center_offset[2],
      expected_offset_z = center_offset[3]
    )

    validation_results <- rbind(validation_results, result)
  }

  return(validation_results)
}

# Run validation for all test cases
all_results <- list()
result_id <- 1

for (volume_name in names(volumes)) {
  message(sprintf("\nValidating: %s", volume_name))

  volume_info <- volumes[[volume_name]]

  # Load NIfTI data using gazeneuro
  # Create temporary NIfTI file for testing
  temp_file <- tempfile(fileext = ".nii.gz")
  writeNifti(volume_info$nifti, temp_file)

  # Load with gazeneuro's preload function
  nifti_data <- preload_nifti_data(temp_file)

  # Get gaze patterns for this volume
  volume_patterns <- gaze_patterns %>%
    filter(volume_name == !!volume_name)

  # Test 1: Round-trip transformation accuracy
  message("  Testing round-trip transformations...")

  # Sample gaze points for testing (every 100th point to reduce computation)
  test_points <- volume_patterns %>%
    group_by(pattern_id, plane, slice_idx) %>%
    slice(seq(1, n(), 100)) %>%
    ungroup()

  round_trip_results <- data.frame()

  for (i in 1:nrow(test_points)) {
    if (i %% 100 == 0) cat(".")

    point <- test_points[i, ]
    result <- test_round_trip(point, nifti_data, point$plane, point$slice_idx)
    result$pattern_id <- point$pattern_id
    result$pattern_name <- point$pattern_name
    result$volume_name <- volume_name
    result$plane <- point$plane
    result$slice_idx <- point$slice_idx

    round_trip_results <- rbind(round_trip_results, result)
  }
  cat("\n")

  # Test 2: Target coordinate validation
  message("  Testing known target coordinates...")
  target_patterns <- volume_patterns %>%
    filter(pattern_name == "target_fixations")

  if (nrow(target_patterns) > 0) {
    target_validation <- validate_targets(target_patterns, volume_info)
    target_validation$volume_name <- volume_name
  } else {
    target_validation <- data.frame()
  }

  # Store results
  all_results[[result_id]] <- list(
    volume_name = volume_name,
    test_case = volume_info$test_case,
    round_trip = round_trip_results,
    target_validation = target_validation,
    transformation_stats = list(
      mean_voxel_error = mean(round_trip_results$voxel_error),
      max_voxel_error = max(round_trip_results$voxel_error),
      mean_screen_error = mean(round_trip_results$screen_error_total),
      max_screen_error = max(round_trip_results$screen_error_total),
      sub_voxel_accuracy = mean(round_trip_results$voxel_error < 0.5)
    )
  )

  result_id <- result_id + 1

  # Clean up temp file
  unlink(temp_file)
}

# Aggregate results
message("\nAggregating results...")

# Combine all round-trip results
all_round_trips <- map_df(all_results, ~ .x$round_trip)

# Create summary by test case
summary_by_case <- all_round_trips %>%
  group_by(volume_name) %>%
  summarise(
    n_tests = n(),
    mean_voxel_error = mean(voxel_error),
    max_voxel_error = max(voxel_error),
    p95_voxel_error = quantile(voxel_error, 0.95),
    mean_screen_error = mean(screen_error_total),
    max_screen_error = max(screen_error_total),
    sub_voxel_rate = mean(voxel_error < 0.5),
    sub_pixel_rate = mean(screen_error_total < 1/512),  # Assuming 512x512 display
    .groups = "drop"
  )

# Test oblique vs non-oblique performance
volume_metadata <- map_df(volumes, ~ data.frame(
  volume_name = .x$test_case$name,
  is_oblique = .x$test_case$is_oblique
))

oblique_comparison <- all_round_trips %>%
  left_join(volume_metadata, by = "volume_name") %>%
  group_by(is_oblique) %>%
  summarise(
    n_volumes = n_distinct(volume_name),
    n_tests = n(),
    mean_voxel_error = mean(voxel_error),
    mean_screen_error = mean(screen_error_total),
    sub_voxel_rate = mean(voxel_error < 0.5),
    .groups = "drop"
  )

# Save results
qsave(all_results, "coordinate_validation/transformation_results.qs")
write_csv(summary_by_case, "coordinate_validation/transformation_summary.csv")
write_csv(oblique_comparison, "coordinate_validation/oblique_comparison.csv")

# Save detailed round-trip results
write_csv(all_round_trips %>% sample_n(min(10000, nrow(.))),
          "coordinate_validation/round_trip_sample.csv")

# Print summary
cat("\n=== COORDINATE TRANSFORMATION VALIDATION SUMMARY ===\n\n")

cat("Overall Performance:\n")
cat(sprintf("  Total transformations tested: %d\n", nrow(all_round_trips)))
cat(sprintf("  Mean voxel error: %.6f\n", mean(all_round_trips$voxel_error)))
cat(sprintf("  Sub-voxel accuracy: %.1f%%\n",
            100 * mean(all_round_trips$voxel_error < 0.5)))
cat(sprintf("  Mean screen error: %.6f (normalized units)\n",
            mean(all_round_trips$screen_error_total)))

cat("\nOblique vs Non-oblique:\n")
print(oblique_comparison)

cat("\nWorst performing cases:\n")
worst_cases <- summary_by_case %>%
  arrange(desc(mean_voxel_error)) %>%
  head(5)
print(worst_cases)

cat("\nResults saved to coordinate_validation/\n")
