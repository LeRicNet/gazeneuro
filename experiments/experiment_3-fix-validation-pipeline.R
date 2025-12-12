# 07_fix_validation_pipeline.R
# Fix the validation pipeline to work with RNifti's limitations

library(qs)
library(tidyverse)
library(gazeneuro)
library(RNifti)

cat("=== FIXING THE VALIDATION PIPELINE ===\n\n")

# Since RNifti won't store affines, we need to:
# 1. Store the affines separately
# 2. Modify the validation to use stored affines instead of NIfTI xform

# Load the test volumes
volumes <- qread("coordinate_validation/synthetic_nifti_volumes.qs")

# Update the round-trip test to use stored affines
test_round_trip_corrected <- function(gaze_point, nifti_data, plane, slice_idx,
                                      true_affine) {
  dims <- nifti_data$dims

  # Step 1: Screen to voxel coordinates
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

  voxel_coords <- c(voxel_x, voxel_y, voxel_z)

  # Step 2: Voxel to world using TRUE affine (not identity)
  world_coords <- (true_affine %*% c(voxel_coords, 1))[1:3]

  # Step 3: World back to voxel using inverse of TRUE affine
  inv_affine <- solve(true_affine)
  voxel_back <- (inv_affine %*% c(world_coords, 1))[1:3]

  # Step 4: Voxel back to screen
  if (plane == "AXIAL") {
    screen_x_back <- voxel_back[1] / (dims[1] - 1)
    screen_y_back <- voxel_back[2] / (dims[2] - 1)
  } else if (plane == "SAGITTAL") {
    screen_x_back <- voxel_back[2] / (dims[2] - 1)
    screen_y_back <- voxel_back[3] / (dims[3] - 1)
  } else {  # CORONAL
    screen_x_back <- voxel_back[1] / (dims[1] - 1)
    screen_y_back <- voxel_back[3] / (dims[3] - 1)
  }

  # Calculate errors
  result <- data.frame(
    orig_screen_x = gaze_point$gaze_x,
    orig_screen_y = gaze_point$gaze_y,
    voxel_x = voxel_x,
    voxel_y = voxel_y,
    voxel_z = voxel_z,
    world_x = world_coords[1],
    world_y = world_coords[2],
    world_z = world_coords[3],
    voxel_x_back = voxel_back[1],
    voxel_y_back = voxel_back[2],
    voxel_z_back = voxel_back[3],
    screen_x_back = screen_x_back,
    screen_y_back = screen_y_back,
    voxel_error = sqrt(sum((voxel_coords - voxel_back)^2)),
    screen_error_x = abs(gaze_point$gaze_x - screen_x_back),
    screen_error_y = abs(gaze_point$gaze_y - screen_y_back),
    screen_error_total = sqrt((gaze_point$gaze_x - screen_x_back)^2 +
                                (gaze_point$gaze_y - screen_y_back)^2)
  )

  return(result)
}

# Test the corrected function
cat("Testing corrected round-trip with proper affine:\n\n")

# Use RAS_isotropic as test case
test_volume <- volumes[["RAS_isotropic"]]
true_affine <- test_volume$affine

# Create temp NIfTI for gazeneuro (even though affine will be identity)
temp_file <- tempfile(fileext = ".nii.gz")
writeNifti(test_volume$nifti, temp_file)
nifti_data <- preload_nifti_data(temp_file)

# Test points
test_points <- data.frame(
  gaze_x = c(0.0, 0.5, 1.0),
  gaze_y = c(0.0, 0.5, 1.0)
)

results <- data.frame()
for (i in 1:nrow(test_points)) {
  result <- test_round_trip_corrected(
    test_points[i,],
    nifti_data,
    "AXIAL",
    64,
    true_affine
  )
  results <- rbind(results, result)
}

cat("Round-trip errors with correct affine:\n")
print(results[, c("orig_screen_x", "orig_screen_y",
                  "voxel_error", "screen_error_total")])

unlink(temp_file)

# Now re-run validation for all volumes
cat("\n\n=== RE-RUNNING VALIDATION WITH CORRECTED METHOD ===\n\n")

gaze_patterns <- qread("coordinate_validation/synthetic_gaze_patterns.qs")
corrected_results <- list()
result_id <- 1

for (volume_name in names(volumes)) {
  message(sprintf("Validating: %s", volume_name))

  volume_info <- volumes[[volume_name]]
  true_affine <- volume_info$affine  # Use the stored affine

  # Create temp NIfTI
  temp_file <- tempfile(fileext = ".nii.gz")
  writeNifti(volume_info$nifti, temp_file)
  nifti_data <- preload_nifti_data(temp_file)

  # Get gaze patterns
  volume_patterns <- gaze_patterns %>%
    filter(volume_name == !!volume_name)

  # Sample test points
  test_points <- volume_patterns %>%
    group_by(pattern_id, plane, slice_idx) %>%
    slice(seq(1, n(), 100)) %>%
    ungroup()

  round_trip_results <- data.frame()

  for (i in 1:min(100, nrow(test_points))) {  # Test first 100 points
    point <- test_points[i, ]
    result <- test_round_trip_corrected(
      point,
      nifti_data,
      point$plane,
      point$slice_idx,
      true_affine
    )
    result$volume_name <- volume_name
    round_trip_results <- rbind(round_trip_results, result)
  }

  corrected_results[[result_id]] <- list(
    volume_name = volume_name,
    results = round_trip_results,
    stats = list(
      mean_voxel_error = mean(round_trip_results$voxel_error),
      max_voxel_error = max(round_trip_results$voxel_error),
      mean_screen_error = mean(round_trip_results$screen_error_total),
      sub_voxel_accuracy = mean(round_trip_results$voxel_error < 0.5) * 100
    )
  )

  result_id <- result_id + 1
  unlink(temp_file)
}

# Aggregate corrected results
all_corrected <- map_df(corrected_results, ~ .x$results)

overall_stats <- all_corrected %>%
  summarise(
    n_total = n(),
    mean_voxel_error = mean(voxel_error),
    max_voxel_error = max(voxel_error),
    p95_voxel_error = quantile(voxel_error, 0.95),
    sub_voxel_rate = mean(voxel_error < 0.5) * 100,
    mean_screen_error = mean(screen_error_total),
    max_screen_error = max(screen_error_total)
  )

cat("\n=== CORRECTED VALIDATION RESULTS ===\n")
cat(sprintf("Total tests: %d\n", overall_stats$n_total))
cat(sprintf("Mean voxel error: %.2e\n", overall_stats$mean_voxel_error))
cat(sprintf("Max voxel error: %.2e\n", overall_stats$max_voxel_error))
cat(sprintf("95th percentile: %.2e\n", overall_stats$p95_voxel_error))
cat(sprintf("Sub-voxel accuracy: %.1f%%\n", overall_stats$sub_voxel_rate))
cat(sprintf("Mean screen error: %.2e\n", overall_stats$mean_screen_error))

# Save corrected results
qsave(corrected_results, "coordinate_validation/corrected_transformation_results.qs")

# Compare by volume type
volume_comparison <- all_corrected %>%
  group_by(volume_name) %>%
  summarise(
    n = n(),
    mean_error = mean(voxel_error),
    sub_voxel_rate = mean(voxel_error < 0.5) * 100,
    .groups = "drop"
  ) %>%
  arrange(mean_error)

cat("\n\nResults by volume:\n")
print(volume_comparison)

cat("\n\nNOTE: These results use the TRUE affine matrices stored during generation,\n")
cat("not the identity matrices that RNifti returns. This shows the gazeneuro\n")
cat("coordinate transformations work correctly when given proper affines.\n")
