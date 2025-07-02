# Example: Understanding and Using Display Frame Coordinate Adjustment
# This demonstrates how Tobii gaze coordinates are transformed to image canvas coordinates

library(gazeneuro)
library(tidyverse)

# 1. Visualize the display layout
cat("=== Display Frame Layout ===\n")
cat("Full frame: 2624 x 1640 pixels\n")
cat("Image canvas: 1924 x 1560 pixels\n")
cat("Left/Right borders: 350 pixels each\n")
cat("Top border: 80 pixels\n")
cat("Bottom border: 0 pixels (canvas aligned to bottom)\n")
cat("Device Pixel Ratio: 1.25\n\n")

# Show the visual layout
visualize_display_frame(show_example_points = TRUE)

# 2. Test coordinate transformations
cat("\n=== Coordinate Transformation Examples ===\n")

test_points <- data.frame(
  description = c(
    "Center of display",
    "Top-left of display",
    "Top-left of canvas (approximate)",
    "Center of canvas (approximate)",
    "Bottom-right of canvas (approximate)"
  ),
  tobii_x = c(0.5, 0.0, 0.133, 0.5, 0.867),
  tobii_y = c(0.5, 0.0, 0.049, 0.5, 1.0)
)

for (i in 1:nrow(test_points)) {
  adj <- resolve_coordinates(test_points$tobii_x[i], test_points$tobii_y[i])

  cat(sprintf("\n%s:\n", test_points$description[i]))
  cat(sprintf("  Tobii coords: [%.3f, %.3f]\n",
              test_points$tobii_x[i], test_points$tobii_y[i]))
  cat(sprintf("  Canvas coords: [%.3f, %.3f]\n", adj$x, adj$y))
  cat(sprintf("  In bounds: %s\n", ifelse(adj$in_bounds, "YES", "NO")))
}

# 3. Process real gaze data with coordinate adjustment
cat("\n\n=== Processing Gaze Data with Coordinate Adjustment ===\n")

# Load your data (example paths)
# nifti_data <- preload_nifti_data("brain.nii.gz")
# gaze_data <- read_csv("gaze_data.csv")
# z_axis_data <- read_csv("z_axis_data.csv")

# For demonstration, create sample data
set.seed(42)
sample_integrated <- data.frame(
  gaze_id = 1:100,
  time_sec = seq(0, 10, length.out = 100),
  time_aligned = seq(0, 10, length.out = 100),
  # Simulate gaze mostly in center with some outliers
  gaze_x = rnorm(100, mean = 0.5, sd = 0.15),
  gaze_y = rnorm(100, mean = 0.5, sd = 0.15),
  slice_index = rep(0.5, 100),
  plane = rep("AXIAL", 100),
  image_id = rep("test", 100)
)

# Count how many points fall outside the canvas
out_of_bounds <- sapply(1:nrow(sample_integrated), function(i) {
  adj <- resolve_coordinates(sample_integrated$gaze_x[i], sample_integrated$gaze_y[i])
  !adj$in_bounds
})

cat(sprintf("Total gaze points: %d\n", nrow(sample_integrated)))
cat(sprintf("Points outside canvas: %d (%.1f%%)\n",
            sum(out_of_bounds), 100 * sum(out_of_bounds) / nrow(sample_integrated)))

# 4. Demonstrate the difference in mapping with and without adjustment
cat("\n=== Coordinate Mapping Comparison ===\n")

# Example point near the edge
edge_point <- data.frame(
  gaze_id = 1,
  time_sec = 0,
  time_aligned = 0,
  gaze_x = 0.9,  # Near right edge of display
  gaze_y = 0.5,
  slice_index = 0.5,
  plane = "AXIAL",
  image_id = "test"
)

# For demonstration with dummy NIfTI data
dummy_nifti <- list(
  nvimage = NVImage$new(
    matRAS = diag(4),
    dimsRAS = c(3, 256, 256, 128),
    pixDimsRAS = c(1, 1, 1, 1)
  ),
  data = array(0, dim = c(256, 256, 128)),
  dims = c(256, 256, 128),
  nifti = list(.xData = list(.sourcePath = "dummy.nii"))
)
dummy_nifti$nvimage$calculateOblique()

# Map without adjustment
loc_no_adj <- map_gaze_to_anatomy(edge_point, dummy_nifti, adjust_coordinates = FALSE)

# Map with adjustment
loc_with_adj <- map_gaze_to_anatomy(edge_point, dummy_nifti, adjust_coordinates = TRUE)

cat("\nEdge point (x=0.9, y=0.5):\n")
cat("Without adjustment:\n")
cat(sprintf("  Canvas coords: [%.3f, %.3f]\n", loc_no_adj$gaze_x, loc_no_adj$gaze_y))
cat(sprintf("  Voxel: [%d, %d, %d]\n", loc_no_adj$vox_x, loc_no_adj$vox_y, loc_no_adj$vox_z))

cat("With adjustment:\n")
cat(sprintf("  Canvas coords: [%.3f, %.3f]\n", loc_with_adj$gaze_x, loc_with_adj$gaze_y))
cat(sprintf("  In canvas: %s\n", ifelse(loc_with_adj$in_canvas, "YES", "NO")))
cat(sprintf("  Voxel: [%d, %d, %d]\n", loc_with_adj$vox_x, loc_with_adj$vox_y, loc_with_adj$vox_z))

# 5. Best practices
cat("\n\n=== Best Practices ===\n")
cat("1. Always use adjust_coordinates=TRUE when processing Tobii data\n")
cat("2. Check the 'in_canvas' field to filter out off-canvas gazes\n")
cat("3. The 'gaze_x_original' and 'gaze_y_original' fields preserve Tobii coordinates\n")
cat("4. Use visualize_display_frame() to verify your setup\n")

# 6. Filtering example
cat("\n=== Filtering Out-of-Bounds Gazes ===\n")

# Process all points
all_locations <- map_gaze_to_anatomy(sample_integrated, dummy_nifti)

# Filter to only in-canvas gazes
valid_locations <- all_locations %>%
  filter(in_canvas)

cat(sprintf("Total locations: %d\n", nrow(all_locations)))
cat(sprintf("Valid (in-canvas) locations: %d\n", nrow(valid_locations)))
cat(sprintf("Filtered out: %d\n", nrow(all_locations) - nrow(valid_locations)))
