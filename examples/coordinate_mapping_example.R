# Example: Using Coordinate Mapping in gazeneuro
# This example demonstrates how to map gaze coordinates to anatomical locations

library(gazeneuro)
library(tidyverse)

# 1. Load your NIfTI data
nifti_data <- preload_nifti_data("path/to/your/brain.nii.gz")

# 2. Load and integrate gaze tracking data
gaze_data <- read_csv("path/to/gaze_data.csv")
z_axis_data <- read_csv("path/to/z_axis_data.csv")

integrated <- integrate_all_gaze_points(gaze_data, z_axis_data)

# 3. Map gaze points to anatomical coordinates
cat("Mapping gaze points to anatomical space...\n")
anatomical_locations <- map_gaze_to_anatomy(integrated, nifti_data)

# 4. Example: Find which voxel was looked at most
most_viewed <- anatomical_locations %>%
  group_by(vox_x, vox_y, vox_z) %>%
  summarise(
    n_views = n(),
    total_time = sum(time_aligned),
    mean_intensity = mean(intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_views)) %>%
  slice(1)

cat("\nMost viewed voxel:\n")
cat(sprintf("  Voxel: [%d, %d, %d]\n",
            most_viewed$vox_x, most_viewed$vox_y, most_viewed$vox_z))
cat(sprintf("  Number of gazes: %d\n", most_viewed$n_views))
cat(sprintf("  Total viewing time: %.2f seconds\n", most_viewed$total_time))
cat(sprintf("  Mean intensity: %.2f\n", most_viewed$mean_intensity))

# 5. Example: Create a simple report of gaze patterns by slice
slice_report <- anatomical_locations %>%
  group_by(plane, slice_num) %>%
  summarise(
    n_gazes = n(),
    unique_voxels = n_distinct(paste(vox_x, vox_y)),
    mean_x_mm = mean(mm_x),
    mean_y_mm = mean(mm_y),
    spread_mm = sqrt(var(mm_x) + var(mm_y)),
    .groups = "drop"
  ) %>%
  arrange(plane, slice_num)

cat("\nGaze patterns by slice (top 10):\n")
print(head(slice_report, 10))

# 6. Example: Visualize gaze locations on a specific slice
# Choose a slice with many gaze points
target_slice <- slice_report %>%
  filter(plane == "AXIAL") %>%
  arrange(desc(n_gazes)) %>%
  slice(1) %>%
  pull(slice_num)

# Get gaze points for that slice
slice_gazes <- anatomical_locations %>%
  filter(plane == "AXIAL", slice_num == target_slice)

# Create visualization with the actual slice image
slice_data <- nifti_data$data[,, target_slice]

par(mfrow = c(1, 2))

# Left: Gaze overlay on brain slice
image(slice_data, col = gray((0:255)/255),
      main = sprintf("Gaze on Brain - Slice %d", target_slice),
      xlab = "X", ylab = "Y")
points(slice_gazes$gaze_x, slice_gazes$gaze_y,
       pch = 19, col = rgb(1, 0, 0, 0.3))

# Add mean gaze position
mean_x <- mean(slice_gazes$gaze_x)
mean_y <- mean(slice_gazes$gaze_y)
points(mean_x, mean_y, col = "blue", pch = 3, cex = 2, lwd = 2)

# Right: Anatomical coordinates
plot(slice_gazes$mm_x, slice_gazes$mm_y,
     pch = 19, col = rgb(0, 0, 1, 0.3),
     main = "Anatomical Coordinates",
     xlab = "X (mm)", ylab = "Y (mm)")

# Add coordinate grid
grid()

par(mfrow = c(1, 1))

# 7. Example: Time-based analysis
# How does gaze focus change over time?
time_windows <- anatomical_locations %>%
  mutate(
    time_window = cut(time_sec, breaks = 20)
  ) %>%
  group_by(time_window) %>%
  summarise(
    n_gazes = n(),
    spatial_variance = var(mm_x) + var(mm_y) + var(mm_z),
    mean_intensity = mean(intensity, na.rm = TRUE),
    .groups = "drop"
  )

# Plot spatial variance over time
plot(1:nrow(time_windows), time_windows$spatial_variance,
     type = "b", pch = 19,
     xlab = "Time Window", ylab = "Spatial Variance (mm²)",
     main = "Gaze Focus Over Time")

# 8. Save results for further analysis
cat("\nSaving anatomical location data...\n")
write_csv(anatomical_locations, "gaze_anatomical_locations.csv")

# Create summary report
summary_report <- list(
  total_gazes = nrow(anatomical_locations),
  unique_voxels = anatomical_locations %>%
    distinct(vox_x, vox_y, vox_z) %>%
    nrow(),
  spatial_extent_mm = list(
    x_range = range(anatomical_locations$mm_x),
    y_range = range(anatomical_locations$mm_y),
    z_range = range(anatomical_locations$mm_z)
  ),
  planes_viewed = unique(anatomical_locations$plane),
  total_duration = max(anatomical_locations$time_sec) -
    min(anatomical_locations$time_sec)
)

cat("\nSummary Report:\n")
cat(sprintf("  Total gaze points mapped: %d\n", summary_report$total_gazes))
cat(sprintf("  Unique voxels viewed: %d\n", summary_report$unique_voxels))
cat(sprintf("  Total duration: %.2f seconds\n", summary_report$total_duration))
cat(sprintf("  Spatial extent: X[%.1f, %.1f], Y[%.1f, %.1f], Z[%.1f, %.1f] mm\n",
            summary_report$spatial_extent_mm$x_range[1],
            summary_report$spatial_extent_mm$x_range[2],
            summary_report$spatial_extent_mm$y_range[1],
            summary_report$spatial_extent_mm$y_range[2],
            summary_report$spatial_extent_mm$z_range[1],
            summary_report$spatial_extent_mm$z_range[2]))

# 9. Optional: Export for visualization in other tools
# Create a simple point cloud file
point_cloud <- anatomical_locations %>%
  select(mm_x, mm_y, mm_z, intensity) %>%
  na.omit()

write.table(point_cloud, "gaze_point_cloud.txt",
            row.names = FALSE, col.names = FALSE)
cat("\nPoint cloud exported to gaze_point_cloud.txt\n")

cat("\nCoordinate mapping complete!\n")
