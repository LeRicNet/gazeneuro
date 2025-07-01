library(gazeneuro)
library(tidyverse)

# 1. Load NIfTI data
nifti_data <- preload_nifti_data("../2025_ATPC_Study/data/img/1.3.12.2.1107.5.2.30.26451.2007113015563474695781235.0.0.0.nii.gz")

# 2. Load gaze and z-axis data
gaze_data <- read_csv('../2025_ATPC_Study/data/gaze_data/records-drzmexhret-dab29f3dcfc6a92938d61db1d396e492c4590042625149988a148b9a95ccddf1.csv')
z_axis <- read_csv('../2025_ATPC_Study/data/scroll_index_atpc.csv')

# Filter z-axis to match gaze data session
z_axis <- z_axis %>%
  filter(narrative_session == gaze_data$tracking_session[1])

# 3. Integrate the data
integrated <- integrate_all_gaze_points(gaze_data, z_axis)

# 4. Original visualization
plot_slice_with_all_gaze(nifti_data, integrated, slice_num = 13)

# 5. NEW: Map gaze points to anatomical coordinates
cat("\n=== Coordinate Mapping ===\n")
anatomical_locations <- map_gaze_to_anatomy(integrated, nifti_data)

# Show first few anatomical locations
cat("\nFirst 5 gaze points in anatomical space:\n")
anatomical_locations %>%
  select(time_sec, gaze_x, gaze_y, mm_x, mm_y, mm_z, vox_x, vox_y, vox_z, intensity) %>%
  head(5) %>%
  print()

# 6. Find most-viewed voxels
most_viewed <- batch_locate_gaze(integrated, nifti_data, summarize = TRUE)
cat("\nTop 5 most-viewed voxels:\n")
print(head(most_viewed, 5))

# 7. Example: Single point location (like JavaScript locate())
# Get a specific gaze point
example_point <- integrated[100, ]  # 100th gaze point
location <- locate_single_point(
  x = example_point$gaze_x,
  y = example_point$gaze_y,
  z = example_point$slice_index,
  nifti_data = nifti_data,
  plane = example_point$plane
)

cat("\nExample single point location:\n")
cat(sprintf("  Gaze: [%.3f, %.3f]\n", location$xy[1], location$xy[2]))
cat(sprintf("  MM: [%.1f, %.1f, %.1f]\n", location$mm[1], location$mm[2], location$mm[3]))
cat(sprintf("  Voxel: [%d, %d, %d]\n", location$vox[1], location$vox[2], location$vox[3]))
cat(sprintf("  Intensity: %.1f\n", location$values[[1]]$value))

# 8. Visualize coordinate mapping
plot_coordinate_mapping(anatomical_locations, nifti_data)

# 9. Export results
# CSV format (for statistical analysis)
export_anatomical_locations(anatomical_locations, "gaze_anatomical_coords.csv")

# NIfTI heatmap (for visualization in neuroimaging software)
export_anatomical_locations(anatomical_locations, "gaze_heatmap.nii.gz",
                            format = "nifti", nifti_template = nifti_data)

cat("\nCoordinate mapping complete. Results exported.\n")
