# gazeneuro

<!-- badges: start -->
[![R-CMD-check](https://github.comlericnet/gazeneuro/workflows/R-CMD-check/badge.svg)](https://github.com/lericnet/gazeneuro/actions)
<!-- badges: end -->

The `gazeneuro` package provides tools for integrating eye tracking gaze data with NIfTI neuroimaging files. It enables temporal alignment of gaze tracking data with brain slice viewing events, coordinate transformations between voxel and world spaces, and comprehensive visualization of gaze patterns on neuroimaging data.

## Features

- **Temporal Integration**: Automatically align gaze tracking timestamps with neuroimaging slice viewing events
- **Image Bounds Support**: Handle cases where the neuroimaging slice doesn't fill the entire display area
- **Coordinate Transformation**: Convert between display, image, voxel, and world coordinates using gl-matrix style transformations
- **Visualization Tools**: 
  - Overlay gaze data on individual brain slices
  - Create grid visualizations of all viewed slices
  - Generate gaze density heatmaps
  - Export animation frames for video creation
- **Analysis Functions**: Calculate gaze statistics by slice, viewing duration, and fixation patterns

## Installation

You can install the development version of gazeneuro from GitHub:

```r
# install.packages("devtools")
devtools::install_github("lericnet/gazeneuro")
```

## Quick Start

```r
library(gazeneuro)

# Load NIfTI data
nifti_data <- preload_nifti_data("path/to/your/brain.nii.gz")

# Define image bounds if the image doesn't fill the display
# (e.g., image occupies center 80% of display)
image_bounds <- list(
  left = 0.1,    # 10% margin on left
  right = 0.9,   # 10% margin on right
  top = 0.1,     # 10% margin on top
  bottom = 0.9   # 10% margin on bottom
)

# Integrate gaze tracking data with slice viewing events
# Image bounds are applied during integration
integrated <- integrate_all_gaze_points(
  gaze_data, 
  z_axis_data,
  image_bounds = image_bounds
)

# Visualize gaze on a specific slice (uses adjusted coordinates automatically)
plot_slice_with_all_gaze(nifti_data, integrated, slice_num = 12)

# Create a grid view of all viewed slices
summary_stats <- plot_all_slices_with_gaze(nifti_data, integrated)

# Map gaze to anatomical coordinates
anatomical <- map_gaze_to_anatomy(nifti_data, integrated)
```

## Image Bounds and Coordinate Transformations

### When to Use Image Bounds

Image bounds are necessary when:
- The neuroimaging slice has margins or padding in the display
- The image is letterboxed/pillarboxed due to aspect ratio differences
- The display includes UI elements (headers, footers, sidebars)
- The image only occupies a portion of the screen

### Coordinate Transformation Pipeline

The package handles multiple coordinate systems:

1. **Display Coordinates** (0-1): Raw gaze coordinates from eye tracker
2. **Image Coordinates** (0-1): Adjusted for image position within display
3. **Voxel Coordinates**: Integer indices into the NIfTI data array
4. **World Coordinates** (mm): Physical positions in scanner space

```r
# The transformation pipeline:
# Display coords → Image coords → Voxel coords → World coords

# This happens automatically when you:
integrated <- integrate_all_gaze_points(gaze_data, z_axis_data, 
                                       image_bounds = image_bounds)

# View the coordinate transformations:
anatomical <- map_gaze_to_anatomy(integrated, nifti_data)
anatomical %>%
  select(gaze_x, gaze_y,              # Original display coords
         gaze_x_adjusted, gaze_y_adjusted,  # Adjusted for image bounds
         voxel_x, voxel_y, voxel_z,    # Voxel indices
         world_x, world_y, world_z)     # World coordinates (mm)
```

### Automatic Bounds Detection

If you don't know the exact image bounds, you can estimate them from the gaze data:

```r
# Function to detect bounds from gaze density
detect_image_bounds <- function(gaze_data, threshold_pct = 2) {
  x_bounds <- quantile(gaze_data$gaze_x, c(threshold_pct/100, 1-threshold_pct/100))
  y_bounds <- quantile(gaze_data$gaze_y, c(threshold_pct/100, 1-threshold_pct/100))
  
  list(
    left = x_bounds[1],
    right = x_bounds[2], 
    top = y_bounds[1],
    bottom = y_bounds[2]
  )
}

auto_bounds <- detect_image_bounds(gaze_data)
integrated <- integrate_all_gaze_points(gaze_data, z_axis_data, 
                                       image_bounds = auto_bounds)
```

## Data Format

### Gaze Data
The gaze tracking data should be a data frame with the following columns:
- `device_time_stamp`: Timestamp in microseconds
- `gaze_point_on_display_area_x`: X coordinate (0-1)
- `gaze_point_on_display_area_y`: Y coordinate (0-1)

### Z-axis Data
The slice viewing event data should contain:
- `client_timestamp`: Timestamp in milliseconds
- `plane`: Imaging plane (AXIAL, SAGITTAL, or CORONAL)
- `index`: Normalized slice index (0-1)
- `image_id`: Image identifier

## Main Functions

### Data Loading and Integration
- `preload_nifti_data()`: Load NIfTI file with pre-extracted data
- `integrate_all_gaze_points()`: Temporally align gaze with slice viewing, apply image bounds

### Visualization
- `plot_slice_with_all_gaze()`: Plot single slice with gaze overlay
- `plot_all_slices_with_gaze()`: Grid view of all viewed slices
- `create_gaze_animation()`: Generate animation frames
- `create_gaze_summary_plot()`: Create summary visualizations

### Analysis
- `analyze_gaze_by_slice()`: Calculate gaze statistics per slice
- `map_gaze_to_anatomy()`: Convert gaze to anatomical coordinates
- `check_integration_quality()`: Diagnostic plots for integration quality

### Interactive Tools
- `safe_coordinate_picker()`: Click on slices to get coordinates
- `safe_quick_view()`: Quick slice viewer with coordinate info

## Example Workflow

```r
library(gazeneuro)
library(tidyverse)

# 1. Load your data
nifti_data <- preload_nifti_data("brain_scan.nii.gz")

# 2. Define image bounds (if needed)
# Example: letterboxed square image in 16:9 display
display_aspect <- 16/9
image_aspect <- 1
image_width <- image_aspect / display_aspect

image_bounds <- list(
  left = (1 - image_width) / 2,
  right = 1 - (1 - image_width) / 2,
  top = 0,
  bottom = 1
)

# 3. Integrate with image bounds
integrated <- integrate_all_gaze_points(
  gaze_data = your_gaze_data,
  z_axis = your_slice_events,
  image_bounds = image_bounds
)

# 4. Check integration quality
check_integration_quality(your_gaze_data, your_slice_events, integrated)

# 5. Analyze gaze patterns
# Only analyzes points within image bounds
gaze_stats <- analyze_gaze_by_slice(integrated, nifti_data)

# Check how many points were outside bounds
integrated %>%
  summarise(
    total = n(),
    within_bounds = sum(within_bounds),
    outside_bounds = sum(!within_bounds),
    pct_within = 100 * mean(within_bounds)
  )

# 6. Create visualizations
# All functions automatically use adjusted coordinates
plot_slice_with_all_gaze(nifti_data, integrated, slice_num = 15)
plot_all_slices_with_gaze(nifti_data, integrated, plane = "AXIAL")
create_gaze_summary_plot(nifti_data, integrated)

# 7. Map to anatomical space
anatomical <- map_gaze_to_anatomy(integrated, nifti_data)

# 8. Export results
write.csv(gaze_stats, "gaze_statistics.csv")
write.csv(anatomical, "gaze_anatomical_coordinates.csv")
```

## Advanced Features

### Coordinate Transformations with NVImage

The package includes a complete implementation of gl-matrix operations for neuroimaging coordinate transformations:

```r
# Access the NVImage object
nvimg <- nifti_data$nvimage

# Convert mm to voxel coordinates
voxel_coords <- nvimg$mm2vox(c(10, -20, 30))

# Convert fractional to mm coordinates
mm_coords <- nvimg$convertFrac2MM(c(0.5, 0.5, 0.5))
```

### Anatomical Region Lookup

If you have an atlas file, you can identify anatomical regions:

```r
# Load atlas
atlas_data <- preload_nifti_data("atlas.nii.gz")

# Get anatomical regions for gaze points
regions <- get_anatomical_regions(anatomical_coords, atlas_data)
```

## Tips for Image Bounds

1. **Measure actual bounds**: Use a calibration image or known markers to determine exact bounds
2. **Consider UI elements**: Account for any headers, footers, or sidebars in your display
3. **Check aspect ratios**: Calculate letterbox/pillarbox bounds for mismatched aspect ratios
4. **Validate with visualization**: Use `check_integration_quality()` to verify bounds are correct
5. **Store bounds metadata**: Save image bounds with your experimental data for reproducibility

## Citation

If you use this package in your research, please cite:

```
@software{gazeneuro,
  title = {gazeneuro: Gaze Tracking and Neuroimaging Integration},
  author = {Eric Prince},
  year = {2025},
  url = {https://github.com/lericnet/gazeneuro}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Issues

If you encounter any problems, please [file an issue](https://github.com/lericnet/gazeneuro/issues) along with a detailed description.
