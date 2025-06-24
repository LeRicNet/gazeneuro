# gazeneuro

<!-- badges: start -->
[![R-CMD-check](https://github.comlericnet/gazeneuro/workflows/R-CMD-check/badge.svg)](https://github.com/lericnet/gazeneuro/actions)
<!-- badges: end -->

The `gazeneuro` package provides tools for integrating eye tracking gaze data with NIfTI neuroimaging files. It enables temporal alignment of gaze tracking data with brain slice viewing events, coordinate transformations between voxel and world spaces, and comprehensive visualization of gaze patterns on neuroimaging data.

## Features

- **Temporal Integration**: Automatically align gaze tracking timestamps with neuroimaging slice viewing events
- **Coordinate Transformation**: Convert between voxel, fractional, and world coordinates using gl-matrix style transformations
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

# Integrate gaze tracking data with slice viewing events
integrated <- integrate_all_gaze_points(gaze_data, z_axis_data)

# Visualize gaze on a specific slice
plot_slice_with_all_gaze(nifti_data, integrated, slice_num = 12)

# Create a grid view of all viewed slices
summary_stats <- plot_all_slices_with_gaze(nifti_data, integrated)

# Generate animation frames
create_gaze_animation(nifti_data, integrated, output_dir = "gaze_animation")
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
- `integrate_all_gaze_points()`: Temporally align gaze with slice viewing

### Visualization
- `plot_slice_with_all_gaze()`: Plot single slice with gaze overlay
- `plot_all_slices_with_gaze()`: Grid view of all viewed slices
- `create_gaze_animation()`: Generate animation frames
- `create_gaze_summary_plot()`: Create summary visualizations

### Analysis
- `analyze_gaze_by_slice()`: Calculate gaze statistics per slice
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

# 2. Integrate gaze tracking with slice viewing
integrated <- integrate_all_gaze_points(
  gaze_data = your_gaze_data,
  z_axis = your_slice_events
)

# 3. Check integration quality
check_integration_quality(your_gaze_data, your_slice_events, integrated)

# 4. Analyze gaze patterns
gaze_stats <- analyze_gaze_by_slice(integrated, nifti_data)

# 5. Create visualizations
# Single slice with gaze
plot_slice_with_all_gaze(nifti_data, integrated, slice_num = 15)

# All viewed slices
plot_all_slices_with_gaze(nifti_data, integrated, plane = "AXIAL")

# Summary plots
create_gaze_summary_plot(nifti_data, integrated)

# 6. Export for further analysis
write.csv(gaze_stats, "gaze_statistics.csv")
```

## Advanced Features

### Coordinate Transformations
The package includes a complete implementation of gl-matrix operations for neuroimaging coordinate transformations:

```r
# Create NVImage object for coordinate transformations
nvimg <- nifti_data$nvimage

# Convert mm to voxel coordinates
voxel_coords <- nvimg$mm2vox(c(10, -20, 30))

# Convert fractional to mm coordinates
mm_coords <- nvimg$convertFrac2MM(c(0.5, 0.5, 0.5))
```

### Interactive Coordinate Picker
```r
# Click on image to get voxel and world coordinates
clicked_points <- safe_coordinate_picker(nifti_data, slice = 20)
```

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
