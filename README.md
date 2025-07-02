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

## Coordinate Mapping Functions

### Anatomical Coordinate Conversion
- `map_gaze_to_anatomy()`: Convert all gaze points to anatomical coordinates
- `locate_single_point()`: Map a single gaze point (equivalent to JavaScript locate())
- `batch_locate_gaze()`: Process and summarize gaze locations
- `get_slice_number()`: Convert normalized slice indices to slice numbers

### Utilities
- `resolve_coordinates()`: Adjust coordinates for device pixel ratio and display frame
- `safe_get_value()`: Safely extract intensity values from image data
- `export_anatomical_locations()`: Export locations as CSV, JSON, or NIfTI heatmap

### Display Frame Handling
- `visualize_display_frame()`: Visualize the display layout
- `plot_slice_with_validity()`: Show valid vs out-of-bounds gaze points

## Coordinate Mapping Example

The package now includes functions to map gaze coordinates to anatomical locations in neuroimaging space:

```r
# Map gaze points to anatomical coordinates
anatomical_locations <- map_gaze_to_anatomy(integrated, nifti_data)

# Find most-viewed brain regions
most_viewed <- batch_locate_gaze(integrated, nifti_data, summarize = TRUE)

# Map a single point (like JavaScript locate())
location <- locate_single_point(
  x = 0.5,  # gaze x (0-1)
  y = 0.5,  # gaze y (0-1)
  z = 0.5,  # slice position (0-1)
  nifti_data = nifti_data
)

# Access coordinates
print(location$mm)   # World coordinates in mm
print(location$vox)  # Voxel indices
print(location$values[[1]]$value)  # Intensity value

# Export as NIfTI heatmap
export_anatomical_locations(
  anatomical_locations, 
  "gaze_heatmap.nii.gz",
  format = "nifti",
  nifti_template = nifti_data
)
```

### Display Frame Configuration

The package handles a specific display setup where the image canvas (1924×1560) is positioned within a larger frame (2624×1640):

```r
# Visualize the display layout
visualize_display_frame()

# Check which gazes fall within the image canvas
plot_slice_with_validity(nifti_data, integrated, slice_num = 12)
```

**Important**: When plotting gaze data, the package handles Y-axis inversion between Tobii coordinates (Y=0 at top) and plot coordinates (Y=0 at bottom) automatically.

This allows you to:
- Track which brain regions participants focused on
- Filter out gazes that fall on screen borders
- Create heatmaps of gaze density in anatomical space
- Export data for analysis in other neuroimaging tools
- Integrate eye tracking with brain imaging analysis pipelines

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
