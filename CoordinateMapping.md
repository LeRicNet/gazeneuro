# Coordinate System Documentation for gazeneuro

## Display Layout

The eye tracking setup has a specific display configuration:

```
Full Display Frame: 2624 x 1640 pixels
├── Left Border: 350 pixels
├── Image Canvas: 1924 x 1560 pixels  
├── Right Border: 350 pixels
├── Top Border: 80 pixels
└── Bottom Border: 0 pixels (canvas aligned to bottom)

Device Pixel Ratio (DPR): 1.25
```

The canvas is **aligned to the bottom** of the display frame.

## Coordinate Systems

### 1. Tobii Coordinates (Input)
- **Range**: 0-1 normalized to the full display frame (2624 x 1640)
- **Origin**: Top-left of the display (Y=0 at top, Y=1 at bottom)
- **What it represents**: Where the user is looking on the entire screen

### 2. Canvas Coordinates (Intermediate)
- **Range**: 0-1 normalized to the image canvas (1924 x 1560)
- **Origin**: Top-left of the canvas area (Y=0 at top, Y=1 at bottom)
- **What it represents**: Where the user is looking within the actual image area

### 3. Plot Coordinates (Display)
- **Range**: 0-1 in plot space
- **Origin**: Bottom-left (Y=0 at bottom, Y=1 at top) - standard R plot orientation
- **Conversion**: `plot_y = 1 - tobii_y` or `plot_y = 1 - canvas_y`

### 4. Anatomical Coordinates (Output)
- **Voxel coordinates**: Integer indices into the 3D image array
- **MM coordinates**: Physical location in millimeters
- **What it represents**: The anatomical location being viewed

## Important: Y-Axis Inversion

When plotting gaze data on brain slices, you must invert the Y coordinate:

```r
# Tobii/Canvas Y: 0=top, 1=bottom
# Plot Y: 0=bottom, 1=top

# Always invert Y for plotting:
plot_y <- 1 - tobii_y
```

This ensures that:
- Top-left gazes appear in the top-left of the plot
- Bottom-right gazes appear in the bottom-right of the plot
- The canvas boundaries are correctly oriented

## Coordinate Transformation

The `resolve_coordinates()` function handles the transformation:

```r
# Input: Tobii coordinates (0-1 on full display)
tobii_x = 0.5  # Center of display
tobii_y = 0.5

# Transform to canvas coordinates
adjusted <- resolve_coordinates(tobii_x, tobii_y)
# Returns: list(x = 0.5, y = 0.476, in_bounds = TRUE)
```

### Key Boundaries

| Location | Tobii X | Tobii Y | Canvas X | Canvas Y | In Bounds |
|----------|---------|---------|----------|----------|-----------|
| Display center | 0.500 | 0.500 | 0.500 | 0.474 | ✓ |
| Canvas top-left | 0.133 | 0.049 | 0.000 | 0.000 | ✓ |
| Canvas bottom-right | 0.867 | 1.000 | 1.000 | 1.000 | ✓ |
| Display top-left | 0.000 | 0.000 | -0.182 | -0.051 | ✗ |
| Display bottom-left | 0.000 | 1.000 | -0.182 | 1.026 | ✗ |

## Usage in the Package

### Default Behavior
All coordinate mapping functions use `adjust_coordinates = TRUE` by default:

```r
# Automatically adjusts Tobii to canvas coordinates
locations <- map_gaze_to_anatomy(integrated_data, nifti_data)

# Access both coordinate systems
locations$gaze_x_original  # Original Tobii coordinates
locations$gaze_x          # Adjusted canvas coordinates
locations$in_canvas       # Whether point is within image area
```

### Filtering Out-of-Bounds Gazes

```r
# Only analyze gazes that fall within the image canvas
valid_locations <- locations %>%
  filter(in_canvas)
```

### Manual Coordinate Adjustment

```r
# For custom processing
adj <- resolve_coordinates(
  x = 0.7,              # Tobii X
  y = 0.3,              # Tobii Y
  frame_width = 2624,   # Can override defaults
  frame_height = 1640,
  canvas_width = 1924,
  canvas_height = 1560,
  left_border = 350,
  top_border = 80,
  dpr = 1.25
)
```

## Visualization

Use `visualize_display_frame()` to see the layout:

```r
# Shows frame, canvas, and example coordinate mappings
visualize_display_frame(show_example_points = TRUE)
```

## Important Notes

1. **Always use adjusted coordinates** for anatomical mapping to ensure accuracy
2. **Check `in_canvas`** to filter out gazes that fall on the borders
3. **Original coordinates preserved** in `gaze_x_original` and `gaze_y_original`
4. **DPR handling** is automatic - the 1.25 device pixel ratio is applied internally

## Troubleshooting

### Common Issues

1. **Too many out-of-bounds points**: Check if Tobii calibration includes full screen
2. **Unexpected mapping**: Verify display dimensions match your setup
3. **Edge effects**: Gazes near canvas edges may have reduced accuracy

### Validation

```r
# Check your setup
test_points <- data.frame(
  x = c(0.133, 0.5, 0.867),  # Left edge, center, right edge
  y = c(0.5, 0.5, 0.5)
)

for(i in 1:nrow(test_points)) {
  adj <- resolve_coordinates(test_points$x[i], test_points$y[i])
  cat(sprintf("Tobii [%.3f, %.3f] -> Canvas [%.3f, %.3f] %s\n",
              test_points$x[i], test_points$y[i],
              adj$x, adj$y,
              ifelse(adj$in_bounds, "✓", "✗")))
}
```

Expected output:
```
Tobii [0.133, 0.500] -> Canvas [0.000, 0.476] ✓
Tobii [0.500, 0.500] -> Canvas [0.500, 0.476] ✓
Tobii [0.867, 0.500] -> Canvas [1.000, 0.476] ✓
```
