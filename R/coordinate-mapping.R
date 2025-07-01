#' Map gaze coordinates to anatomical locations
#'
#' @description
#' Convert gaze tracking coordinates to anatomical locations in mm and voxel space.
#' This provides similar functionality to the JavaScript locate() method.
#'
#' @param integrated_data Data frame from integrate_all_gaze_points()
#' @param nifti_data List from preload_nifti_data()
#' @param row_indices Which rows to process (NULL for all)
#' @param include_values Whether to extract intensity values at locations
#' @param adjust_coordinates Whether to adjust coordinates for display frame (default TRUE)
#'
#' @return Data frame with anatomical locations for each gaze point
#' @export
#' @importFrom dplyr mutate select
#' @examples
#' \dontrun{
#' # Map all gaze points to anatomical locations
#' locations <- map_gaze_to_anatomy(integrated, nifti_data)
#'
#' # Map only specific rows
#' locations <- map_gaze_to_anatomy(integrated, nifti_data, row_indices = 1:100)
#' }
map_gaze_to_anatomy <- function(integrated_data, nifti_data,
                                row_indices = NULL,
                                include_values = TRUE,
                                adjust_coordinates = TRUE) {

  # Validate input
  required_cols <- c("gaze_x", "gaze_y", "slice_index", "time_sec", "time_aligned", "plane")
  missing_cols <- setdiff(required_cols, names(integrated_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in integrated_data: ", paste(missing_cols, collapse = ", "))
  }

  if (is.null(row_indices)) {
    row_indices <- 1:nrow(integrated_data)
  }

  # Get NVImage object
  nvimg <- nifti_data$nvimage
  img_data <- nifti_data$data
  dims <- nifti_data$dims

  # Process each gaze point
  locations <- lapply(row_indices, function(i) {
    row <- integrated_data[i, ]

    # Adjust coordinates if needed (for display frame)
    if (adjust_coordinates) {
      adjusted <- resolve_coordinates(row$gaze_x, row$gaze_y)
      gaze_x_adj <- adjusted$x
      gaze_y_adj <- adjusted$y
      in_bounds <- adjusted$in_bounds
    } else {
      gaze_x_adj <- row$gaze_x
      gaze_y_adj <- row$gaze_y
      in_bounds <- TRUE
    }

    # Create fractional coordinates
    frac <- c(gaze_x_adj, gaze_y_adj, row$slice_index)

    # Convert to mm coordinates
    mm_coords <- nvimg$convertFrac2MM(frac, isForceSliceMM = TRUE)
    mm <- mm_coords[1:3]

    # Convert to voxel coordinates
    voxel <- nvimg$mm2vox(mm, frac = FALSE)
    voxel_frac <- nvimg$mm2vox(mm, frac = TRUE)

    # Get actual slice number based on plane
    slice_num <- get_slice_number(row$slice_index, row$plane, dims)

    # Extract intensity value if requested
    value <- NA
    if (include_values && isTRUE(in_bounds)) {
      value <- safe_get_value(img_data, voxel, dims)
    }

    # Create location record
    list(
      # Original gaze data
      gaze_id = row$gaze_id,
      time_sec = row$time_sec,
      time_aligned = row$time_aligned,

      # Original coordinates (Tobii frame)
      gaze_x_original = row$gaze_x,
      gaze_y_original = row$gaze_y,

      # Adjusted coordinates (canvas)
      gaze_x = gaze_x_adj,
      gaze_y = gaze_y_adj,
      in_canvas = in_bounds,

      slice_index = row$slice_index,

      # Fractional coordinates (0-1)
      frac_x = frac[1],
      frac_y = frac[2],
      frac_z = frac[3],

      # MM coordinates
      mm_x = mm[1],
      mm_y = mm[2],
      mm_z = mm[3],

      # Voxel coordinates (integer)
      vox_x = voxel[1],
      vox_y = voxel[2],
      vox_z = voxel[3],

      # Voxel coordinates (fractional)
      vox_frac_x = voxel_frac[1],
      vox_frac_y = voxel_frac[2],
      vox_frac_z = voxel_frac[3],

      # Metadata
      plane = row$plane,
      slice_num = slice_num,
      image_id = row$image_id,
      intensity = value
    )
  })

  # Convert to data frame
  do.call(rbind.data.frame, locations)
}

#' Convert normalized slice index to actual slice number
#'
#' @param slice_index Normalized slice index (0-1)
#' @param plane Imaging plane (AXIAL, SAGITTAL, CORONAL)
#' @param dims Image dimensions
#' @return Slice number (1-based)
#' @export
get_slice_number <- function(slice_index, plane, dims) {
  if (plane == "AXIAL") {
    round(slice_index * (dims[3] - 1)) + 1
  } else if (plane == "SAGITTAL") {
    round(slice_index * (dims[1] - 1)) + 1
  } else if (plane == "CORONAL") {
    round(slice_index * (dims[2] - 1)) + 1
  } else {
    NA
  }
}

#' Safely extract intensity value from image data
#'
#' @param img_data 3D array of image data
#' @param voxel Voxel coordinates (1-based)
#' @param dims Image dimensions
#' @return Intensity value or NA if out of bounds
#' @export
safe_get_value <- function(img_data, voxel, dims) {
  # Check bounds
  if (any(voxel < 0) ||
      voxel[1] >= dims[1] ||
      voxel[2] >= dims[2] ||
      voxel[3] >= dims[3]) {
    return(NA)
  }

  # R uses 1-based indexing
  tryCatch({
    img_data[voxel[1] + 1, voxel[2] + 1, voxel[3] + 1]
  }, error = function(e) {
    NA
  })
}

#' Create anatomical location for a single gaze point
#'
#' @description
#' This is equivalent to the JavaScript locate() function for a single point.
#' Takes gaze coordinates and returns anatomical location information.
#'
#' @param x Gaze x coordinate (0-1 normalized to full display)
#' @param y Gaze y coordinate (0-1 normalized to full display)
#' @param z Normalized slice position (0-1)
#' @param nifti_data List from preload_nifti_data()
#' @param plane Imaging plane
#' @param timestamp_us Timestamp in microseconds
#' @param adjust_coordinates Whether to adjust for display frame (default TRUE)
#'
#' @return List with location information
#' @export
locate_single_point <- function(x, y, z, nifti_data,
                                plane = "AXIAL",
                                timestamp_us = NULL,
                                adjust_coordinates = TRUE) {

  nvimg <- nifti_data$nvimage
  img_data <- nifti_data$data
  dims <- nifti_data$dims

  # Adjust coordinates if needed
  if (adjust_coordinates) {
    adjusted <- resolve_coordinates(x, y)
    x_adj <- adjusted$x
    y_adj <- adjusted$y
    in_bounds <- adjusted$in_bounds
  } else {
    x_adj <- x
    y_adj <- y
    in_bounds <- TRUE
  }

  # Create fractional coordinates
  frac <- c(x_adj, y_adj, z)

  # Convert to mm
  mm_coords <- nvimg$convertFrac2MM(frac, isForceSliceMM = TRUE)
  mm <- mm_coords[1:3]

  # Convert to voxel
  voxel <- nvimg$mm2vox(mm, frac = FALSE)
  voxel_frac <- nvimg$mm2vox(mm, frac = TRUE)

  # Get intensity value
  value <- ifelse(in_bounds, safe_get_value(img_data, voxel, dims), NA)

  # Calculate slice visibility (simplified version)
  visibility <- ifelse(!is.na(value) && value > 0, 1, 0)

  # Create location object similar to JavaScript
  location <- list(
    mm = mm,
    vox = voxel,
    frac = frac,
    xy = c(x, y),  # Original coordinates
    xy_adjusted = c(x_adj, y_adj),  # Adjusted coordinates
    in_canvas = in_bounds,
    values = list(
      list(
        image_id = 'nifti',
        image_plane = plane,
        timestamp_us = timestamp_us %||% as.numeric(Sys.time()) * 1e6,
        filename = nifti_data$nifti@.xData$.sourcePath,
        value = value,
        tracking_mode = 'gaze',
        visibility = visibility,
        mm = mm,
        vox = voxel,
        mm_image = mm,
        vox_image = voxel
      )
    )
  )

  return(location)
}

#' Batch process gaze points to anatomical locations
#'
#' @description
#' Process multiple gaze points and create a summary of anatomical locations visited.
#'
#' @param integrated_data Data frame from integrate_all_gaze_points()
#' @param nifti_data List from preload_nifti_data()
#' @param summarize Whether to summarize by unique locations
#'
#' @return Data frame with anatomical locations
#' @export
#' @importFrom dplyr group_by summarise n distinct
batch_locate_gaze <- function(integrated_data, nifti_data, summarize = FALSE) {

  message("Mapping ", nrow(integrated_data), " gaze points to anatomical locations...")

  # Map all points
  locations <- map_gaze_to_anatomy(integrated_data, nifti_data)

  if (summarize) {
    # Summarize by unique voxel locations
    summary <- locations %>%
      group_by(vox_x, vox_y, vox_z, plane) %>%
      summarise(
        n_gazes = n(),
        duration = max(time_sec) - min(time_sec),
        mean_intensity = mean(intensity, na.rm = TRUE),
        first_visit = min(time_sec),
        last_visit = max(time_sec),
        mm_x = first(mm_x),
        mm_y = first(mm_y),
        mm_z = first(mm_z),
        .groups = "drop"
      ) %>%
      arrange(desc(n_gazes))

    message("Found ", nrow(summary), " unique voxel locations visited")
    return(summary)
  }

  return(locations)
}

#' Resolve coordinates with display frame adjustment
#'
#' @description
#' Converts Tobii gaze coordinates (0-1 normalized to full display) to
#' image canvas coordinates (0-1 normalized to actual image area).
#' Handles the frame borders and device pixel ratio.
#'
#' @param x X coordinate from Tobii (0-1 normalized to full display)
#' @param y Y coordinate from Tobii (0-1 normalized to full display)
#' @param frame_width Full frame width in pixels (default 2624)
#' @param frame_height Full frame height in pixels (default 1640)
#' @param canvas_width Image canvas width in pixels (default 1924)
#' @param canvas_height Image canvas height in pixels (default 1560)
#' @param left_border Left border width in pixels (default 350)
#' @param top_border Top border height in pixels (default 80)
#' @param dpr Device pixel ratio (default 1.25)
#' @return List with adjusted x and y coordinates (0-1 normalized to canvas)
#' @export
resolve_coordinates <- function(x, y,
                                frame_width = 2624,
                                frame_height = 1640,
                                canvas_width = 1924,
                                canvas_height = 1560,
                                left_border = 350,
                                top_border = 80,
                                dpr = 1.25) {

  # Convert from normalized (0-1) to pixel coordinates on full frame
  pixel_x <- x * frame_width * dpr
  pixel_y <- y * frame_height * dpr

  # Adjust for borders to get coordinates relative to canvas
  canvas_pixel_x <- pixel_x - (left_border * dpr)
  canvas_pixel_y <- pixel_y - (top_border * dpr)

  # Normalize to canvas dimensions (0-1)
  normalized_x <- canvas_pixel_x / (canvas_width * dpr)
  normalized_y <- canvas_pixel_y / (canvas_height * dpr)

  # Clamp to valid range [0, 1]
  normalized_x <- pmax(0, pmin(1, normalized_x))
  normalized_y <- pmax(0, pmin(1, normalized_y))

  list(
    x = normalized_x,
    y = normalized_y,
    in_bounds = (canvas_pixel_x >= 0 &&
                   canvas_pixel_x <= canvas_width * dpr &&
                   canvas_pixel_y >= 0 &&
                   canvas_pixel_y <= canvas_height * dpr)
  )
}

#' Create coordinate mapping summary plot
#'
#' @param locations Data frame from map_gaze_to_anatomy()
#' @param nifti_data List from preload_nifti_data()
#' @export
#' @importFrom graphics image points text
plot_coordinate_mapping <- function(locations, nifti_data) {

  # Get middle slice for visualization
  dims <- nifti_data$dims
  mid_slice <- ceiling(dims[3] / 2)

  # Filter locations for middle axial slice
  slice_locs <- locations %>%
    filter(plane == "AXIAL",
           abs(vox_z - mid_slice) < 2)

  if (nrow(slice_locs) == 0) {
    message("No locations found for middle slice")
    return(invisible(NULL))
  }

  # Plot the slice
  slice_data <- nifti_data$data[,, mid_slice]

  par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))

  # Left panel: Voxel space
  image(slice_data, col = gray((0:255)/255),
        main = "Voxel Space Mapping",
        xlab = "X voxel", ylab = "Y voxel")

  # Add gaze points in voxel space
  points(slice_locs$vox_x / dims[1],
         slice_locs$vox_y / dims[2],
         col = rgb(1, 0, 0, 0.5), pch = 19, cex = 0.8)

  # Right panel: MM space
  image(slice_data, col = gray((0:255)/255),
        main = "MM Space Mapping",
        xlab = "X mm", ylab = "Y mm")

  # Show mm coordinates
  text(0.5, 0.95,
       sprintf("MM range: X[%.1f, %.1f] Y[%.1f, %.1f]",
               min(slice_locs$mm_x), max(slice_locs$mm_x),
               min(slice_locs$mm_y), max(slice_locs$mm_y)),
       cex = 0.8)

  par(mfrow = c(1, 1))
}

# Helper function for NULL default
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Visualize display frame and canvas layout
#'
#' @description
#' Shows the relationship between the full display frame and the image canvas,
#' helping to understand the coordinate transformation.
#'
#' @param show_example_points Whether to show example gaze points
#' @export
#' @importFrom graphics rect text points legend
visualize_display_frame <- function(show_example_points = TRUE) {
  # Display parameters
  frame_width <- 2624
  frame_height <- 1640
  canvas_width <- 1924
  canvas_height <- 1560
  left_border <- 350
  top_border <- 80
  dpr <- 1.25

  # Set up plot
  par(mar = c(4, 4, 3, 2))
  plot(0, 0, type = "n",
       xlim = c(0, frame_width),
       ylim = c(0, frame_height),
       xlab = "X (pixels)",
       ylab = "Y (pixels)",
       main = "Display Frame and Image Canvas Layout",
       asp = 1)

  # Draw frame
  rect(0, 0, frame_width, frame_height,
       border = "black", lwd = 2)

  # Draw canvas
  rect(left_border, top_border,
       left_border + canvas_width,
       top_border + canvas_height,
       border = "blue", lwd = 2, lty = 2)

  # Add labels
  text(frame_width/2, frame_height - 20,
       sprintf("Full Frame: %d x %d px", frame_width, frame_height),
       cex = 1.2)

  text(left_border + canvas_width/2, top_border + canvas_height/2,
       sprintf("Image Canvas: %d x %d px", canvas_width, canvas_height),
       cex = 1.2, col = "blue")

  # Show borders
  text(left_border/2, frame_height/2,
       sprintf("Left: %dpx", left_border),
       srt = 90, cex = 0.8)

  text(frame_width - left_border/2, frame_height/2,
       sprintf("Right: %dpx", left_border),
       srt = 90, cex = 0.8)

  text(frame_width/2, top_border/2,
       sprintf("Top: %dpx", top_border),
       cex = 0.8)

  text(frame_width/2, frame_height - (frame_height - top_border - canvas_height)/2,
       sprintf("Bottom: %dpx", frame_height - top_border - canvas_height),
       cex = 0.8)

  if (show_example_points) {
    # Example gaze points in Tobii coordinates
    example_points <- data.frame(
      x = c(0.25, 0.5, 0.75, 0.1, 0.9),
      y = c(0.25, 0.5, 0.75, 0.5, 0.5),
      label = c("TL", "Center", "BR", "Left", "Right")
    )

    for (i in 1:nrow(example_points)) {
      # Original position (Tobii coordinates)
      orig_x <- example_points$x[i] * frame_width
      orig_y <- example_points$y[i] * frame_height

      # Adjusted position
      adj <- resolve_coordinates(example_points$x[i], example_points$y[i])
      adj_x <- left_border + adj$x * canvas_width
      adj_y <- top_border + adj$y * canvas_height

      # Plot original
      points(orig_x, orig_y, pch = 19, col = "red", cex = 1.5)

      # Plot adjusted if in bounds
      if (adj$in_bounds) {
        points(adj_x, adj_y, pch = 19, col = "green", cex = 1.5)
        # Connect with arrow
        arrows(orig_x, orig_y, adj_x, adj_y,
               length = 0.1, col = "gray", lty = 2)
      }

      # Label
      text(orig_x, orig_y - 30, example_points$label[i],
           cex = 0.8, col = "red")
    }

    legend("bottomright",
           legend = c("Tobii coordinates", "Canvas coordinates", "Out of bounds"),
           pch = c(19, 19, 19),
           col = c("red", "green", "gray"),
           cex = 0.8)
  }

  # Add DPR note
  text(10, 30, sprintf("Device Pixel Ratio: %.2f", dpr),
       adj = 0, cex = 0.8, font = 3)
}

#' Export gaze anatomical locations for external analysis
#'
#' @param locations Data frame from map_gaze_to_anatomy()
#' @param output_file Path to output file
#' @param format Output format ("csv", "json", or "nifti")
#' @param nifti_template Template NIfTI file (required for nifti format)
#' @export
#' @examples
#' \dontrun{
#' # Export as CSV
#' export_anatomical_locations(locations, "gaze_locations.csv")
#'
#' # Export as NIfTI heatmap
#' export_anatomical_locations(locations, "gaze_heatmap.nii.gz",
#'                           format = "nifti", nifti_template = nifti_data)
#' }
export_anatomical_locations <- function(locations, output_file,
                                        format = "csv",
                                        nifti_template = NULL) {

  if (format == "csv") {
    write.csv(locations, output_file, row.names = FALSE)
    message("Exported ", nrow(locations), " locations to ", output_file)

  } else if (format == "json") {
    # Convert to list format for JSON
    json_data <- lapply(1:nrow(locations), function(i) {
      as.list(locations[i, ])
    })
    jsonlite::write_json(json_data, output_file)
    message("Exported ", nrow(locations), " locations to ", output_file)

  } else if (format == "nifti") {
    if (is.null(nifti_template)) {
      stop("nifti_template required for NIfTI output")
    }

    # Create 3D heatmap
    dims <- dim(nifti_template$data)
    heatmap <- array(0, dim = dims)

    # Accumulate gaze counts
    for (i in 1:nrow(locations)) {
      vox <- c(locations$vox_x[i], locations$vox_y[i], locations$vox_z[i])
      if (all(vox >= 0) && all(vox < dims)) {
        heatmap[vox[1]+1, vox[2]+1, vox[3]+1] <-
          heatmap[vox[1]+1, vox[2]+1, vox[3]+1] + 1
      }
    }

    # Write using RNifti
    RNifti::writeNifti(heatmap, output_file, template = nifti_template$nifti)
    message("Exported gaze heatmap to ", output_file)
  } else {
    stop("Unsupported format: ", format)
  }
}
