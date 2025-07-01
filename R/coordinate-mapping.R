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
                                include_values = TRUE) {

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

    # Create fractional coordinates
    # gaze_x and gaze_y are already 0-1
    # slice_index is the normalized z coordinate
    frac <- c(row$gaze_x, row$gaze_y, row$slice_index)

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
    if (include_values) {
      value <- safe_get_value(img_data, voxel, dims)
    }

    # Create location record
    list(
      # Original gaze data
      gaze_id = row$gaze_id,
      time_sec = row$time_sec,
      time_aligned = row$time_aligned,

      # Coordinates
      gaze_x = row$gaze_x,
      gaze_y = row$gaze_y,
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
#' @param x Gaze x coordinate (0-1)
#' @param y Gaze y coordinate (0-1)
#' @param z Normalized slice position (0-1)
#' @param nifti_data List from preload_nifti_data()
#' @param plane Imaging plane
#' @param timestamp_us Timestamp in microseconds
#'
#' @return List with location information
#' @export
locate_single_point <- function(x, y, z, nifti_data,
                                plane = "AXIAL",
                                timestamp_us = NULL) {

  nvimg <- nifti_data$nvimage
  img_data <- nifti_data$data
  dims <- nifti_data$dims

  # Create fractional coordinates
  frac <- c(x, y, z)

  # Convert to mm
  mm_coords <- nvimg$convertFrac2MM(frac, isForceSliceMM = TRUE)
  mm <- mm_coords[1:3]

  # Convert to voxel
  voxel <- nvimg$mm2vox(mm, frac = FALSE)
  voxel_frac <- nvimg$mm2vox(mm, frac = TRUE)

  # Get intensity value
  value <- safe_get_value(img_data, voxel, dims)

  # Calculate slice visibility (simplified version)
  # In the JS code, this extracts a slice and calculates max value
  # Here we'll check if the current voxel is non-zero
  visibility <- ifelse(!is.na(value) && value > 0, 1, 0)

  # Create location object similar to JavaScript
  location <- list(
    mm = mm,
    vox = voxel,
    frac = frac,
    xy = c(x, y),
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

#' Resolve coordinates with device pixel ratio adjustment
#'
#' @description
#' R equivalent of resolveCoordinates() - adjusts coordinates for display.
#' In R context, this is mainly for ensuring coordinates are properly scaled.
#'
#' @param x X coordinate
#' @param y Y coordinate
#' @param dpr Device pixel ratio (default 1)
#' @return List with adjusted x and y coordinates
#' @export
resolve_coordinates <- function(x, y, dpr = 1) {
  list(
    x = x * dpr,
    y = y * dpr
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
