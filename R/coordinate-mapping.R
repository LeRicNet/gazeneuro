#' Coordinate Mapping Functions for Gaze-Neuroimaging Integration
#'
#' Functions for mapping gaze coordinates to brain slice coordinates
#' and determining slice membership.

#' Locate a single gaze point on a brain slice
#'
#' @description
#' Maps a gaze coordinate (in normalized display space) to the corresponding
#' voxel location in the brain image, given the current slice being viewed.
#'
#' @param gaze_x Gaze X coordinate (0-1, normalized display width)
#' @param gaze_y Gaze Y coordinate (0-1, normalized display height)
#' @param slice_num Current slice number being viewed (1-based)
#' @param nifti_data List from preload_nifti_data() containing nvimage object
#' @param plane Imaging plane (AXIAL, SAGITTAL, or CORONAL)
#' @param display_bounds Optional list with xmin, xmax, ymin, ymax of display area
#'
#' @return Named list with:
#'   \item{voxel}{Voxel coordinates (i, j, k)}
#'   \item{mm}{World coordinates in mm (x, y, z)}
#'   \item{frac}{Fractional coordinates (0-1)}
#'   \item{valid}{Logical indicating if coordinates are within image bounds}
#'
#' @export
locate_single_point <- function(gaze_x, gaze_y, slice_num, nifti_data,
                                plane = "AXIAL",
                                display_bounds = NULL) {

  # Validate inputs
  if (is.null(nifti_data) || is.null(nifti_data$nvimage)) {
    stop("nifti_data must be a valid object from preload_nifti_data()")
  }

  nvimg <- nifti_data$nvimage
  dims <- nifti_data$dims

  # Default display bounds (full image)
  if (is.null(display_bounds)) {
    display_bounds <- list(xmin = 0, xmax = 1, ymin = 0, ymax = 1)
  }

  # Map display coordinates to image fractional coordinates
  # Account for potential display cropping/scaling
  img_frac_x <- (gaze_x - display_bounds$xmin) /
    (display_bounds$xmax - display_bounds$xmin)
  img_frac_y <- (gaze_y - display_bounds$ymin) /
    (display_bounds$ymax - display_bounds$ymin)

  # Clamp to [0, 1]
  img_frac_x <- max(0, min(1, img_frac_x))
  img_frac_y <- max(0, min(1, img_frac_y))

  # Convert slice number to fractional z coordinate (1-based to 0-1)
  # Use 0-based indexing internally, then normalize
  slice_0based <- slice_num - 1

  # Build fractional coordinate based on plane
  if (plane == "AXIAL") {
    # AXIAL: X-Y plane, slice along Z
    frac <- c(img_frac_x, img_frac_y, slice_0based / (dims[3] - 1))
    voxel <- c(
      round(img_frac_x * (dims[1] - 1)),
      round(img_frac_y * (dims[2] - 1)),
      slice_0based
    )
  } else if (plane == "SAGITTAL") {
    # SAGITTAL: Y-Z plane, slice along X
    frac <- c(slice_0based / (dims[1] - 1), img_frac_x, img_frac_y)
    voxel <- c(
      slice_0based,
      round(img_frac_x * (dims[2] - 1)),
      round(img_frac_y * (dims[3] - 1))
    )
  } else if (plane == "CORONAL") {
    # CORONAL: X-Z plane, slice along Y
    frac <- c(img_frac_x, slice_0based / (dims[2] - 1), img_frac_y)
    voxel <- c(
      round(img_frac_x * (dims[1] - 1)),
      slice_0based,
      round(img_frac_y * (dims[3] - 1))
    )
  } else {
    stop("plane must be AXIAL, SAGITTAL, or CORONAL")
  }

  # Convert to mm coordinates
  mm <- nvimg$convertFrac2MM(frac, isForceSliceMM = TRUE)[1:3]

  # Check if within bounds
  valid <- all(voxel >= 0) &&
    voxel[1] < dims[1] &&
    voxel[2] < dims[2] &&
    voxel[3] < dims[3]

  return(list(
    voxel = voxel + 1,  # Return 1-based voxel indices
    mm = mm,
    frac = frac,
    valid = valid
  ))
}


#' Get slice number from normalized index
#'
#' @description
#' Converts a normalized slice index (0-1) to an integer slice number (1-based).
#' Handles boundary cases and rounding consistently.
#'
#' @param normalized_index Slice index as fraction (0-1)
#' @param n_slices Total number of slices
#' @param method Rounding method: "round" (default), "floor", or "ceiling"
#'
#' @return Integer slice number (1-based)
#'
#' @details
#' The normalized index maps linearly to slice numbers:
#' - index = 0 → slice 1 (first slice)
#' - index = 1 → slice n_slices (last slice)
#' - index = 0.5 → middle slice
#'
#' @export
get_slice_number <- function(normalized_index, n_slices, method = "round") {

  # Validate inputs
  if (n_slices < 1) {
    stop("n_slices must be >= 1")
  }

  # Handle single-slice case

  if (n_slices == 1) {
    return(1L)
  }

  # Clamp normalized index to [0, 1]
  normalized_index <- max(0, min(1, normalized_index))

  # Map [0, 1] to [0, n_slices-1] (0-based)
  slice_0based <- normalized_index * (n_slices - 1)

  # Apply rounding method
  slice_0based <- switch(method,
                         "round" = round(slice_0based),
                         "floor" = floor(slice_0based),
                         "ceiling" = ceiling(slice_0based),
                         stop("method must be 'round', 'floor', or 'ceiling'")
  )

  # Convert to 1-based and ensure bounds
  slice_num <- as.integer(slice_0based + 1)
  slice_num <- max(1L, min(n_slices, slice_num))

  return(slice_num)
}


#' Get normalized index from slice number
#'
#' @description
#' Converts an integer slice number (1-based) to a normalized index (0-1).
#' Inverse of get_slice_number().
#'
#' @param slice_num Integer slice number (1-based)
#' @param n_slices Total number of slices
#'
#' @return Normalized index (0-1)
#'
#' @export
get_normalized_index <- function(slice_num, n_slices) {

  if (n_slices < 1) {
    stop("n_slices must be >= 1")
  }

  if (n_slices == 1) {
    return(0.5)
  }

  # Clamp slice_num to valid range
  slice_num <- max(1, min(n_slices, slice_num))

  # Map [1, n_slices] to [0, 1]
  normalized <- (slice_num - 1) / (n_slices - 1)

  return(normalized)
}


#' Map multiple gaze points to slice coordinates
#'
#' @description
#' Batch version of locate_single_point for processing multiple gaze points.
#'
#' @param gaze_data Data frame with gaze_x, gaze_y columns
#' @param slice_data Data frame with slice_num column (matching rows in gaze_data)
#' @param nifti_data List from preload_nifti_data()
#' @param plane Imaging plane
#'
#' @return Data frame with added voxel_i, voxel_j, voxel_k, mm_x, mm_y, mm_z columns
#'
#' @export
#' @importFrom dplyr mutate
map_gaze_to_voxels <- function(gaze_data, slice_data, nifti_data, plane = "AXIAL") {

  if (nrow(gaze_data) != nrow(slice_data)) {
    stop("gaze_data and slice_data must have the same number of rows")
  }

  # Initialize output columns
  n <- nrow(gaze_data)
  voxel_i <- integer(n)
  voxel_j <- integer(n)
  voxel_k <- integer(n)
  mm_x <- numeric(n)
  mm_y <- numeric(n)
  mm_z <- numeric(n)
  valid <- logical(n)

  # Process each point
  for (i in seq_len(n)) {
    result <- locate_single_point(
      gaze_x = gaze_data$gaze_x[i],
      gaze_y = gaze_data$gaze_y[i],
      slice_num = slice_data$slice_num[i],
      nifti_data = nifti_data,
      plane = plane
    )

    voxel_i[i] <- result$voxel[1]
    voxel_j[i] <- result$voxel[2]
    voxel_k[i] <- result$voxel[3]
    mm_x[i] <- result$mm[1]
    mm_y[i] <- result$mm[2]
    mm_z[i] <- result$mm[3]
    valid[i] <- result$valid
  }

  # Add to data frame
  gaze_data$voxel_i <- voxel_i
  gaze_data$voxel_j <- voxel_j
  gaze_data$voxel_k <- voxel_k
  gaze_data$mm_x <- mm_x
  gaze_data$mm_y <- mm_y
  gaze_data$mm_z <- mm_z
  gaze_data$valid <- valid

  return(gaze_data)
}


#' Calculate slice transition boundaries
#'
#' @description
#' Given z-axis events, calculates the time boundaries where slice transitions occur.
#'
#' @param z_axis_data Data frame with time_sec and slice_num columns
#'
#' @return Data frame with transition_time, slice_before, slice_after columns
#'
#' @export
get_slice_transitions <- function(z_axis_data) {

  if (nrow(z_axis_data) < 2) {
    return(data.frame(
      transition_time = numeric(0),
      slice_before = integer(0),
      slice_after = integer(0)
    ))
  }

  # Sort by time
  z_sorted <- z_axis_data[order(z_axis_data$time_sec), ]

  # Find transitions (where slice changes)
  n <- nrow(z_sorted)
  transitions <- data.frame(
    transition_time = z_sorted$time_sec[-1],
    slice_before = z_sorted$slice_num[-n],
    slice_after = z_sorted$slice_num[-1]
  )

  # Keep only actual transitions
  transitions <- transitions[transitions$slice_before != transitions$slice_after, ]

  return(transitions)
}
