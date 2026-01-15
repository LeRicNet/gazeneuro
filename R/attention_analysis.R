#' Attention Analysis Functions for gazeneuro (Base R version)
#'
#' @description
#' Functions for generating attention heatmaps and analyzing gaze patterns.
#' Includes fixation detection, 2D slice heatmaps, and 3D volumetric attention maps.
#' This version uses only base R functions for maximum compatibility.

# =============================================================================
# Fixation Detection
# =============================================================================

#' Detect fixations using velocity threshold
#'
#' @description
#' Identifies fixation periods in gaze data using a velocity-based threshold.
#' Samples below the threshold are classified as fixations; above as saccades.
#'
#' @param integrated_data Data frame from integrate_all_gaze_points()
#' @param velocity_threshold Velocity threshold in normalized units/second (default: 0.1)
#' @param min_fixation_duration Minimum fixation duration in seconds (default: 0.1)
#' @param screen_diagonal_deg Visual angle of screen diagonal in degrees (default: 30)
#'
#' @return Data frame with added columns for velocity, is_fixation, fixation_id, fixation_duration
#'
#' @export
detect_fixations <- function(integrated_data,
                             velocity_threshold = 0.5,
                             min_fixation_duration = 0.05,
                             screen_diagonal_deg = 30) {

  # Sort by time
  data <- integrated_data[order(integrated_data$time_sec), ]
  n <- nrow(data)

  # Calculate velocity (normalized units per second)
  dt <- c(0, diff(data$time_sec))
  dx <- c(0, diff(data$gaze_x))
  dy <- c(0, diff(data$gaze_y))

  # Avoid division by zero
  dt[dt < 1e-6] <- 1e-6

  data$velocity <- sqrt(dx^2 + dy^2) / dt
  data$velocity_deg <- data$velocity * screen_diagonal_deg

  # Classify as fixation or saccade
  data$is_fixation <- data$velocity < velocity_threshold

  # Assign fixation IDs (contiguous fixation periods)
  data$fixation_id <- NA_integer_
  current_id <- 0
  in_fixation <- FALSE

  for (i in 1:n) {
    if (data$is_fixation[i]) {
      if (!in_fixation) {
        current_id <- current_id + 1
        in_fixation <- TRUE
      }
      data$fixation_id[i] <- current_id
    } else {
      in_fixation <- FALSE
    }
  }

  # Calculate fixation durations and filter by minimum
  data$fixation_duration <- NA_real_

  if (current_id > 0) {
    for (fid in 1:current_id) {
      idx <- which(data$fixation_id == fid)
      if (length(idx) > 0) {
        duration <- max(data$time_sec[idx]) - min(data$time_sec[idx])
        if (duration >= min_fixation_duration) {
          data$fixation_duration[idx] <- duration
        } else {
          # Invalidate short fixations
          data$is_fixation[idx] <- FALSE
          data$fixation_id[idx] <- NA_integer_
        }
      }
    }
  }

  # Summary message
  n_fixations <- length(unique(na.omit(data$fixation_id)))
  pct_fixation <- 100 * mean(data$is_fixation, na.rm = TRUE)

  message(sprintf("Fixation detection: %d fixations identified (%.1f%% of samples)",
                  n_fixations, pct_fixation))

  return(data)
}

# =============================================================================
# Attention Weights
# =============================================================================

#' Calculate attention weights for gaze samples
#'
#' @description
#' Computes weights for each gaze sample based on fixation duration and
#' velocity-mediated uncertainty.
#'
#' @param data Data frame with fixation detection (from detect_fixations())
#' @param weight_by_duration Weight samples by fixation duration
#' @param uncertainty_weighted Down-weight high-velocity samples
#' @param velocity_halfweight Velocity at which weight is halved (default: 0.05)
#'
#' @return Data frame with added 'attention_weight' column
#'
#' @export
calculate_attention_weights <- function(data,
                                        weight_by_duration = TRUE,
                                        uncertainty_weighted = TRUE,
                                        velocity_halfweight = 0.05) {

  n <- nrow(data)
  data$attention_weight <- rep(1.0, n)

  # Duration weighting: longer fixations contribute more
  if (weight_by_duration) {
    duration_weight <- ifelse(data$is_fixation & !is.na(data$fixation_duration),
                              data$fixation_duration, 0)
    mean_dur <- mean(duration_weight, na.rm = TRUE)
    if (mean_dur > 1e-6) {
      duration_weight <- duration_weight / mean_dur
    }
    data$attention_weight <- data$attention_weight * duration_weight
  }

  # Uncertainty weighting: high velocity = more uncertainty = lower weight
  if (uncertainty_weighted) {
    uncertainty_weight <- 1 / (1 + data$velocity / velocity_halfweight)
    data$attention_weight <- data$attention_weight * uncertainty_weight
  }

  # Normalize weights to sum to number of samples
  total_weight <- sum(data$attention_weight, na.rm = TRUE)
  if (total_weight > 0) {
    data$attention_weight <- data$attention_weight * n / total_weight
  }

  return(data)
}

# =============================================================================
# Weighted 2D KDE
# =============================================================================

#' Weighted 2D kernel density estimation
#'
#' @description
#' Computes a weighted kernel density estimate on a 2D grid.
#'
#' @param x X coordinates
#' @param y Y coordinates
#' @param weights Sample weights (default: equal weights)
#' @param n Grid resolution (scalar or length-2 vector)
#' @param bandwidth Kernel bandwidth (scalar or length-2 vector)
#' @param xlim X limits
#' @param ylim Y limits
#'
#' @return List with x, y coordinates and z density matrix
#'
#' @export
weighted_kde2d <- function(x, y, weights = NULL, n = 50,
                           bandwidth = NULL, xlim = NULL, ylim = NULL) {

  # Handle inputs
  if (is.null(weights)) {
    weights <- rep(1, length(x))
  }

  if (length(n) == 1) {
    n <- c(n, n)
  }

  if (is.null(xlim)) {
    xlim <- range(x, na.rm = TRUE)
  }
  if (is.null(ylim)) {
    ylim <- range(y, na.rm = TRUE)
  }

  # Bandwidth selection (Silverman's rule if not specified)
  if (is.null(bandwidth)) {
    bandwidth <- c(
      1.06 * sd(x, na.rm = TRUE) * length(x)^(-1/5),
      1.06 * sd(y, na.rm = TRUE) * length(y)^(-1/5)
    )
  }
  if (length(bandwidth) == 1) {
    bandwidth <- c(bandwidth, bandwidth)
  }

  # Create output grid
  gx <- seq(xlim[1], xlim[2], length.out = n[1])
  gy <- seq(ylim[1], ylim[2], length.out = n[2])

  # Initialize density matrix
  z <- matrix(0, nrow = n[1], ncol = n[2])

  # Compute weighted KDE using Gaussian kernel
  for (i in seq_along(x)) {
    if (!is.na(weights[i]) && weights[i] > 0 && !is.na(x[i]) && !is.na(y[i])) {
      kx <- dnorm(gx, mean = x[i], sd = bandwidth[1])
      ky <- dnorm(gy, mean = y[i], sd = bandwidth[2])
      z <- z + weights[i] * outer(kx, ky)
    }
  }

  # Normalize by total weight
  valid <- !is.na(weights) & weights > 0
  total_weight <- sum(weights[valid])
  if (total_weight > 0) {
    z <- z / total_weight
  }

  return(list(x = gx, y = gy, z = z))
}

# =============================================================================
# 2D Slice Heatmaps
# =============================================================================

#' Generate attention heatmap for a single slice
#'
#' @param integrated_data Data frame from integrate_all_gaze_points()
#' @param nifti_data List from preload_nifti_data()
#' @param slice_num Slice number (1-based)
#' @param plane Imaging plane (AXIAL, SAGITTAL, or CORONAL)
#' @param bandwidth Kernel bandwidth in voxels (default: 5)
#' @param grid_resolution Output resolution; NULL uses NIfTI dimensions
#' @param weight_by_duration Weight by fixation duration
#' @param uncertainty_weighted Down-weight uncertain samples
#' @param normalize Normalize output to [0, 1]
#'
#' @return List containing heatmap matrix and metadata
#'
#' @export
generate_slice_heatmap <- function(integrated_data,
                                   nifti_data,
                                   slice_num,
                                   plane = "AXIAL",
                                   bandwidth = 5,
                                   grid_resolution = NULL,
                                   weight_by_duration = TRUE,
                                   uncertainty_weighted = TRUE,
                                   normalize = TRUE) {

  dims <- nifti_data$dims

  # Determine slice dimensions based on plane
  if (plane == "AXIAL") {
    slice_dims <- c(dims[1], dims[2])
    max_slice <- dims[3]
  } else if (plane == "SAGITTAL") {
    slice_dims <- c(dims[2], dims[3])
    max_slice <- dims[1]
  } else if (plane == "CORONAL") {
    slice_dims <- c(dims[1], dims[3])
    max_slice <- dims[2]
  } else {
    stop("plane must be AXIAL, SAGITTAL, or CORONAL")
  }

  # Validate slice number
  if (slice_num < 1 || slice_num > max_slice) {
    stop(sprintf("slice_num must be between 1 and %d for %s plane", max_slice, plane))
  }

  # Get gaze data for this slice
  gaze_for_slice <- get_all_gaze_for_slice(integrated_data, slice_num, plane, dims)

  if (nrow(gaze_for_slice) == 0) {
    message(sprintf("No gaze data for slice %d (%s)", slice_num, plane))
    return(list(
      heatmap = matrix(0, nrow = slice_dims[1], ncol = slice_dims[2]),
      x_coords = seq(0, 1, length.out = slice_dims[1]),
      y_coords = seq(0, 1, length.out = slice_dims[2]),
      slice_num = slice_num,
      plane = plane,
      n_samples = 0
    ))
  }

  # Detect fixations and calculate weights
  gaze_weighted <- detect_fixations(gaze_for_slice)
  gaze_weighted <- calculate_attention_weights(
    gaze_weighted,
    weight_by_duration = weight_by_duration,
    uncertainty_weighted = uncertainty_weighted
  )

  # Set grid resolution
  if (is.null(grid_resolution)) {
    grid_resolution <- slice_dims
  }

  # Convert normalized gaze coordinates to voxel coordinates
  gaze_weighted$voxel_x <- gaze_weighted$gaze_x * (slice_dims[1] - 1)
  gaze_weighted$voxel_y <- gaze_weighted$gaze_y * (slice_dims[2] - 1)

  # Generate weighted kernel density estimate
  heatmap <- weighted_kde2d(
    x = gaze_weighted$voxel_x,
    y = gaze_weighted$voxel_y,
    weights = gaze_weighted$attention_weight,
    n = grid_resolution,
    bandwidth = bandwidth,
    xlim = c(0, slice_dims[1] - 1),
    ylim = c(0, slice_dims[2] - 1)
  )

  # Normalize to [0, 1] if requested
  if (normalize && max(heatmap$z) > 0) {
    heatmap$z <- heatmap$z / max(heatmap$z)
  }

  message(sprintf("Generated heatmap for slice %d (%s): %d samples, %.1f%% fixations",
                  slice_num, plane, nrow(gaze_weighted),
                  100 * mean(gaze_weighted$is_fixation, na.rm = TRUE)))

  return(list(
    heatmap = heatmap$z,
    x_coords = heatmap$x,
    y_coords = heatmap$y,
    slice_num = slice_num,
    plane = plane,
    n_samples = nrow(gaze_weighted),
    pct_fixation = 100 * mean(gaze_weighted$is_fixation, na.rm = TRUE)
  ))
}

# =============================================================================
# 3D Volume Heatmaps
# =============================================================================

#' Generate 3D attention volume
#'
#' @param integrated_data Data frame from integrate_all_gaze_points()
#' @param nifti_data List from preload_nifti_data()
#' @param plane Imaging plane (AXIAL, SAGITTAL, or CORONAL)
#' @param bandwidth Kernel bandwidth in voxels (default: 5)
#' @param weight_by_duration Weight by fixation duration
#' @param uncertainty_weighted Down-weight uncertain samples
#' @param normalize Normalize output to [0, 1]
#' @param verbose Print progress messages
#'
#' @return List containing 3D volume array and metadata
#'
#' @export
generate_volume_heatmap <- function(integrated_data,
                                    nifti_data,
                                    plane = "AXIAL",
                                    bandwidth = 5,
                                    weight_by_duration = TRUE,
                                    uncertainty_weighted = TRUE,
                                    normalize = TRUE,
                                    verbose = TRUE) {

  dims <- nifti_data$dims

  # Initialize volume
  attention_volume <- array(0, dim = dims)

  # Determine slice count based on plane
  if (plane == "AXIAL") {
    n_slices <- dims[3]
  } else if (plane == "SAGITTAL") {
    n_slices <- dims[1]
  } else if (plane == "CORONAL") {
    n_slices <- dims[2]
  } else {
    stop("plane must be AXIAL, SAGITTAL, or CORONAL")
  }

  # Find which slices have gaze data
  plane_data <- integrated_data[integrated_data$plane == plane, ]

  if (plane == "AXIAL") {
    slice_nums <- round(plane_data$slice_index * (dims[3] - 1)) + 1
  } else if (plane == "SAGITTAL") {
    slice_nums <- round(plane_data$slice_index * (dims[1] - 1)) + 1
  } else {
    slice_nums <- round(plane_data$slice_index * (dims[2] - 1)) + 1
  }

  slices_with_data <- sort(unique(slice_nums))

  if (verbose) {
    message(sprintf("Generating volume heatmap: %d slices with gaze data",
                    length(slices_with_data)))
  }

  # Process each slice
  n_slices_processed <- 0
  for (slice_num in slices_with_data) {

    if (verbose && n_slices_processed %% 5 == 0) {
      message(sprintf("  Processing slice %d/%d...",
                      n_slices_processed + 1, length(slices_with_data)))
    }

    # Generate slice heatmap
    slice_result <- generate_slice_heatmap(
      integrated_data = integrated_data,
      nifti_data = nifti_data,
      slice_num = slice_num,
      plane = plane,
      bandwidth = bandwidth,
      weight_by_duration = weight_by_duration,
      uncertainty_weighted = uncertainty_weighted,
      normalize = FALSE
    )

    # Insert into volume at correct position
    if (slice_result$n_samples > 0) {
      if (plane == "AXIAL") {
        attention_volume[, , slice_num] <- slice_result$heatmap
      } else if (plane == "SAGITTAL") {
        attention_volume[slice_num, , ] <- slice_result$heatmap
      } else if (plane == "CORONAL") {
        attention_volume[, slice_num, ] <- slice_result$heatmap
      }
      n_slices_processed <- n_slices_processed + 1
    }
  }

  # Normalize entire volume to [0, 1]
  if (normalize && max(attention_volume) > 0) {
    attention_volume <- attention_volume / max(attention_volume)
  }

  if (verbose) {
    message(sprintf("Volume heatmap complete: %d slices processed, max attention = %.3f",
                    n_slices_processed, max(attention_volume)))
  }

  return(list(
    volume = attention_volume,
    affine = nifti_data$nvimage$matRAS,
    dims = dims,
    plane = plane,
    n_slices_with_data = n_slices_processed
  ))
}

# =============================================================================
# Export Functions
# =============================================================================

#' Export attention data as CSV
#'
#' @param volume_result Result from generate_volume_heatmap()
#' @param output_path Path for output CSV file
#' @param threshold Minimum attention value to include (default: 0.01)
#'
#' @return Invisible path to created file
#'
#' @export
export_attention_csv <- function(volume_result, output_path, threshold = 0.01) {

  vol <- volume_result$volume
  dims <- volume_result$dims

  # Get indices of voxels above threshold
  above_threshold <- which(vol > threshold, arr.ind = TRUE)

  if (nrow(above_threshold) == 0) {
    message("No voxels above threshold - no file created")
    return(invisible(NULL))
  }

  # Build output data frame
  output_df <- data.frame(
    voxel_x = above_threshold[, 1] - 1,
    voxel_y = above_threshold[, 2] - 1,
    voxel_z = above_threshold[, 3] - 1,
    attention = vol[above_threshold]
  )

  # Sort by attention (descending)
  output_df <- output_df[order(-output_df$attention), ]

  # Write CSV
  write.csv(output_df, output_path, row.names = FALSE)

  message(sprintf("Exported %d voxels to: %s", nrow(output_df), output_path))

  return(invisible(output_path))
}

# =============================================================================
# Visualization
# =============================================================================

#' Plot attention heatmap overlaid on brain slice
#'
#' @param slice_heatmap Result from generate_slice_heatmap()
#' @param nifti_data List from preload_nifti_data()
#' @param alpha Transparency of heatmap overlay (0-1)
#' @param show_colorbar Display colorbar
#'
#' @export
plot_attention_overlay <- function(slice_heatmap,
                                   nifti_data,
                                   alpha = 0.6,
                                   show_colorbar = TRUE) {

  slice_num <- slice_heatmap$slice_num
  plane <- slice_heatmap$plane
  heatmap <- slice_heatmap$heatmap

  # Get anatomical slice
  if (plane == "AXIAL") {
    anat_slice <- nifti_data$data[, , slice_num]
  } else if (plane == "SAGITTAL") {
    anat_slice <- nifti_data$data[slice_num, , ]
  } else {
    anat_slice <- nifti_data$data[, slice_num, ]
  }

  # Set up color palette
  colormap <- colorRampPalette(c("transparent", "yellow", "orange", "red"))(256)

  # Set up plot layout
  if (show_colorbar) {
    layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))
  }

  par(mar = c(4, 4, 3, 1))

  # Plot anatomical background
  image(anat_slice,
        col = gray((0:255)/255),
        main = sprintf("%s Slice %d - Attention Heatmap\n(%d samples)",
                       plane, slice_num, slice_heatmap$n_samples),
        xlab = "X", ylab = "Y",
        axes = FALSE)

  # Overlay attention heatmap
  if (max(heatmap) > 0) {
    heatmap_masked <- heatmap
    heatmap_masked[heatmap < 0.05 * max(heatmap)] <- NA

    image(heatmap_masked,
          col = colormap,
          add = TRUE)
  }

  # Add colorbar if requested
  if (show_colorbar) {
    par(mar = c(4, 1, 3, 3))
    colorbar_vals <- seq(0, 1, length.out = 256)
    image(1, colorbar_vals, t(as.matrix(colorbar_vals)),
          col = colormap,
          xaxt = "n", yaxt = "n",
          xlab = "", ylab = "")
    axis(4, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"))
    mtext("Attention", side = 4, line = 2)
  }

  # Reset layout
  layout(1)
  par(mar = c(5.1, 4.1, 4.1, 2.1))
}
