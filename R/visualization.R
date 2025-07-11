#' Plot NIfTI slice with synchronized gaze data
#'
#' @description
#' Displays a NIfTI slice overlaid with gaze tracking data that was
#' recorded while viewing that slice. Uses pre-adjusted coordinates
#' from integrate_all_gaze_points() if available.
#'
#' @param nifti_data List from preload_nifti_data()
#' @param integrated_data Data frame from integrate_all_gaze_points()
#' @param slice_num Slice number to display
#' @param plane Imaging plane (AXIAL, SAGITTAL, or CORONAL)
#' @param show_density Whether to show density contours for many points
#'
#' @return Data frame of gaze points for the displayed slice
#' @export
#' @importFrom graphics image points lines legend abline text contour
#' @importFrom MASS kde2d
plot_slice_with_all_gaze <- function(nifti_data, integrated_data,
                                     slice_num = 12, plane = "AXIAL",
                                     show_density = TRUE) {

  # Get all gaze points for this slice
  gaze_points <- get_all_gaze_for_slice(integrated_data, slice_num, plane,
                                        dims = dim(nifti_data$data))

  message(sprintf("\nSlice %d: %d gaze points over %.1f seconds",
                  slice_num, nrow(gaze_points),
                  ifelse(nrow(gaze_points) > 0,
                         max(gaze_points$time_aligned) - min(gaze_points$time_aligned),
                         0)))

  # Get slice data
  if (plane == "AXIAL") {
    slice_data <- nifti_data$data[,, slice_num]
  } else if (plane == "SAGITTAL") {
    slice_data <- nifti_data$data[slice_num, , ]
  } else {
    slice_data <- nifti_data$data[, slice_num, ]
  }

  # Create plot
  par(mar = c(4, 4, 3, 2))

  # Base image
  image(slice_data, col = gray((0:255)/255),
        main = sprintf("%s Slice %d - %d gaze points",
                       plane, slice_num, nrow(gaze_points)))

  if (nrow(gaze_points) > 0) {
    # Use adjusted coordinates if available, otherwise use original
    if ("gaze_x_adjusted" %in% names(gaze_points)) {
      gaze_x_plot <- gaze_points$gaze_x_adjusted
      gaze_y_plot <- gaze_points$gaze_y_adjusted
    } else {
      gaze_x_plot <- gaze_points$gaze_x
      gaze_y_plot <- gaze_points$gaze_y
    }

    if (show_density && nrow(gaze_points) > 50) {
      # Create density overlay
      kde <- MASS::kde2d(gaze_x_plot, gaze_y_plot,
                         n = 50, lims = c(0, 1, 0, 1))
      contour(kde, add = TRUE, col = heat.colors(10), lwd = 2)
    } else {
      # Plot individual points
      points(gaze_x_plot, gaze_y_plot,
             col = rgb(1, 0, 0, 0.1), pch = 19, cex = 0.5)
    }

    # Show temporal progression with color
    n_points <- nrow(gaze_points)
    if (n_points > 1 && n_points < 500) {
      colors <- colorRampPalette(c("green", "yellow", "red"))(n_points)
      points(gaze_x_plot, gaze_y_plot,
             col = colors, pch = 19, cex = 0.7)
    }

    # Add mean position
    mean_x <- mean(gaze_x_plot)
    mean_y <- mean(gaze_y_plot)
    points(mean_x, mean_y, col = "blue", pch = 3, cex = 2, lwd = 3)

    # Add stats to plot
    text(0.05, 0.95, sprintf("n = %d", nrow(gaze_points)),
         adj = c(0, 1), cex = 0.8)
    text(0.05, 0.90, sprintf("%.1f sec",
                             max(gaze_points$time_aligned) - min(gaze_points$time_aligned)),
         adj = c(0, 1), cex = 0.8)
  }

  return(invisible(gaze_points))
}

#' Plot all viewed slices in a grid with gaze data
#'
#' @description
#' Creates a grid visualization showing all slices that were viewed,
#' each overlaid with its corresponding gaze data. Uses pre-adjusted
#' coordinates from integrate_all_gaze_points() if available.
#'
#' @param nifti_data List from preload_nifti_data()
#' @param integrated_data Data frame from integrate_all_gaze_points()
#' @param plane Imaging plane to display
#' @param max_slices Maximum number of slices to display
#' @param ncol Number of columns in grid (NULL for auto)
#'
#' @return Data frame with summary statistics per slice
#' @export
#' @importFrom graphics par image points lines box mtext
#' @importFrom dplyr filter mutate pull unique case_when
plot_all_slices_with_gaze <- function(nifti_data, integrated_data,
                                      plane = "AXIAL",
                                      max_slices = 25,
                                      ncol = NULL) {

  # Get unique slices that were viewed
  viewed_slices <- integrated_data %>%
    filter(plane == !!plane) %>%
    mutate(
      # Convert normalized index to slice number
      slice_num = case_when(
        plane == "AXIAL" ~ round(slice_index * (nifti_data$nvimage$dimsRAS[4] - 1)) + 1,
        plane == "SAGITTAL" ~ round(slice_index * (nifti_data$nvimage$dimsRAS[2] - 1)) + 1,
        plane == "CORONAL" ~ round(slice_index * (nifti_data$nvimage$dimsRAS[3] - 1)) + 1,
        TRUE ~ NA_real_
      )
    ) %>%
    filter(!is.na(slice_num)) %>%
    pull(slice_num) %>%
    unique() %>%
    sort()

  n_slices <- length(viewed_slices)
  message(sprintf("\nFound %d %s slices with gaze data", n_slices, plane))

  if (n_slices == 0) {
    message("No slices found for this plane")
    return(invisible(NULL))
  }

  # Limit number of slices if requested
  if (n_slices > max_slices) {
    message(sprintf("Limiting to first %d slices", max_slices))
    viewed_slices <- viewed_slices[1:max_slices]
    n_slices <- max_slices
  }

  # Calculate grid dimensions
  if (is.null(ncol)) {
    ncol <- ceiling(sqrt(n_slices))
  }
  nrow <- ceiling(n_slices / ncol)

  # Set up the plotting grid
  par(mfrow = c(nrow, ncol), mar = c(1, 1, 2, 1), oma = c(2, 2, 4, 2))

  # Plot each slice
  for (i in 1:n_slices) {
    slice_num <- viewed_slices[i]

    # Get slice data
    if (plane == "AXIAL") {
      slice_data <- nifti_data$data[,, slice_num]
    } else if (plane == "SAGITTAL") {
      slice_data <- nifti_data$data[slice_num, , ]
    } else if (plane == "CORONAL") {
      slice_data <- nifti_data$data[, slice_num, ]
    }

    # Get gaze points for this slice
    gaze_points <- get_all_gaze_for_slice(integrated_data, slice_num, plane,
                                          dims = dim(nifti_data$data))

    # Plot the slice
    image(slice_data,
          col = gray((0:255)/255),
          axes = FALSE,
          main = sprintf("Slice %d (n=%d)", slice_num, nrow(gaze_points)),
          cex.main = 0.8)

    # Add gaze points if any
    if (nrow(gaze_points) > 0) {
      # Use adjusted coordinates if available
      if ("gaze_x_adjusted" %in% names(gaze_points)) {
        gaze_x_plot <- gaze_points$gaze_x_adjusted
        gaze_y_plot <- gaze_points$gaze_y_adjusted
      } else {
        gaze_x_plot <- gaze_points$gaze_x
        gaze_y_plot <- gaze_points$gaze_y
      }

      # Add all gaze points
      points(gaze_x_plot,
             gaze_y_plot,
             col = rgb(1, 0, 0, 0.4),
             pch = 19,
             cex = 0.3)

      # Add trajectory lines if more than 5 points
      if (nrow(gaze_points) > 5) {
        lines(gaze_x_plot,
              gaze_y_plot,
              col = rgb(1, 0, 0, 0.2),
              lwd = 1)
      }

      # Mark center of gaze
      mean_x <- mean(gaze_x_plot)
      mean_y <- mean(gaze_y_plot)
      points(mean_x, mean_y,
             col = "blue", pch = 3, cex = 1, lwd = 2)
    }

    # Add slice border
    box()
  }

  # Add empty plots if needed to fill grid
  if (n_slices < nrow * ncol) {
    for (i in (n_slices + 1):(nrow * ncol)) {
      plot.new()
    }
  }

  # Add overall title
  mtext(sprintf("%s Slices with Gaze Overlay", plane),
        outer = TRUE, cex = 1.2, line = 2)
  mtext(sprintf("Red: gaze points, Blue cross: mean gaze position"),
        outer = TRUE, cex = 0.8, line = 0.5)

  # Reset par
  par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0))

  # Return summary statistics
  summary_stats <- analyze_gaze_by_slice(integrated_data, nifti_data, use_adjusted = TRUE)

  return(invisible(summary_stats))
}

#' Create gaze animation frames
#'
#' @description
#' Generates PNG frames showing gaze movement over time on NIfTI slices.
#' Frames can be combined into a video using ffmpeg. Uses pre-adjusted
#' coordinates from integrate_all_gaze_points() if available.
#'
#' @param nifti_data List from preload_nifti_data()
#' @param integrated_data Data frame from integrate_all_gaze_points()
#' @param output_dir Directory to save animation frames
#'
#' @export
#' @importFrom grDevices png dev.off
create_gaze_animation <- function(nifti_data, integrated_data,
                                  output_dir = "gaze_animation") {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Determine which coordinates to use
  use_adjusted <- "gaze_x_adjusted" %in% names(integrated_data)

  # Get unique time points where slices changed
  slice_changes <- integrated_data

  # Filter to within bounds if using adjusted coordinates
  if (use_adjusted && "within_bounds" %in% names(slice_changes)) {
    slice_changes <- slice_changes %>%
      filter(within_bounds)
  }

  slice_changes <- slice_changes %>%
    group_by(time_sec, plane, slice_index) %>%
    summarise(
      gaze_x = if (use_adjusted) mean(gaze_x_adjusted) else mean(gaze_x),
      gaze_y = if (use_adjusted) mean(gaze_y_adjusted) else mean(gaze_y),
      .groups = "drop"
    ) %>%
    arrange(time_sec)

  message(sprintf("\nCreating animation frames for %d time points...",
                  nrow(slice_changes)))

  # Create frames
  for (i in 1:nrow(slice_changes)) {
    row <- slice_changes[i,]

    # Determine slice number from normalized index
    if (row$plane == "AXIAL") {
      slice_num <- round(row$slice_index * 24) + 1  # For 25 slices
      slice_data <- nifti_data$data[,, slice_num]
    } else if (row$plane == "SAGITTAL") {
      slice_num <- round(row$slice_index * 511) + 1
      slice_data <- nifti_data$data[slice_num, , ]
    } else {
      slice_num <- round(row$slice_index * 511) + 1
      slice_data <- nifti_data$data[, slice_num, ]
    }

    # Create filename
    filename <- sprintf("%s/frame_%04d_%.2fs_%s_slice%02d.png",
                        output_dir, i, row$time_sec, row$plane, slice_num)

    png(filename, width = 600, height = 600)

    # Plot slice
    par(mar = c(2, 2, 3, 1))
    image(slice_data, col = gray((0:255)/255),
          main = sprintf("Time: %.1fs - %s Slice %d",
                         row$time_sec, row$plane, slice_num))

    # Transform y-coordinate (invert y-axis)
    gaze_y_plot <- row$gaze_y

    # Add current gaze point
    points(row$gaze_x, gaze_y_plot,
           col = "red", pch = 19, cex = 2)

    # Add crosshairs
    abline(v = row$gaze_x, col = "red", lty = 2, lwd = 0.5)
    abline(h = gaze_y_plot, col = "red", lty = 2, lwd = 0.5)

    dev.off()
  }

  message(sprintf("Animation frames saved to %s/", output_dir))
  message("To create video: ffmpeg -framerate 10 -pattern_type glob -i '*.png' -c:v libx264 gaze_video.mp4")
}

#' Map gaze coordinates to anatomical locations
#'
#' @description
#' Converts gaze coordinates to anatomical coordinates in the NIfTI space.
#' Uses pre-adjusted coordinates from integrate_all_gaze_points() if available.
#'
#' @param integrated_data Data frame from integrate_all_gaze_points()
#' @param nifti_data List from preload_nifti_data()
#' @param plane Imaging plane (AXIAL, SAGITTAL, or CORONAL)
#'
#' @return Data frame with anatomical coordinates for each gaze point
#' @export
#' @importFrom dplyr mutate select filter case_when
map_gaze_to_anatomy <- function(integrated_data, nifti_data, plane = NULL) {

  nvimg <- nifti_data$nvimage
  dims <- nifti_data$dims

  # Filter by plane if specified
  if (!is.null(plane)) {
    integrated_data <- integrated_data %>%
      filter(plane == !!plane)
  }

  # Use adjusted coordinates if available
  if ("gaze_x_adjusted" %in% names(integrated_data)) {
    # Filter to only points within bounds
    if ("within_bounds" %in% names(integrated_data)) {
      integrated_data <- integrated_data %>%
        filter(within_bounds)
    }
    coord_x <- "gaze_x_adjusted"
    coord_y <- "gaze_y_adjusted"
  } else {
    coord_x <- "gaze_x"
    coord_y <- "gaze_y"
  }

  # Map gaze coordinates to anatomical space
  anatomical_coords <- integrated_data %>%
    mutate(
      # Invert y-axis for image coordinate system
      gaze_y_transformed = !!sym(coord_y),
      gaze_x_transformed = !!sym(coord_x),

      # Convert normalized slice index to actual slice number
      slice_num = case_when(
        plane == "AXIAL" ~ round(slice_index * (dims[3] - 1)) + 1,
        plane == "SAGITTAL" ~ round(slice_index * (dims[1] - 1)) + 1,
        plane == "CORONAL" ~ round(slice_index * (dims[2] - 1)) + 1,
        TRUE ~ NA_real_
      ),

      # Convert gaze coordinates to voxel coordinates
      voxel_x = case_when(
        plane == "AXIAL" ~ round(gaze_x_transformed * (dims[1] - 1)),
        plane == "SAGITTAL" ~ round(gaze_x_transformed * (dims[2] - 1)),
        plane == "CORONAL" ~ round(gaze_x_transformed * (dims[1] - 1)),
        TRUE ~ NA_real_
      ),

      voxel_y = case_when(
        plane == "AXIAL" ~ round(gaze_y_transformed * (dims[2] - 1)),
        plane == "SAGITTAL" ~ round(gaze_y_transformed * (dims[3] - 1)),
        plane == "CORONAL" ~ round(gaze_y_transformed * (dims[3] - 1)),
        TRUE ~ NA_real_
      ),

      voxel_z = case_when(
        plane == "AXIAL" ~ slice_num - 1,
        plane == "SAGITTAL" ~ slice_num - 1,
        plane == "CORONAL" ~ slice_num - 1,
        TRUE ~ NA_real_
      )
    )

  # Convert voxel coordinates to world coordinates
  anatomical_coords <- anatomical_coords %>%
    mutate(
      # Create fractional coordinates
      frac_x = voxel_x / (dims[1] - 1),
      frac_y = voxel_y / (dims[2] - 1),
      frac_z = voxel_z / (dims[3] - 1)
    )

  # Apply coordinate transformation for each point
  world_coords <- matrix(NA, nrow = nrow(anatomical_coords), ncol = 3)

  for (i in 1:nrow(anatomical_coords)) {
    if (!is.na(anatomical_coords$frac_x[i])) {
      frac <- c(anatomical_coords$frac_x[i],
                anatomical_coords$frac_y[i],
                anatomical_coords$frac_z[i])

      # Convert to world coordinates
      world <- nvimg$convertFrac2MM(frac, isForceSliceMM = TRUE)
      world_coords[i, ] <- world[1:3]
    }
  }

  # Add world coordinates to the data frame
  anatomical_coords$world_x <- world_coords[, 1]
  anatomical_coords$world_y <- world_coords[, 2]
  anatomical_coords$world_z <- world_coords[, 3]

  # Select relevant columns
  result <- anatomical_coords %>%
    select(
      time_sec, time_aligned,
      plane, slice_num,
      gaze_x, gaze_y,
      gaze_x_transformed, gaze_y_transformed,
      voxel_x, voxel_y, voxel_z,
      world_x, world_y, world_z,
      everything()
    )

  return(result)
}

#' Get anatomical region for gaze points
#'
#' @description
#' Determines the anatomical region (if available from atlas) for each gaze point
#' based on its world coordinates.
#'
#' @param anatomical_coords Data frame from map_gaze_to_anatomy()
#' @param atlas_data Optional atlas NIfTI data for region lookup
#'
#' @return Data frame with anatomical region labels
#' @export
get_anatomical_regions <- function(anatomical_coords, atlas_data = NULL) {

  if (is.null(atlas_data)) {
    message("No atlas data provided. Returning coordinates without region labels.")
    anatomical_coords$region <- NA
    return(anatomical_coords)
  }

  # For each gaze point, look up the region in the atlas
  regions <- numeric(nrow(anatomical_coords))

  for (i in 1:nrow(anatomical_coords)) {
    if (!is.na(anatomical_coords$voxel_x[i])) {
      # Get voxel indices (add 1 for R's 1-based indexing)
      x <- anatomical_coords$voxel_x[i] + 1
      y <- anatomical_coords$voxel_y[i] + 1
      z <- anatomical_coords$voxel_z[i] + 1

      # Bounds check
      if (x >= 1 && x <= dim(atlas_data$data)[1] &&
          y >= 1 && y <= dim(atlas_data$data)[2] &&
          z >= 1 && z <= dim(atlas_data$data)[3]) {

        regions[i] <- atlas_data$data[x, y, z]
      } else {
        regions[i] <- NA
      }
    } else {
      regions[i] <- NA
    }
  }

  anatomical_coords$region_id <- regions

  # If you have a lookup table for region names, apply it here
  # anatomical_coords$region_name <- lookup_region_name(regions)

  return(anatomical_coords)
}

#' Create summary visualization of gaze patterns
#'
#' @param nifti_data List from preload_nifti_data()
#' @param integrated_data Data frame from integrate_all_gaze_points()
#' @export
#' @importFrom graphics plot barplot contour
#' @importFrom ggplot2 ggplot aes geom_point geom_bar labs theme_minimal coord_flip
#' @importFrom dplyr group_by summarise count top_n
create_gaze_summary_plot <- function(nifti_data, integrated_data) {

  # Set up 2x2 plot
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

  # 1. Timeline of slice viewing
  integrated_data %>%
    mutate(slice_label = paste(plane, round(slice_index * 100))) %>%
    ggplot(aes(x = time_sec, y = slice_label, color = plane)) +
    geom_point() +
    labs(title = "Slice Viewing Timeline", x = "Time (s)", y = "Slice") +
    theme_minimal() -> p1
  plot(p1)

  # 2. Gaze density across all slices
  # Use adjusted coordinates if available
  if ("gaze_x_adjusted" %in% names(integrated_data)) {
    valid_data <- integrated_data
    if ("within_bounds" %in% names(integrated_data)) {
      valid_data <- valid_data %>% filter(within_bounds)
    }
    gaze_x_plot <- valid_data$gaze_x_adjusted
    gaze_y_plot <- valid_data$gaze_y_adjusted
  } else {
    gaze_x_plot <- integrated_data$gaze_x
    gaze_y_plot <- integrated_data$gaze_y
  }

  plot(gaze_x_plot, gaze_y_plot,
       pch = 19, cex = 0.3, col = rgb(0, 0, 0, 0.1),
       main = "Overall Gaze Distribution",
       xlab = "Gaze X", ylab = "Gaze Y (transformed)",
       xlim = c(0, 1), ylim = c(0, 1))

  # Add density contours
  if (length(gaze_x_plot) > 50) {
    kde <- MASS::kde2d(gaze_x_plot, gaze_y_plot, n = 50)
    contour(kde, add = TRUE, col = "blue", lwd = 2)
  }

  # 3. Viewing duration by plane
  integrated_data %>%
    group_by(plane) %>%
    summarise(duration = max(time_sec) - min(time_sec)) %>%
    barplot(duration ~ plane, data = .,
            main = "Viewing Duration by Plane",
            ylab = "Duration (seconds)",
            col = c("lightblue", "lightgreen", "lightcoral"))

  # 4. Number of gazes per slice
  integrated_data %>%
    mutate(slice_id = paste(plane, round(slice_index * 100))) %>%
    count(slice_id) %>%
    top_n(10, n) %>%
    ggplot(aes(x = reorder(slice_id, n), y = n)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = "Top 10 Most Viewed Slices", x = "", y = "Number of Gazes") +
    theme_minimal() -> p4
  plot(p4)

  par(mfrow = c(1, 1))
}
