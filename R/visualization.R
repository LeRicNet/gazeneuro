#' Plot NIfTI slice with synchronized gaze data - IMAGE TRANSFORMATION
#'
#' @description
#' Displays a NIfTI slice overlaid with gaze tracking data.
#' The image is transformed to match where it was displayed during eye tracking.
#'
#' @param nifti_data List from preload_nifti_data()
#' @param integrated_data Data frame from integrate_all_gaze_points()
#' @param slice_num Slice number to display
#' @param plane Imaging plane (AXIAL, SAGITTAL, or CORONAL)
#' @param show_density Whether to show density contours for many points
#' @param image_transform List with x_min, x_max defining where image was displayed
#' @param auto_detect_transform Auto-detect image bounds from gaze distribution
#'
#' @return Data frame of gaze points for the displayed slice
#' @export
plot_slice_with_all_gaze <- function(nifti_data, integrated_data,
                                     slice_num = 12, plane = "AXIAL",
                                     show_density = TRUE,
                                     image_transform = NULL,
                                     auto_detect_transform = TRUE) {

  # Get all gaze points for this slice
  gaze_points <- get_all_gaze_for_slice(integrated_data, slice_num, plane,
                                        dims = dim(nifti_data$data))

  gaze_points$gaze_y <- 1 - gaze_points$gaze_y

  message(sprintf("\nSlice %d: %d gaze points over %.1f seconds",
                  slice_num, nrow(gaze_points),
                  ifelse(nrow(gaze_points) > 0,
                         max(gaze_points$time_aligned) - min(gaze_points$time_aligned),
                         0)))

  # AUTO-DETECT IMAGE BOUNDS from gaze distribution
  if (is.null(image_transform) && auto_detect_transform && nrow(gaze_points) > 10) {
    # Find the actual bounds where gaze points cluster
    # This assumes most gaze is on the image, not the background

    # Use density or percentiles to find image bounds
    x_density <- density(gaze_points$gaze_x, na.rm = TRUE)
    y_density <- density(gaze_points$gaze_y, na.rm = TRUE)

    # Simple approach: use 2nd and 98th percentiles
    x_bounds <- quantile(gaze_points$gaze_x, c(0.02, 0.98), na.rm = TRUE)
    y_bounds <- quantile(gaze_points$gaze_y, c(0.02, 0.98), na.rm = TRUE)

    # Add small padding
    x_range <- x_bounds[2] - x_bounds[1]
    y_range <- y_bounds[2] - y_bounds[1]

    image_transform <- list(
      x_min = x_bounds[1] - 0.02 * x_range,
      x_max = x_bounds[2] + 0.02 * x_range,
      y_min = y_bounds[1] - 0.02 * y_range,
      y_max = y_bounds[2] + 0.02 * y_range
    )

    message("\nAuto-detected image bounds from gaze distribution:")
    message(sprintf("  X: [%.3f, %.3f]", image_transform$x_min, image_transform$x_max))
    message(sprintf("  Y: [%.3f, %.3f]", image_transform$y_min, image_transform$y_max))
  }

  # Default to full canvas if no transform specified
  if (is.null(image_transform)) {
    image_transform <- list(x_min = 0, x_max = 1, y_min = 0, y_max = 1)
  }

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

  # Set up the plot area to show full gaze coordinate space [0,1]
  plot(NULL, NULL,
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Gaze X", ylab = "Gaze Y",
       main = sprintf("%s Slice %d - %d gaze points",
                      plane, slice_num, nrow(gaze_points)))

  # Draw the transformed NIfTI image at the detected/specified location
  # Use rasterImage to place the image at specific coordinates
  if (!is.null(slice_data)) {

    slice_data <- apply(t(slice_data), 2, rev)

    # Convert slice data to raster
    slice_raster <- as.raster(slice_data / max(slice_data, na.rm = TRUE))

    # Draw image at transformed coordinates
    rasterImage(slice_raster,
                xleft = image_transform$x_min,
                ybottom = image_transform$y_min,
                xright = image_transform$x_max,
                ytop = image_transform$y_max)

    # Draw image boundary
    rect(image_transform$x_min, image_transform$y_min,
         image_transform$x_max, image_transform$y_max,
         border = "blue", lwd = 2)
  }

  # Now overlay gaze points (unchanged)
  if (nrow(gaze_points) > 0) {
    if (show_density && nrow(gaze_points) > 50) {
      # Create density overlay
      kde <- MASS::kde2d(gaze_points$gaze_x, gaze_points$gaze_y,
                         n = 50, lims = c(0, 1, 0, 1))
      contour(kde, add = TRUE, col = heat.colors(10), lwd = 2)
    } else {
      # Plot individual points
      points(gaze_points$gaze_x, gaze_points$gaze_y,
             col = rgb(1, 0, 0, 0.3), pch = 19, cex = 0.5)
    }

    # Show temporal progression with color
    n_points <- nrow(gaze_points)
    if (n_points > 1 && n_points < 500) {
      colors <- colorRampPalette(c("green", "yellow", "red"))(n_points)
      points(gaze_points$gaze_x, gaze_points$gaze_y,
             col = colors, pch = 19, cex = 0.7)
    }

    # Add mean position
    mean_x <- mean(gaze_points$gaze_x)
    mean_y <- mean(gaze_points$gaze_y)
    points(mean_x, mean_y, col = "blue", pch = 3, cex = 2, lwd = 3)

    # Add stats
    text(0.05, 0.95, sprintf("n = %d", nrow(gaze_points)),
         adj = c(0, 1), cex = 0.8)
    text(0.05, 0.90, sprintf("%.1f sec",
                             max(gaze_points$time_aligned) - min(gaze_points$time_aligned)),
         adj = c(0, 1), cex = 0.8)

    # Show transform info
    if (!all(c(image_transform$x_min, image_transform$x_max) == c(0, 1))) {
      text(0.95, 0.05, "Image transformed",
           adj = c(1, 0), cex = 0.8, col = "darkgreen")
    }
  }

  return(invisible(gaze_points))
}

#' Detect image display bounds from gaze data
#'
#' @description
#' Analyzes gaze distribution across all slices to determine where
#' the NIfTI image was displayed on screen during eye tracking.
#'
#' @param integrated_data Data from integrate_all_gaze_points()
#' @param plane Imaging plane to analyze
#' @param method Detection method: "percentile", "density", or "convex_hull"
#' @return List with x_min, x_max, y_min, y_max
#' @export
detect_image_bounds <- function(integrated_data, plane = "SAGITTAL",
                                method = "percentile") {

  # Get all gaze points for this plane
  plane_data <- integrated_data %>%
    filter(plane == !!plane)

  if (nrow(plane_data) < 100) {
    warning("Insufficient gaze data for reliable bound detection")
    return(list(x_min = 0, x_max = 1, y_min = 0, y_max = 1))
  }

  bounds <- switch(method,
                   "percentile" = {
                     # Use percentiles to exclude outliers
                     x_bounds <- quantile(plane_data$gaze_x, c(0.01, 0.99), na.rm = TRUE)
                     y_bounds <- quantile(plane_data$gaze_y, c(0.01, 0.99), na.rm = TRUE)

                     # Add small padding
                     x_pad <- 0.02 * (x_bounds[2] - x_bounds[1])
                     y_pad <- 0.02 * (y_bounds[2] - y_bounds[1])

                     list(
                       x_min = x_bounds[1] - x_pad,
                       x_max = x_bounds[2] + x_pad,
                       y_min = y_bounds[1] - y_pad,
                       y_max = y_bounds[2] + y_pad
                     )
                   },

                   "density" = {
                     # Find high-density region
                     kde <- MASS::kde2d(plane_data$gaze_x, plane_data$gaze_y, n = 100)

                     # Find contour containing 95% of points
                     # This is more complex - simplified version:
                     x_density <- density(plane_data$gaze_x)
                     y_density <- density(plane_data$gaze_y)

                     # Find peaks and estimate bounds
                     x_peak <- x_density$x[which.max(x_density$y)]
                     y_peak <- y_density$x[which.max(y_density$y)]

                     # Estimate spread
                     x_spread <- sd(plane_data$gaze_x) * 3
                     y_spread <- sd(plane_data$gaze_y) * 3

                     list(
                       x_min = max(0, x_peak - x_spread),
                       x_max = min(1, x_peak + x_spread),
                       y_min = max(0, y_peak - y_spread),
                       y_max = min(1, y_peak + y_spread)
                     )
                   },

                   "convex_hull" = {
                     # Use convex hull of gaze points
                     hull_indices <- chull(plane_data$gaze_x, plane_data$gaze_y)
                     hull_x <- plane_data$gaze_x[hull_indices]
                     hull_y <- plane_data$gaze_y[hull_indices]

                     list(
                       x_min = min(hull_x),
                       x_max = max(hull_x),
                       y_min = min(hull_y),
                       y_max = max(hull_y)
                     )
                   }
  )

  # Ensure bounds are within [0,1]
  bounds$x_min <- max(0, bounds$x_min)
  bounds$x_max <- min(1, bounds$x_max)
  bounds$y_min <- max(0, bounds$y_min)
  bounds$y_max <- min(1, bounds$y_max)

  # Report findings
  message(sprintf("\nDetected image bounds (%s method):", method))
  message(sprintf("  X: [%.3f, %.3f] (width: %.3f)",
                  bounds$x_min, bounds$x_max, bounds$x_max - bounds$x_min))
  message(sprintf("  Y: [%.3f, %.3f] (height: %.3f)",
                  bounds$y_min, bounds$y_max, bounds$y_max - bounds$y_min))

  return(bounds)
}

#' Visualize detected image bounds
#'
#' @param integrated_data Data from integrate_all_gaze_points()
#' @param image_bounds Bounds from detect_image_bounds()
#' @param plane Imaging plane
#' @export
visualize_image_bounds <- function(integrated_data, image_bounds = NULL,
                                   plane = "SAGITTAL") {

  plane_data <- integrated_data %>%
    filter(plane == !!plane)

  # Auto-detect if not provided
  if (is.null(image_bounds)) {
    image_bounds <- detect_image_bounds(integrated_data, plane)
  }

  # Need to invert the y-axis
  plane_data$gaze_y <- 1 - plane_data$gaze_y

  # Create visualization
  plot(plane_data$gaze_x, plane_data$gaze_y,
       pch = 19, cex = 0.2, col = rgb(0, 0, 0, 0.1),
       main = sprintf("Detected Image Display Area - %s Plane", plane),
       xlab = "Gaze X", ylab = "Gaze Y",
       xlim = c(0, 1), ylim = c(0, 1))

  # Show full canvas
  rect(0, 0, 1, 1, border = "gray", lwd = 1, lty = 2)
  text(0.5, -0.02, "Full canvas [0,1]", cex = 0.8, col = "gray", xpd = TRUE)

  # Show detected image area
  rect(image_bounds$x_min, image_bounds$y_min,
       image_bounds$x_max, image_bounds$y_max,
       border = "blue", lwd = 3)

  # Add labels
  text(image_bounds$x_min, image_bounds$y_min - 0.02,
       sprintf("(%.2f, %.2f)", image_bounds$x_min, image_bounds$y_min),
       cex = 0.7, col = "blue", adj = c(0, 1), xpd = TRUE)

  text(image_bounds$x_max, image_bounds$y_max + 0.02,
       sprintf("(%.2f, %.2f)", image_bounds$x_max, image_bounds$y_max),
       cex = 0.7, col = "blue", adj = c(1, 0), xpd = TRUE)

  # Show center
  cx <- mean(c(image_bounds$x_min, image_bounds$x_max))
  cy <- mean(c(image_bounds$y_min, image_bounds$y_max))
  points(cx, cy, pch = 3, cex = 2, col = "blue", lwd = 2)

  # Add density contour
  if (nrow(plane_data) > 50) {
    kde <- MASS::kde2d(plane_data$gaze_x, plane_data$gaze_y, n = 50)
    contour(kde, add = TRUE, col = "red", levels = 5)
  }

  legend("topright",
         legend = c("Gaze points", "Detected image area", "Density contours"),
         col = c("black", "blue", "red"),
         pch = c(19, NA, NA),
         lty = c(NA, 1, 1),
         lwd = c(NA, 3, 1),
         cex = 0.8)
}
#' Plot all viewed slices in a grid with gaze data
#'
#' @description
#' Creates a grid visualization showing all slices that were viewed,
#' each overlaid with its corresponding gaze data.
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
    slice_normalized <- (slice_num - 1) / (nifti_data$nvimage$dimsRAS[
      ifelse(plane == "AXIAL", 4, ifelse(plane == "SAGITTAL", 2, 3))
    ] - 1)

    gaze_points <- integrated_data %>%
      filter(
        plane == !!plane,
        abs(slice_index - slice_normalized) < 0.02
      )

    # Plot the slice
    image(slice_data,
          col = gray((0:255)/255),
          axes = FALSE,
          main = sprintf("Slice %d (n=%d)", slice_num, nrow(gaze_points)),
          cex.main = 0.8)

    # Add gaze points if any
    if (nrow(gaze_points) > 0) {
      # Add all gaze points
      points(gaze_points$gaze_x,
             gaze_points$gaze_y,
             col = rgb(1, 0, 0, 0.4),
             pch = 19,
             cex = 0.3)

      # Add trajectory lines if more than 5 points
      if (nrow(gaze_points) > 5) {
        lines(gaze_points$gaze_x,
              gaze_points$gaze_y,
              col = rgb(1, 0, 0, 0.2),
              lwd = 1)
      }

      # Mark center of gaze
      mean_x <- mean(gaze_points$gaze_x)
      mean_y <- mean(gaze_points$gaze_y)
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
  summary_stats <- integrated_data %>%
    filter(plane == !!plane) %>%
    mutate(
      slice_num = case_when(
        plane == "AXIAL" ~ round(slice_index * (nifti_data$nvimage$dimsRAS[4] - 1)) + 1,
        plane == "SAGITTAL" ~ round(slice_index * (nifti_data$nvimage$dimsRAS[2] - 1)) + 1,
        plane == "CORONAL" ~ round(slice_index * (nifti_data$nvimage$dimsRAS[3] - 1)) + 1
      )
    ) %>%
    group_by(slice_num) %>%
    summarise(
      n_gazes = n(),
      duration = max(time_sec) - min(time_sec),
      mean_x = mean(gaze_x),
      mean_y = mean(gaze_y),
      sd_x = sd(gaze_x),
      sd_y = sd(gaze_y),
      .groups = "drop"
    )

  return(invisible(summary_stats))
}

#' Create gaze animation frames
#'
#' @description
#' Generates PNG frames showing gaze movement over time on NIfTI slices.
#' Frames can be combined into a video using ffmpeg.
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

  # Get unique time points where slices changed
  slice_changes <- integrated_data %>%
    group_by(time_sec, plane, slice_index) %>%
    summarise(
      gaze_x = mean(gaze_x),
      gaze_y = mean(gaze_y),
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

    # Add current gaze point
    points(row$gaze_x, row$gaze_y,
           col = "red", pch = 19, cex = 2)

    # Add crosshairs
    abline(v = row$gaze_x, col = "red", lty = 2, lwd = 0.5)
    abline(h = row$gaze_y, col = "red", lty = 2, lwd = 0.5)

    dev.off()
  }

  message(sprintf("Animation frames saved to %s/", output_dir))
  message("To create video: ffmpeg -framerate 10 -pattern_type glob -i '*.png' -c:v libx264 gaze_video.mp4")
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
  plot(integrated_data$gaze_x, integrated_data$gaze_y,
       pch = 19, cex = 0.3, col = rgb(0, 0, 0, 0.1),
       main = "Overall Gaze Distribution",
       xlab = "Gaze X", ylab = "Gaze Y",
       xlim = c(0, 1), ylim = c(0, 1))

  # Add density contours
  if (nrow(integrated_data) > 50) {
    kde <- MASS::kde2d(integrated_data$gaze_x, integrated_data$gaze_y, n = 50)
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
