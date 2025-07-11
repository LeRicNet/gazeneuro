#' Integrate gaze tracking data with z-axis slice information
#'
#' @description
#' Aligns gaze tracking data temporally with z-axis events (slice changes)
#' and assigns each gaze point to the slice being viewed at that time.
#' Optionally transforms gaze coordinates based on image bounds.
#'
#' @param gaze_data Data frame containing gaze tracking data with columns:
#'   \describe{
#'     \item{device_time_stamp}{Timestamp in microseconds}
#'     \item{gaze_point_on_display_area_x}{X coordinate (0-1)}
#'     \item{gaze_point_on_display_area_y}{Y coordinate (0-1)}
#'   }
#' @param z_axis Data frame containing z-axis slice information with columns:
#'   \describe{
#'     \item{client_timestamp}{Timestamp in milliseconds}
#'     \item{plane}{Imaging plane (AXIAL, SAGITTAL, CORONAL)}
#'     \item{index}{Normalized slice index (0-1)}
#'     \item{image_id}{Image identifier}
#'   }
#' @param latency_correction Optional manual latency correction in seconds
#' @param image_bounds Optional list with image position in display area:
#'   \describe{
#'     \item{left}{Left edge of image (0-1)}
#'     \item{right}{Right edge of image (0-1)}
#'     \item{top}{Top edge of image (0-1)}
#'     \item{bottom}{Bottom edge of image (0-1)}
#'   }
#'
#' @return Data frame with integrated gaze and slice information, including
#'   adjusted gaze coordinates if image_bounds provided
#' @export
#' @importFrom dplyr mutate select filter arrange left_join group_by summarise
integrate_all_gaze_points <- function(gaze_data, z_axis,
                                      latency_correction = NULL,
                                      image_bounds = NULL) {

  # 1. Normalize timestamps to seconds from start
  gaze_clean <- gaze_data %>%
    mutate(
      time_sec = (device_time_stamp - min(device_time_stamp)) / 1e6,
      gaze_id = row_number()  # Keep track of original gaze points
    ) %>%
    select(gaze_id, time_sec, gaze_x = gaze_point_on_display_area_x,
           gaze_y = gaze_point_on_display_area_y, everything())

  z_clean <- z_axis %>%
    mutate(
      time_sec = (client_timestamp - min(client_timestamp)) / 1000,
      z_event_id = row_number()
    ) %>%
    select(z_event_id, time_sec, plane, slice_index = index, image_id, everything())

  # 2. Calculate latency if not provided
  if (is.null(latency_correction)) {
    gaze_duration <- max(gaze_clean$time_sec)
    z_duration <- max(z_clean$time_sec)
    latency <- z_duration - gaze_duration
  } else {
    latency <- latency_correction
  }

  message("Integration Summary:")
  message("────────────────────────────────────────")
  message("Gaze points: ", nrow(gaze_clean))
  message("Z-axis events: ", nrow(z_clean))
  message("Gaze duration: ", round(max(gaze_clean$time_sec), 3), " sec")
  message("Z-axis duration: ", round(max(z_clean$time_sec), 3), " sec")
  message("Estimated latency: ", round(latency, 3), " sec")

  if (!is.null(image_bounds)) {
    message("Image bounds provided: left=", round(image_bounds$left, 3),
            ", right=", round(image_bounds$right, 3),
            ", top=", round(image_bounds$top, 3),
            ", bottom=", round(image_bounds$bottom, 3))
  }
  message("")

  # Apply latency correction to gaze
  gaze_clean$time_aligned <- gaze_clean$time_sec + latency

  # 3. For each gaze point, find which slice was being viewed
  # Create intervals from z-axis events
  z_intervals <- z_clean %>%
    arrange(time_sec) %>%
    mutate(
      time_start = time_sec,
      time_end = lead(time_sec, default = max(time_sec) + 1)
    )

  # Method 1: Using a join approach (faster for large data)
  integrated <- gaze_clean %>%
    mutate(dummy = 1) %>%
    left_join(z_intervals %>% mutate(dummy = 1), by = "dummy") %>%
    filter(
      time_aligned >= time_start,
      time_aligned < time_end
    ) %>%
    select(-dummy, -time_start, -time_end)

  # Alternative Method 2: Using cut() (if Method 1 has issues)
  if (nrow(integrated) == 0) {
    message("Using alternative integration method...")

    # Create breaks from z-axis events
    breaks <- c(-Inf, z_intervals$time_sec, Inf)

    # Assign each gaze point to an interval
    gaze_clean$interval <- cut(gaze_clean$time_aligned,
                               breaks = breaks,
                               labels = FALSE,
                               right = FALSE)

    # Merge with z-axis data
    integrated <- gaze_clean %>%
      mutate(z_event_id = pmin(interval, nrow(z_intervals))) %>%
      left_join(z_intervals, by = "z_event_id")
  }

  # 4. Apply image bounds transformation if provided
  if (!is.null(image_bounds)) {
    # Calculate image dimensions within display
    image_width <- image_bounds$right - image_bounds$left
    image_height <- image_bounds$bottom - image_bounds$top

    integrated <- integrated %>%
      mutate(
        # Transform display coordinates to image space (0-1)
        gaze_x_adjusted = (gaze_x - image_bounds$left) / image_width,
        gaze_y_adjusted = (gaze_y - image_bounds$top) / image_height,

        # Flag points that are within image bounds
        within_bounds = gaze_x_adjusted >= 0 & gaze_x_adjusted <= 1 &
          gaze_y_adjusted >= 0 & gaze_y_adjusted <= 1
      )

    # Report how many points are outside bounds
    n_outside <- sum(!integrated$within_bounds)
    if (n_outside > 0) {
      message(sprintf("Note: %d gaze points (%.1f%%) fall outside image bounds",
                      n_outside, 100 * n_outside / nrow(integrated)))
    }
  } else {
    # No image bounds - use original coordinates
    integrated <- integrated %>%
      mutate(
        gaze_x_adjusted = gaze_x,
        gaze_y_adjusted = gaze_y,
        within_bounds = TRUE
      )
  }

  message("Integration Results:")
  message("────────────────────────────────────────")
  message("Gaze points matched: ", nrow(integrated),
          sprintf(" (%.1f%%)", 100 * nrow(integrated) / nrow(gaze_clean)))

  if (!is.null(image_bounds)) {
    n_within <- sum(integrated$within_bounds)
    message("Gaze points within image bounds: ", n_within,
            sprintf(" (%.1f%%)", 100 * n_within / nrow(integrated)))
  }

  # Check for unmatched points
  unmatched <- nrow(gaze_clean) - nrow(integrated)
  if (unmatched > 0) {
    message("Unmatched gaze points: ", unmatched)
    message("(These may be before first or after last z-axis event)")
  }

  return(integrated)
}

#' Get all gaze points for a specific slice
#'
#' @param integrated_data Data frame from integrate_all_gaze_points()
#' @param slice_num Slice number (1-based)
#' @param plane Imaging plane (AXIAL, SAGITTAL, or CORONAL)
#' @param dims Dimensions of the image (x, y, z)
#' @param filter_bounds Whether to filter out points outside image bounds
#'
#' @return Data frame of gaze points for the specified slice
#' @export
get_all_gaze_for_slice <- function(integrated_data, slice_num, plane = "AXIAL",
                                   dims = c(512, 512, 25), filter_bounds = TRUE) {

  # Convert slice number to normalized index
  if (plane == "AXIAL") {
    slice_normalized <- (slice_num - 1) / (dims[3] - 1)
  } else if (plane == "SAGITTAL") {
    slice_normalized <- (slice_num - 1) / (dims[1] - 1)
  } else {
    slice_normalized <- (slice_num - 1) / (dims[2] - 1)
  }

  # Get all gaze points for this slice
  gaze_points <- integrated_data %>%
    filter(
      plane == !!plane,
      abs(slice_index - slice_normalized) < 0.02
    )

  # Filter out points outside bounds if requested
  if (filter_bounds && "within_bounds" %in% names(gaze_points)) {
    gaze_points <- gaze_points %>%
      filter(within_bounds)
  }

  gaze_points %>%
    arrange(time_aligned)
}

#' Check integration quality and visualize alignment
#'
#' @description
#' Creates diagnostic plots to verify the quality of gaze-slice integration.
#'
#' @param gaze_data Original gaze data
#' @param z_axis Original z-axis data
#' @param integrated Integrated data from integrate_all_gaze_points()
#' @export
#' @importFrom graphics plot points axis legend barplot rect
#' @importFrom dplyr count arrange group_by summarise
check_integration_quality <- function(gaze_data, z_axis, integrated) {

  # Set up plots
  par(mfrow = c(3, 1), mar = c(4, 4, 3, 2))

  # 1. Original gaze timeline
  gaze_times <- (gaze_data$device_time_stamp - min(gaze_data$device_time_stamp)) / 1e6
  plot(gaze_times, rep(1, length(gaze_times)),
       pch = ".", col = "blue", ylim = c(0.5, 3.5),
       main = "Data Timeline Comparison",
       xlab = "Time (seconds)", ylab = "",
       yaxt = "n")

  # Add z-axis events
  z_times <- (z_axis$client_timestamp - min(z_axis$client_timestamp)) / 1000
  points(z_times, rep(2, length(z_times)), pch = 19, col = "red", cex = 1)

  # Add integrated (aligned) gaze
  if (nrow(integrated) > 0) {
    points(integrated$time_aligned, rep(3, nrow(integrated)),
           pch = ".", col = "green")
  }

  axis(2, at = 1:3, labels = c("Original Gaze", "Z-axis Events", "Aligned Gaze"))
  legend("topright",
         legend = c("Gaze (original)", "Z-axis events", "Gaze (aligned)"),
         col = c("blue", "red", "green"),
         pch = c(46, 19, 46))

  # 2. Gaze count per z-axis event
  event_counts <- integrated %>%
    count(z_event_id) %>%
    arrange(z_event_id)

  barplot(event_counts$n,
          names.arg = event_counts$z_event_id,
          main = "Gaze Points per Z-axis Event",
          xlab = "Z-axis Event ID", ylab = "Number of Gaze Points",
          col = "lightblue")

  # 3. Gaze distribution with image bounds
  if ("gaze_x_adjusted" %in% names(integrated)) {
    plot(integrated$gaze_x, integrated$gaze_y,
         pch = 19, cex = 0.3, col = rgb(0, 0, 0, 0.3),
         xlim = c(0, 1), ylim = c(0, 1),
         main = "Gaze Distribution (Display Space)",
         xlab = "Display X", ylab = "Display Y")

    # Show image bounds if they exist
    if ("within_bounds" %in% names(integrated)) {
      # Infer bounds from adjusted coordinates
      valid_points <- integrated %>% filter(within_bounds)
      if (nrow(valid_points) > 0) {
        # This is approximate - better to pass bounds directly
        x_range <- range(integrated$gaze_x[integrated$within_bounds])
        y_range <- range(integrated$gaze_y[integrated$within_bounds])
        rect(x_range[1], y_range[1], x_range[2], y_range[2],
             border = "red", lwd = 2)
        legend("topright", legend = "Image bounds", col = "red", lty = 1, lwd = 2)
      }
    }
  } else {
    # Time coverage by plane
    coverage <- integrated %>%
      group_by(plane) %>%
      summarise(
        n_gazes = n(),
        time_span = max(time_aligned) - min(time_aligned),
        .groups = "drop"
      )

    barplot(coverage$n_gazes,
            names.arg = coverage$plane,
            main = "Gaze Points by Plane",
            ylab = "Number of Gaze Points",
            col = c("lightcoral", "lightgreen", "lightblue"))
  }

  par(mfrow = c(1, 1))

  # Print summary
  message("\n\nDetailed Summary:")
  message("────────────────────────────────────────")

  if ("within_bounds" %in% names(integrated)) {
    bounds_summary <- integrated %>%
      group_by(plane) %>%
      summarise(
        n_total = n(),
        n_within_bounds = sum(within_bounds),
        pct_within = round(100 * mean(within_bounds), 1),
        .groups = "drop"
      )
    print(bounds_summary)
  } else {
    coverage <- integrated %>%
      group_by(plane) %>%
      summarise(
        n_gazes = n(),
        time_span = max(time_aligned) - min(time_aligned),
        .groups = "drop"
      )
    print(coverage)
  }

  # Check for gaps
  gaps <- integrated %>%
    arrange(time_aligned) %>%
    mutate(time_diff = c(0, diff(time_aligned))) %>%
    filter(time_diff > 0.5)  # Gaps > 500ms

  if (nrow(gaps) > 0) {
    message("\nLarge time gaps found:")
    print(gaps %>% select(time_aligned, time_diff))
  }
}

#' Calculate gaze statistics by slice
#'
#' @param integrated_data Data frame from integrate_all_gaze_points()
#' @param nifti_data Optional NIfTI data for dimension information
#' @param use_adjusted Whether to use adjusted coordinates (if available)
#' @return Data frame with summary statistics per slice
#' @export
#' @importFrom dplyr mutate group_by summarise arrange case_when
analyze_gaze_by_slice <- function(integrated_data, nifti_data = NULL,
                                  use_adjusted = TRUE) {

  # Determine which coordinates to use
  if (use_adjusted && "gaze_x_adjusted" %in% names(integrated_data)) {
    coord_x <- "gaze_x_adjusted"
    coord_y <- "gaze_y_adjusted"

    # Filter to only within bounds if available
    if ("within_bounds" %in% names(integrated_data)) {
      integrated_data <- integrated_data %>%
        filter(within_bounds)
    }
  } else {
    coord_x <- "gaze_x"
    coord_y <- "gaze_y"
  }

  # Group by plane and slice
  slice_summary <- integrated_data %>%
    mutate(
      # Convert normalized index to slice number
      slice_num = case_when(
        plane == "AXIAL" ~ round(slice_index * 24) + 1,
        plane == "SAGITTAL" ~ round(slice_index * 511) + 1,
        plane == "CORONAL" ~ round(slice_index * 511) + 1
      )
    ) %>%
    group_by(plane, slice_num) %>%
    summarise(
      n_gazes = n(),
      duration_sec = max(time_sec) - min(time_sec),
      mean_x = mean(!!sym(coord_x)),
      mean_y = mean(!!sym(coord_y)),
      sd_x = sd(!!sym(coord_x)),
      sd_y = sd(!!sym(coord_y)),
      .groups = "drop"
    ) %>%
    arrange(plane, slice_num)

  return(slice_summary)
}
