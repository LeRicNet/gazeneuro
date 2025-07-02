# Function to plot slices showing valid (in-canvas) vs invalid (out-of-bounds) gaze points

plot_slice_with_validity <- function(nifti_data, integrated, slice_num,
                                     plane = "AXIAL", show_frame = TRUE) {

  # Get slice data
  if (plane == "AXIAL") {
    slice_data <- nifti_data$data[,, slice_num]
    slice_normalized <- (slice_num - 1) / (nifti_data$dims[3] - 1)
  } else if (plane == "SAGITTAL") {
    slice_data <- nifti_data$data[slice_num, , ]
    slice_normalized <- (slice_num - 1) / (nifti_data$dims[1] - 1)
  } else {
    slice_data <- nifti_data$data[, slice_num, ]
    slice_normalized <- (slice_num - 1) / (nifti_data$dims[2] - 1)
  }

  # Get gaze points for this slice
  slice_gazes <- integrated %>%
    filter(
      plane == !!plane,
      abs(slice_index - slice_normalized) < 0.02
    )

  if (nrow(slice_gazes) == 0) {
    message("No gaze points found for ", plane, " slice ", slice_num)
    return(invisible(NULL))
  }

  # Process coordinates to check validity
  validity_check <- map2_dfr(slice_gazes$gaze_x, slice_gazes$gaze_y, function(x, y) {
    adj <- resolve_coordinates(x, y)
    data.frame(
      x_original = x,
      y_original = y,
      x_adjusted = adj$x,
      y_adjusted = adj$y,
      in_bounds = adj$in_bounds
    )
  })

  # Create the plot with extra margin for legend
  par(mar = c(5, 4, 4, 8))  # Increased right margin for legend

  # Plot the brain slice
  image(slice_data, col = gray((0:255)/255),
        main = sprintf("%s Slice %d - Gaze Validity Check", plane, slice_num),
        xlab = "X", ylab = "Y", xlim = c(0, 1), ylim = c(0, 1))

  # Show frame boundaries if requested
  if (show_frame) {
    # Canvas boundaries in normalized coordinates (canvas aligned to bottom)
    left_edge <- 350 / 2624   # ~0.133
    right_edge <- (350 + 1924) / 2624  # ~0.867
    # With canvas at bottom: top edge is at 80px from top
    top_edge <- 80 / 1640     # ~0.049
    bottom_edge <- 1  # Canvas extends to bottom

    # Draw frame outline (in original Tobii coordinates)
    rect(0, 0, 1, 1, border = "darkgray", lty = 3, lwd = 1)

    # Draw canvas outline
    rect(left_edge, top_edge, right_edge, bottom_edge,
         border = "blue", lty = 2, lwd = 2)

    # Add small labels outside the image area
    text(0.5, -0.08, "Canvas boundaries", cex = 0.7, col = "blue", xpd = TRUE)
  }

  # Plot out-of-bounds points (original positions)
  out_bounds <- validity_check[!validity_check$in_bounds, ]
  if (nrow(out_bounds) > 0) {
    points(out_bounds$x_original, out_bounds$y_original,
           pch = 4, col = rgb(1, 0, 0, 0.7), cex = 1, lwd = 2)  # X marks
  }

  # Plot in-bounds points (adjusted positions)
  in_bounds <- validity_check[validity_check$in_bounds, ]
  if (nrow(in_bounds) > 0) {
    points(in_bounds$x_adjusted, in_bounds$y_adjusted,
           pch = 19, col = rgb(0, 0.8, 0, 0.5), cex = 0.6)  # Green dots

    # Show mean of valid points
    mean_x <- mean(in_bounds$x_adjusted)
    mean_y <- mean(in_bounds$y_adjusted)
    points(mean_x, mean_y, pch = 3, col = "blue", cex = 2, lwd = 2)
  }

  # Add legend outside plot area
  par(xpd = TRUE)  # Allow drawing outside plot region
  legend(1.02, 0.8,  # Position to the right of plot
         legend = c(
           sprintf("Valid: %d", nrow(in_bounds)),
           sprintf("Out: %d", nrow(out_bounds)),
           "Mean"
         ),
         pch = c(19, 4, 3),
         col = c(rgb(0, 0.8, 0, 0.5), rgb(1, 0, 0, 0.7), "blue"),
         cex = 0.7,
         bty = "n")  # No box around legend
  par(xpd = FALSE)

  # Print summary
  cat(sprintf("\nSlice %d summary:\n", slice_num))
  cat(sprintf("  Total gaze points: %d\n", nrow(slice_gazes)))
  cat(sprintf("  Valid (in canvas): %d (%.1f%%)\n",
              nrow(in_bounds), 100 * nrow(in_bounds) / nrow(slice_gazes)))
  cat(sprintf("  Out of bounds: %d (%.1f%%)\n",
              nrow(out_bounds), 100 * nrow(out_bounds) / nrow(slice_gazes)))

  return(invisible(validity_check))
}

# Alternative version with legend at bottom
plot_slice_with_validity_bottom_legend <- function(nifti_data, integrated, slice_num,
                                                   plane = "AXIAL", show_frame = TRUE) {

  # Get slice data
  if (plane == "AXIAL") {
    slice_data <- nifti_data$data[,, slice_num]
    slice_normalized <- (slice_num - 1) / (nifti_data$dims[3] - 1)
  } else if (plane == "SAGITTAL") {
    slice_data <- nifti_data$data[slice_num, , ]
    slice_normalized <- (slice_num - 1) / (nifti_data$dims[1] - 1)
  } else {
    slice_data <- nifti_data$data[, slice_num, ]
    slice_normalized <- (slice_num - 1) / (nifti_data$dims[2] - 1)
  }

  # Get gaze points for this slice
  slice_gazes <- integrated %>%
    filter(
      plane == !!plane,
      abs(slice_index - slice_normalized) < 0.02
    )

  if (nrow(slice_gazes) == 0) {
    message("No gaze points found for ", plane, " slice ", slice_num)
    return(invisible(NULL))
  }

  # Process coordinates to check validity
  validity_check <- map2_dfr(slice_gazes$gaze_x, slice_gazes$gaze_y, function(x, y) {
    adj <- resolve_coordinates(x, y)
    data.frame(
      x_original = x,
      y_original = y,
      x_adjusted = adj$x,
      y_adjusted = adj$y,
      in_bounds = adj$in_bounds
    )
  })

  # Count points
  n_valid <- sum(validity_check$in_bounds)
  n_invalid <- sum(!validity_check$in_bounds)

  # Create the plot
  par(mar = c(6, 4, 4, 2))  # Extra bottom margin

  # Plot the brain slice
  image(slice_data, col = gray((0:255)/255),
        main = sprintf("%s Slice %d - Valid: %d (%.0f%%), Invalid: %d",
                       plane, slice_num, n_valid,
                       100 * n_valid / nrow(slice_gazes), n_invalid),
        xlab = "X", ylab = "Y", xlim = c(0, 1), ylim = c(0, 1))

  # Show frame boundaries if requested
  if (show_frame) {
    # Canvas boundaries
    rect(350/2624, 80/1640, (350+1924)/2624, 1,
         border = rgb(0, 0, 1, 0.5), lty = 2, lwd = 1)
  }

  # Plot points
  if (n_invalid > 0) {
    out_bounds <- validity_check[!validity_check$in_bounds, ]
    points(out_bounds$x_original, out_bounds$y_original,
           pch = 4, col = rgb(1, 0, 0, 0.6), cex = 0.8, lwd = 1)
  }

  if (n_valid > 0) {
    in_bounds <- validity_check[validity_check$in_bounds, ]
    points(in_bounds$x_adjusted, in_bounds$y_adjusted,
           pch = 19, col = rgb(0, 0.8, 0, 0.4), cex = 0.5)

    # Mean position
    points(mean(in_bounds$x_adjusted), mean(in_bounds$y_adjusted),
           pch = 3, col = "blue", cex = 1.5, lwd = 2)
  }

  # Simple text summary at bottom
  mtext(sprintf("Green dots = valid gazes, Red X = out of bounds, Blue + = mean position"),
        side = 1, line = 4.5, cex = 0.7)

  return(invisible(validity_check))
}

# Batch function to check all slices
check_all_slice_validity <- function(nifti_data, integrated, plane = "AXIAL") {

  # Get unique slices
  dims <- nifti_data$dims
  slice_nums <- integrated %>%
    filter(plane == !!plane) %>%
    mutate(
      slice_num = case_when(
        plane == "AXIAL" ~ round(slice_index * (dims[3] - 1)) + 1,
        plane == "SAGITTAL" ~ round(slice_index * (dims[1] - 1)) + 1,
        plane == "CORONAL" ~ round(slice_index * (dims[2] - 1)) + 1
      )
    ) %>%
    distinct(slice_num) %>%
    arrange(slice_num) %>%
    pull(slice_num)

  # Check each slice
  validity_summary <- map_dfr(slice_nums, function(s) {
    slice_gazes <- integrated %>%
      filter(
        plane == !!plane,
        abs(slice_index - ((s - 1) / switch(plane,
                                            "AXIAL" = dims[3] - 1,
                                            "SAGITTAL" = dims[1] - 1,
                                            "CORONAL" = dims[2] - 1))) < 0.02
      )

    validity <- map2_lgl(slice_gazes$gaze_x, slice_gazes$gaze_y,
                         ~resolve_coordinates(.x, .y)$in_bounds)

    data.frame(
      slice_num = s,
      total_gazes = nrow(slice_gazes),
      valid_gazes = sum(validity),
      invalid_gazes = sum(!validity),
      percent_valid = 100 * sum(validity) / nrow(slice_gazes)
    )
  })

  return(validity_summary)
}

# Example usage:
# Plot a single slice showing validity (legend on right)
# plot_slice_with_validity(nifti_data, integrated, slice_num = 12)

# Plot with legend at bottom (cleaner)
# plot_slice_with_validity_bottom_legend(nifti_data, integrated, slice_num = 12)

# Check validity for all slices
# validity_summary <- check_all_slice_validity(nifti_data, integrated)
# print(validity_summary)

# Plot only slices with high validity
# good_slices <- validity_summary %>%
#   filter(percent_valid > 80) %>%
#   pull(slice_num)
#
# for (s in good_slices[1:4]) {
#   plot_slice_with_validity_bottom_legend(nifti_data, integrated, s)
# }
