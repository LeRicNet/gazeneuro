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

  # Create the plot
  par(mar = c(5, 4, 4, 2))

  # Plot the brain slice
  image(slice_data, col = gray((0:255)/255),
        main = sprintf("%s Slice %d - Gaze Validity Check", plane, slice_num),
        xlab = "X", ylab = "Y", xlim = c(0, 1), ylim = c(0, 1))

  # Show frame boundaries if requested
  if (show_frame) {
    # Canvas boundaries in normalized coordinates
    left_edge <- 350 / 2624   # ~0.133
    right_edge <- (350 + 1924) / 2624  # ~0.867
    top_edge <- 80 / 1640     # ~0.049
    bottom_edge <- (80 + 1560) / 1640  # ~0.951

    # Draw frame outline (in original Tobii coordinates)
    rect(0, 0, 1, 1, border = "gray", lty = 2, lwd = 2)

    # Draw canvas outline (mapped to display coordinates)
    # Note: these are approximate visual guides
    rect(left_edge, top_edge, right_edge, bottom_edge,
         border = "blue", lty = 1, lwd = 2)

    # Add labels
    text(0.5, 0.02, "Full Display Frame", cex = 0.8, col = "gray")
    text(0.5, 0.95, "Image Canvas", cex = 0.8, col = "blue")
  }

  # Plot out-of-bounds points (original positions)
  out_bounds <- validity_check[!validity_check$in_bounds, ]
  if (nrow(out_bounds) > 0) {
    points(out_bounds$x_original, out_bounds$y_original,
           pch = 4, col = "red", cex = 1.2, lwd = 2)  # X marks
  }

  # Plot in-bounds points (adjusted positions)
  in_bounds <- validity_check[validity_check$in_bounds, ]
  if (nrow(in_bounds) > 0) {
    points(in_bounds$x_adjusted, in_bounds$y_adjusted,
           pch = 19, col = rgb(0, 1, 0, 0.6), cex = 0.8)  # Green dots

    # Show mean of valid points
    mean_x <- mean(in_bounds$x_adjusted)
    mean_y <- mean(in_bounds$y_adjusted)
    points(mean_x, mean_y, pch = 3, col = "blue", cex = 2, lwd = 3)
  }

  # Add legend
  legend("topright",
         legend = c(
           sprintf("Valid gazes: %d", nrow(in_bounds)),
           sprintf("Out of bounds: %d", nrow(out_bounds)),
           "Mean position"
         ),
         pch = c(19, 4, 3),
         col = c(rgb(0, 1, 0, 0.6), "red", "blue"),
         bg = "white",
         cex = 0.8)

  # Print summary
  cat(sprintf("\nSlice %d summary:\n", slice_num))
  cat(sprintf("  Total gaze points: %d\n", nrow(slice_gazes)))
  cat(sprintf("  Valid (in canvas): %d (%.1f%%)\n",
              nrow(in_bounds), 100 * nrow(in_bounds) / nrow(slice_gazes)))
  cat(sprintf("  Out of bounds: %d (%.1f%%)\n",
              nrow(out_bounds), 100 * nrow(out_bounds) / nrow(slice_gazes)))

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
# Plot a single slice showing validity
# plot_slice_with_validity(nifti_data, integrated, slice_num = 12)

# Check validity for all slices
# validity_summary <- check_all_slice_validity(nifti_data, integrated)
# print(validity_summary)

# Plot only slices with high validity
# good_slices <- validity_summary %>%
#   filter(percent_valid > 80) %>%
#   pull(slice_num)
#
# for (s in good_slices[1:4]) {
#   plot_slice_with_validity(nifti_data, integrated, s)
# }
