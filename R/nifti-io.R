#' Load NIfTI file with pre-extracted data
#'
#' @description
#' Loads a NIfTI file and pre-extracts all image data for faster access.
#' Also creates an NVImage object for coordinate transformations.
#'
#' @param filepath Path to the NIfTI file
#' @return A list containing:
#'   \item{nifti}{The RNifti object}
#'   \item{nvimage}{An NVImage object with transformation matrices}
#'   \item{data}{Pre-extracted image data as array}
#'   \item{dims}{Dimensions of the image}
#'
#' @export
#' @importFrom RNifti readNifti xform pixdim
preload_nifti_data <- function(filepath) {
  message("Loading NIfTI file: ", basename(filepath))

  # Load with RNifti
  nii <- RNifti::readNifti(filepath)

  # Extract all data as array
  message("Extracting image data...")
  img_data <- as.array(nii)

  # Create NVImage
  nvimg <- NVImage$new(
    matRAS = RNifti::xform(nii),
    dimsRAS = c(length(dim(nii)), dim(nii)),
    pixDimsRAS = c(1, RNifti::pixdim(nii))
  )
  nvimg$calculateOblique()

  # Return with pre-extracted data
  return(list(
    nifti = nii,
    nvimage = nvimg,
    data = img_data,
    dims = dim(nii)
  ))
}

#' Safely extract a slice from NIfTI data
#'
#' @param nii NIfTI object
#' @param slice_num Slice number to extract
#' @return Matrix of slice data or NULL if error
#' @export
safe_get_slice <- function(nii, slice_num) {
  tryCatch({
    # Convert to array first to ensure data is accessible
    img_array <- as.array(nii)
    return(img_array[,, slice_num])
  }, error = function(e) {
    message("Error accessing slice data: ", e$message)
    return(NULL)
  })
}

#' Quick view of NIfTI slice
#'
#' @description
#' Displays a single slice from NIfTI data with proper aspect ratio
#' and coordinate information.
#'
#' @param data List from preload_nifti_data()
#' @param slice Slice number to display (NULL for middle slice)
#' @return Invisible slice data matrix
#' @export
#' @importFrom graphics image axis box mtext par
safe_quick_view <- function(data, slice = NULL) {
  nii <- data$nifti
  nvimg <- data$nvimage

  # Get dimensions safely
  dims <- dim(nii)

  if (is.null(slice)) {
    slice <- ceiling(dims[3] / 2)  # Middle slice
  }

  # Validate slice number
  if (slice < 1 || slice > dims[3]) {
    message("Invalid slice number. Must be between 1 and ", dims[3])
    return(invisible(NULL))
  }

  # Get slice data safely
  slice_data <- safe_get_slice(nii, slice)

  if (is.null(slice_data)) {
    message("Could not extract slice data")
    return(invisible(NULL))
  }

  # Set up plot
  par(mar = c(4, 4, 3, 2))

  # Display with proper aspect ratio
  asp <- nvimg$pixDimsRAS[3] / nvimg$pixDimsRAS[2]

  # Use matrix coordinates for display
  image(slice_data,
        asp = asp,
        col = gray((0:255)/255),
        xlab = "X (fraction)", ylab = "Y (fraction)",
        main = sprintf("Slice %d of %d", slice, dims[3]),
        axes = FALSE)

  # Add custom axes showing actual dimensions
  axis(1, at = seq(0, 1, 0.25), labels = round(seq(0, dims[1], dims[1]/4)))
  axis(2, at = seq(0, 1, 0.25), labels = round(seq(0, dims[2], dims[2]/4)))
  box()

  # Show world coordinates for center
  center_voxel <- c(dims[1]/2 - 1, dims[2]/2 - 1, slice - 1)  # 0-based
  center_frac <- center_voxel / nvimg$dimsRAS[2:4]
  center_world <- nvimg$convertFrac2MM(center_frac)[1:3]

  mtext(sprintf("Center: World [%.1f, %.1f, %.1f] mm",
                center_world[1], center_world[2], center_world[3]),
        side = 1, line = 3, cex = 0.8)

  return(invisible(slice_data))
}

#' Fast slice viewer using pre-extracted data
#'
#' @param data List from preload_nifti_data()
#' @param slice Slice number to display
#' @export
#' @importFrom graphics image title par
fast_slice_viewer <- function(data, slice = NULL) {
  if (is.null(data$data)) {
    message("No pre-extracted data found. Use preload_nifti_data() first.")
    return(invisible(NULL))
  }

  nvimg <- data$nvimage
  img_data <- data$data
  dims <- data$dims

  if (is.null(slice)) slice <- ceiling(dims[3] / 2)

  # Get slice directly from array
  slice_data <- img_data[,, slice]

  # Create the plot
  par(mar = c(5, 4, 4, 2))
  image(slice_data,
        col = gray((0:255)/255),
        main = sprintf("Slice %d of %d", slice, dims[3]),
        xlab = "X", ylab = "Y")

  # Add coordinate info
  center_voxel <- c(dims[1]/2 - 1, dims[2]/2 - 1, slice - 1)
  center_frac <- center_voxel / nvimg$dimsRAS[2:4]
  center_world <- nvimg$convertFrac2MM(center_frac)[1:3]

  title(sub = sprintf("Center: [%.1f, %.1f, %.1f] mm",
                      center_world[1], center_world[2], center_world[3]),
        cex.sub = 0.8)
}

#' View multiple NIfTI slices in a grid
#'
#' @param data List from preload_nifti_data()
#' @param slices Vector of slice numbers (NULL for equally spaced)
#' @export
#' @importFrom graphics image box par
view_multiple_slices <- function(data, slices = NULL) {
  if (is.null(data$data)) {
    data <- preload_nifti_data(data$nifti@.xData$.sourcePath)
  }

  dims <- data$dims

  # Default to equally spaced slices
  if (is.null(slices)) {
    n_slices <- min(9, dims[3])
    slices <- round(seq(1, dims[3], length.out = n_slices))
  }

  # Set up multi-panel plot
  n <- length(slices)
  nr <- ceiling(sqrt(n))
  nc <- ceiling(n / nr)

  par(mfrow = c(nr, nc), mar = c(2, 2, 3, 1))

  for (s in slices) {
    if (s >= 1 && s <= dims[3]) {
      slice_data <- data$data[,, s]
      image(slice_data,
            col = gray((0:255)/255),
            main = sprintf("Slice %d", s),
            axes = FALSE)
      box()
    }
  }

  par(mfrow = c(1, 1))
}

#' Interactive coordinate picker for NIfTI slices
#'
#' @description
#' Allows interactive clicking on a NIfTI slice to get voxel and world coordinates.
#' Right-click to exit.
#'
#' @param data List from preload_nifti_data()
#' @param slice Slice number to display
#' @return Data frame of clicked points with coordinates
#' @export
#' @importFrom graphics locator points text
safe_coordinate_picker <- function(data, slice = NULL) {
  nii <- data$nifti
  nvimg <- data$nvimage
  dims <- dim(nii)

  if (is.null(slice)) slice <- ceiling(dims[3] / 2)

  # Display the slice
  slice_data <- safe_quick_view(data, slice)

  if (is.null(slice_data)) {
    return(invisible(NULL))
  }

  message("\nClick on the image to see coordinates (right-click to exit)")
  message("Note: Click positions are in fractional coordinates (0-1)\n")

  points_clicked <- data.frame()

  repeat {
    click <- locator(1)
    if (is.null(click)) break

    # Convert from plot coordinates (0-1) to voxel coordinates
    voxel_x <- round(click$x * (dims[1] - 1))
    voxel_y <- round(click$y * (dims[2] - 1))
    voxel_z <- slice - 1  # 0-based for calculations

    # Bounds check
    if (voxel_x >= 0 && voxel_x < dims[1] &&
        voxel_y >= 0 && voxel_y < dims[2]) {

      # Convert to world coordinates
      frac <- c(voxel_x, voxel_y, voxel_z) / nvimg$dimsRAS[2:4]
      world <- nvimg$convertFrac2MM(frac)[1:3]

      # Get intensity safely
      intensity <- tryCatch({
        slice_data[voxel_x + 1, voxel_y + 1]
      }, error = function(e) NA)

      # Display info
      message(sprintf("Click %d: Voxel [%d,%d,%d] → World [%.1f,%.1f,%.1f] mm",
                      nrow(points_clicked) + 1,
                      voxel_x, voxel_y, voxel_z,
                      world[1], world[2], world[3]), appendLF = FALSE)

      if (!is.na(intensity)) {
        message(sprintf(", Intensity: %.0f", intensity))
      } else {
        message("")
      }

      # Mark the point
      points(click$x, click$y, col = "red", pch = 3, cex = 2, lwd = 2)
      text(click$x, click$y, labels = nrow(points_clicked) + 1,
           col = "yellow", pos = 4, cex = 0.8)

      # Store the point
      points_clicked <- rbind(points_clicked, data.frame(
        click_num = nrow(points_clicked) + 1,
        voxel_x = voxel_x,
        voxel_y = voxel_y,
        voxel_z = voxel_z,
        world_x = world[1],
        world_y = world[2],
        world_z = world[3],
        intensity = intensity
      ))
    } else {
      message("Click outside image bounds")
    }
  }

  if (nrow(points_clicked) > 0) {
    message("\nSummary of clicked points:")
    print(points_clicked)
  }

  return(invisible(points_clicked))
}
