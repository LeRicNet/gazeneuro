# 01_generate_test_nifti.R
# Generate synthetic NIfTI volumes with diverse affine transformations
# for coordinate transformation validation

library(RNifti)
library(qs)
library(tidyverse)

# Create output directory for intermediate files
dir.create("coordinate_validation", showWarnings = FALSE)

#' Generate synthetic brain-like volume
#' @param dims Dimensions of the volume (x, y, z)
#' @return 3D array with synthetic brain-like pattern
generate_synthetic_brain <- function(dims = c(256, 256, 128)) {
  # Create coordinate grids
  x <- seq(-1, 1, length.out = dims[1])
  y <- seq(-1, 1, length.out = dims[2])
  z <- seq(-0.5, 0.5, length.out = dims[3])

  # Initialize volume
  volume <- array(0, dim = dims)

  # Create ellipsoid mask (brain shape)
  for (i in 1:dims[1]) {
    for (j in 1:dims[2]) {
      for (k in 1:dims[3]) {
        # Ellipsoid equation
        if ((x[i]^2 / 0.64) + (y[j]^2 / 0.49) + (z[k]^2 / 0.36) < 1) {
          # Add some structure
          volume[i, j, k] <- 100 + 50 * sin(5 * x[i]) * cos(3 * y[j]) +
            20 * sin(10 * z[k])
        }
      }
    }
  }

  # Add noise
  volume <- volume + rnorm(prod(dims), 0, 5)
  volume[volume < 0] <- 0

  return(volume)
}

#' Create affine matrix for specified orientation
#' @param orientation One of "RAS", "LPS", "LAS", "RPI", etc.
#' @param pixdim Voxel dimensions (x, y, z) in mm
#' @param origin Origin coordinates in mm
#' @return 4x4 affine transformation matrix
create_orientation_matrix <- function(orientation, pixdim = c(1, 1, 2),
                                      origin = c(-128, -128, -64)) {
  # Base matrix (identity)
  mat <- diag(4)

  # Apply voxel dimensions
  mat[1:3, 1:3] <- diag(pixdim)

  # Apply orientation flips
  if (substr(orientation, 1, 1) == "L") mat[1, 1] <- -mat[1, 1]  # Left
  if (substr(orientation, 2, 2) == "P") mat[2, 2] <- -mat[2, 2]  # Posterior
  if (substr(orientation, 3, 3) == "I") mat[3, 3] <- -mat[3, 3]  # Inferior

  # Set origin
  mat[1:3, 4] <- origin

  return(mat)
}

#' Create oblique affine matrix
#' @param base_orientation Base orientation
#' @param rotation_angles Rotation angles (x, y, z) in radians
#' @param pixdim Voxel dimensions
#' @param origin Origin coordinates
#' @return 4x4 oblique affine transformation matrix
create_oblique_matrix <- function(base_orientation = "RAS",
                                  rotation_angles = c(0.1, 0.05, 0.15),
                                  pixdim = c(1, 1, 2),
                                  origin = c(-128, -128, -64)) {
  # Start with base orientation
  mat <- create_orientation_matrix(base_orientation, pixdim, origin)

  # Create rotation matrices
  Rx <- matrix(c(1, 0, 0, 0,
                 0, cos(rotation_angles[1]), -sin(rotation_angles[1]), 0,
                 0, sin(rotation_angles[1]), cos(rotation_angles[1]), 0,
                 0, 0, 0, 1), 4, 4, byrow = TRUE)

  Ry <- matrix(c(cos(rotation_angles[2]), 0, sin(rotation_angles[2]), 0,
                 0, 1, 0, 0,
                 -sin(rotation_angles[2]), 0, cos(rotation_angles[2]), 0,
                 0, 0, 0, 1), 4, 4, byrow = TRUE)

  Rz <- matrix(c(cos(rotation_angles[3]), -sin(rotation_angles[3]), 0, 0,
                 sin(rotation_angles[3]), cos(rotation_angles[3]), 0, 0,
                 0, 0, 1, 0,
                 0, 0, 0, 1), 4, 4, byrow = TRUE)

  # Apply rotations
  mat <- mat %*% Rx %*% Ry %*% Rz

  return(mat)
}

# Generate test cases
test_cases <- list()
case_id <- 1

# Standard clinical orientations
clinical_orientations <- c("RAS", "LPS", "LAS", "RPI")
for (orient in clinical_orientations) {
  # Isotropic voxels
  test_cases[[case_id]] <- list(
    id = case_id,
    name = paste0(orient, "_isotropic"),
    orientation = orient,
    dims = c(256, 256, 128),
    pixdim = c(1, 1, 1),
    is_oblique = FALSE
  )
  case_id <- case_id + 1

  # Anisotropic voxels (common in clinical scans)
  test_cases[[case_id]] <- list(
    id = case_id,
    name = paste0(orient, "_anisotropic"),
    orientation = orient,
    dims = c(256, 256, 64),
    pixdim = c(0.9375, 0.9375, 3.0),  # Common clinical resolution
    is_oblique = FALSE
  )
  case_id <- case_id + 1
}

# Oblique acquisitions
oblique_angles <- list(
  small = c(0.05, 0.05, 0.05),    # ~3 degrees
  moderate = c(0.15, 0.1, 0.2),   # ~9-11 degrees
  large = c(0.3, 0.25, 0.35)      # ~17-20 degrees
)

for (angle_name in names(oblique_angles)) {
  test_cases[[case_id]] <- list(
    id = case_id,
    name = paste0("oblique_", angle_name),
    orientation = "RAS",
    dims = c(256, 256, 128),
    pixdim = c(1, 1, 1.5),
    rotation_angles = oblique_angles[[angle_name]],
    is_oblique = TRUE
  )
  case_id <- case_id + 1
}

# Edge cases
# Very small volume
test_cases[[case_id]] <- list(
  id = case_id,
  name = "small_volume",
  orientation = "RAS",
  dims = c(64, 64, 32),
  pixdim = c(3, 3, 5),
  is_oblique = FALSE
)
case_id <- case_id + 1

# Very anisotropic
test_cases[[case_id]] <- list(
  id = case_id,
  name = "very_anisotropic",
  orientation = "LPS",
  dims = c(512, 512, 25),
  pixdim = c(0.5, 0.5, 7.0),
  is_oblique = FALSE
)
case_id <- case_id + 1

# Generate NIfTI files for each test case
results <- list()
for (tc in test_cases) {
  message(sprintf("Generating test case %d: %s", tc$id, tc$name))

  # Generate volume data
  volume <- generate_synthetic_brain(tc$dims)

  # Create affine matrix
  if (tc$is_oblique) {
    affine <- create_oblique_matrix(
      tc$orientation,
      tc$rotation_angles,
      tc$pixdim,
      origin = -tc$pixdim * tc$dims / 2  # Center at origin
    )
  } else {
    affine <- create_orientation_matrix(
      tc$orientation,
      tc$pixdim,
      origin = -tc$pixdim * tc$dims / 2
    )
  }

  # Create NIfTI object with proper affine transformation
  # RNifti's asNifti can accept an xform parameter
  nii <- asNifti(volume, xform = affine, pixdim = c(1, tc$pixdim))

  # Verify the affine was set correctly
  stored_affine <- xform(nii)
  if (!isTRUE(all.equal(affine, stored_affine, tolerance = 1e-6))) {
    warning(sprintf("Affine matrix mismatch for %s", tc$name))
  }

  # Store test case info
  tc$affine <- affine
  tc$volume_stats <- list(
    min = min(volume),
    max = max(volume),
    mean = mean(volume),
    sd = sd(volume)
  )

  # Save to results
  results[[tc$name]] <- list(
    test_case = tc,
    nifti = nii,
    affine = affine
  )
}

# Save results
qsave(results, "coordinate_validation/synthetic_nifti_volumes.qs")

# Save summary
summary_df <- map_df(test_cases, function(tc) {
  data.frame(
    id = tc$id,
    name = tc$name,
    orientation = tc$orientation,
    dims = paste(tc$dims, collapse = "x"),
    pixdim = paste(tc$pixdim, collapse = "x"),
    is_oblique = tc$is_oblique,
    voxel_count = prod(tc$dims)
  )
})

write_csv(summary_df, "coordinate_validation/test_cases_summary.csv")

# Print summary
cat("\nGenerated", length(test_cases), "test NIfTI volumes:\n")
print(summary_df)

cat("\nTest data saved to coordinate_validation/\n")
cat("- synthetic_nifti_volumes.qs: Complete NIfTI data and metadata\n")
cat("- test_cases_summary.csv: Summary of test cases\n")
