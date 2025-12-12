# 05_debug_transformations.R
# Debug the coordinate transformation validation to identify the 0.866 voxel error

library(qs)
library(tidyverse)
library(gazeneuro)
library(RNifti)

# Load test data
volumes <- qread("coordinate_validation/synthetic_nifti_volumes.qs")

# Start with the simplest case - RAS_isotropic
test_volume <- volumes[["RAS_isotropic"]]
temp_file <- tempfile(fileext = ".nii.gz")
writeNifti(test_volume$nifti, temp_file)
nifti_data <- preload_nifti_data(temp_file)

cat("=== DEBUGGING COORDINATE TRANSFORMATIONS ===\n\n")
cat("Test volume: RAS_isotropic\n")
cat("Dimensions:", test_volume$test_case$dims, "\n")
cat("Pixel dimensions:", test_volume$test_case$pixdim, "\n")
cat("Expected affine:\n")
print(test_volume$affine)
cat("\nActual affine from NIfTI:\n")
print(xform(test_volume$nifti))

# Test 1: Simple center point transformation
cat("\n\n=== TEST 1: Center Point ===\n")
dims <- test_volume$test_case$dims
center_voxel <- (dims - 1) / 2
cat("Center voxel (0-based):", center_voxel, "\n")

# Manual transformation using affine
affine <- xform(test_volume$nifti)
center_world_manual <- affine %*% c(center_voxel, 1)
cat("Center world (manual):", center_world_manual[1:3], "\n")

# Using NVImage
nvimg <- nifti_data$nvimage
center_frac <- center_voxel / nvimg$dimsRAS[2:4]
center_world_nvimage <- nvimg$convertFrac2MM(center_frac)[1:3]
cat("Center world (NVImage):", center_world_nvimage, "\n")
cat("Difference:", center_world_nvimage - center_world_manual[1:3], "\n")

# Test 2: Round-trip for center point
cat("\n\n=== TEST 2: Round-trip Center Point ===\n")
# Forward: voxel -> world
voxel_orig <- center_voxel
world_coords <- nvimg$convertFrac2MM(voxel_orig / nvimg$dimsRAS[2:4])[1:3]
cat("Original voxel:", voxel_orig, "\n")
cat("World coords:", world_coords, "\n")

# Backward: world -> voxel
voxel_back <- nvimg$mm2vox(world_coords, frac = TRUE)
cat("Voxel back:", voxel_back, "\n")
cat("Round-trip error:", sqrt(sum((voxel_orig - voxel_back)^2)), "\n")

# Test 3: Check for off-by-one or half-voxel issues
cat("\n\n=== TEST 3: Grid of Test Points ===\n")
test_points <- expand.grid(
  x = c(0, dims[1]/2 - 0.5, dims[1] - 1),
  y = c(0, dims[2]/2 - 0.5, dims[2] - 1),
  z = c(0, dims[3]/2 - 0.5, dims[3] - 1)
)

errors <- numeric(nrow(test_points))
for (i in 1:nrow(test_points)) {
  voxel <- as.numeric(test_points[i,])
  frac <- voxel / nvimg$dimsRAS[2:4]
  world <- nvimg$convertFrac2MM(frac)[1:3]
  voxel_back <- nvimg$mm2vox(world, frac = TRUE)
  errors[i] <- sqrt(sum((voxel - voxel_back)^2))
}

test_points$error <- errors
print(test_points)
cat("\nMean error:", mean(errors), "\n")
cat("All errors equal?", length(unique(errors)) == 1, "\n")

# Test 4: Investigate the 0.866 pattern
cat("\n\n=== TEST 4: Investigating 0.866 Pattern ===\n")
cat("sqrt(3)/2 =", sqrt(3)/2, "\n")
cat("sqrt(0.75) =", sqrt(0.75), "\n")
cat("Is error exactly sqrt(3)/2?", abs(mean(errors) - sqrt(3)/2) < 1e-10, "\n")

# Test 5: Check NVImage dimensions vs actual dimensions
cat("\n\n=== TEST 5: Dimension Mismatch Check ===\n")
cat("Test case dims:", dims, "\n")
cat("NVImage dimsRAS:", nvimg$dimsRAS, "\n")
cat("Dimension match?", all(dims == nvimg$dimsRAS[2:4]), "\n")

# Test 6: Trace through a single transformation step by step
cat("\n\n=== TEST 6: Step-by-Step Transformation ===\n")
test_voxel <- c(100, 100, 50)
cat("Test voxel:", test_voxel, "\n")

# Step 1: Voxel to fractional
frac <- test_voxel / nvimg$dimsRAS[2:4]
cat("Fractional coords:", frac, "\n")

# Step 2: Check frac2mm matrix
cat("\nfrac2mm matrix:\n")
print(nvimg$frac2mm[1:4, 1:4])

# Step 3: Transform to world
world <- nvimg$convertFrac2MM(frac, isForceSliceMM = TRUE)[1:3]
cat("\nWorld coords (slice):", world, "\n")

world_ortho <- nvimg$convertFrac2MM(frac, isForceSliceMM = FALSE)[1:3]
cat("World coords (ortho):", world_ortho, "\n")

# Step 4: Back to voxel
inv_affine <- solve(nvimg$matRAS)
voxel_back_manual <- inv_affine %*% c(world, 1)
cat("\nVoxel back (manual):", voxel_back_manual[1:3], "\n")

voxel_back_nvimage <- nvimg$mm2vox(world, frac = TRUE)
cat("Voxel back (NVImage):", voxel_back_nvimage, "\n")

cat("\nError:", sqrt(sum((test_voxel - voxel_back_nvimage)^2)), "\n")

# Test 7: Check if the issue is with the gaze pattern generation
cat("\n\n=== TEST 7: Gaze Pattern Validation ===\n")
gaze_patterns <- qread("coordinate_validation/synthetic_gaze_patterns.qs")

# Get one fixation point for RAS_isotropic
test_gaze <- gaze_patterns %>%
  filter(volume_name == "RAS_isotropic",
         pattern_name == "target_fixations",
         plane == "AXIAL") %>%
  slice(1)

if (nrow(test_gaze) > 0) {
  cat("Test gaze point:\n")
  cat("  Screen coords:", test_gaze$gaze_x, test_gaze$gaze_y, "\n")
  cat("  Plane:", test_gaze$plane, "Slice:", test_gaze$slice_idx, "\n")

  # Reconstruct voxel coordinates as done in validation
  if (test_gaze$plane == "AXIAL") {
    voxel_x <- test_gaze$gaze_x * (dims[1] - 1)
    voxel_y <- test_gaze$gaze_y * (dims[2] - 1)
    voxel_z <- test_gaze$slice_idx
  }

  cat("\nReconstructed voxel:", voxel_x, voxel_y, voxel_z, "\n")

  # Check if this matches expected target location
  if (!is.na(test_gaze$target_desc)) {
    cat("Target description:", test_gaze$target_desc, "\n")
  }
}

# Clean up
unlink(temp_file)

cat("\n\n=== DEBUGGING SUMMARY ===\n")
cat("The consistent 0.866 error (sqrt(3)/2) suggests a systematic offset.\n")
cat("Check the output above for dimension mismatches or transformation issues.\n")
