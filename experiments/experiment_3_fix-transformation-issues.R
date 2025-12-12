# 06_fix_transformation_issues.R
# Diagnose and fix the coordinate transformation issues

library(RNifti)
library(qs)
library(tidyverse)
library(gazeneuro)

cat("=== IDENTIFYING THE ROOT CAUSE ===\n\n")

# Load the test volume
volumes <- qread("coordinate_validation/synthetic_nifti_volumes.qs")
test_volume <- volumes[["RAS_isotropic"]]

cat("Expected affine from test generation:\n")
print(test_volume$affine)

# Check what's actually in the NIfTI
cat("\nxform() output:\n")
print(xform(test_volume$nifti))

# The xform is returning identity with attributes
# Let's check if we can access the qform/sform directly
cat("\nChecking NIfTI internals:\n")
cat("qform_code:", test_volume$nifti$qform_code, "\n")
cat("sform_code:", test_volume$nifti$sform_code, "\n")

# The real issue: RNifti might be using a different convention
# Let's test if the problem is in our NIfTI generation

cat("\n\n=== TESTING AFFINE PRESERVATION ===\n")

# Create a simple test
test_data <- array(1:64, dim = c(4, 4, 4))
test_affine <- matrix(c(
  2, 0, 0, -4,
  0, 2, 0, -4,
  0, 0, 2, -4,
  0, 0, 0, 1
), nrow = 4, byrow = TRUE)

# Try different ways to set the affine
cat("\nTest affine:\n")
print(test_affine)

# Method 1: Direct xform parameter
nii1 <- asNifti(test_data, xform = test_affine)
xf1 <- xform(nii1)
cat("\nMethod 1 - xform parameter:\n")
print(xf1)

# Save and reload to see if it persists
temp_file <- tempfile(fileext = ".nii.gz")
writeNifti(nii1, temp_file)
nii1_reload <- readNifti(temp_file)
cat("\nAfter save/reload:\n")
print(xform(nii1_reload))
unlink(temp_file)

cat("\n\n=== THE REAL PROBLEM ===\n")

# The issue is likely that our validation is comparing against the wrong thing
# Let's trace through what's actually happening in the transformation

# 1. In the gaze pattern generation, we used world_to_screen which assumes
#    the affine is properly set
# 2. But when we load the NIfTI, the affine is identity
# 3. So the NVImage is using an identity transform

# Let's check what NVImage is actually doing
temp_file <- tempfile(fileext = ".nii.gz")
writeNifti(test_volume$nifti, temp_file)
nifti_data <- preload_nifti_data(temp_file)

cat("\nNVImage matRAS:\n")
print(nifti_data$nvimage$matRAS[1:4, 1:4])

cat("\nNVImage dimsRAS:", nifti_data$nvimage$dimsRAS, "\n")

cat("\nNVImage frac2mm:\n")
print(nifti_data$nvimage$frac2mm[1:4, 1:4])

unlink(temp_file)

cat("\n\n=== THE SOLUTION ===\n")

# The problem is a mismatch between:
# 1. How we generated the gaze patterns (assuming proper affine)
# 2. How the validation tests them (using identity affine)

# Two possible fixes:
cat("\nOption 1: Fix the NIfTI generation to properly save affines\n")
cat("Option 2: Fix the validation to account for identity affines\n")
cat("Option 3: Fix the gaze pattern generation to match what we actually get\n")

# Let's implement Option 3 - regenerate gaze patterns using actual transforms
cat("\n\nThe 0.866 error is because:\n")
cat("- Expected transform includes translation of (-128, -128, -64)\n")
cat("- Actual transform is identity (0, 0, 0)\n")
cat("- The offset is (128, 128, 64) in voxel space\n")
cat("- But we're testing fractional coordinates around (127.5, 127.5, 63.5)\n")
cat("- The error is sqrt(0.5² + 0.5² + 0.5²) = 0.866\n")

cat("\n\n=== QUICK FIX FOR VALIDATION ===\n")

# Create a function that does round-trip with identity affine
test_round_trip_identity <- function(screen_x, screen_y, plane, slice_idx, dims) {
  # Screen to voxel
  if (plane == "AXIAL") {
    voxel <- c(screen_x * (dims[1] - 1),
               screen_y * (dims[2] - 1),
               slice_idx)
  }

  # With identity affine, voxel coords ARE world coords
  world <- voxel

  # Back to voxel (no change with identity)
  voxel_back <- world

  # Back to screen
  if (plane == "AXIAL") {
    screen_back <- c(voxel_back[1] / (dims[1] - 1),
                     voxel_back[2] / (dims[2] - 1))
  }

  # Errors should be zero with identity transform
  return(list(
    voxel_error = sqrt(sum((voxel - voxel_back)^2)),
    screen_error = sqrt(sum((c(screen_x, screen_y) - screen_back)^2))
  ))
}

# Test it
result <- test_round_trip_identity(0.5, 0.5, "AXIAL", 64, c(256, 256, 128))
cat("\nIdentity transform round-trip errors:\n")
cat("  Voxel error:", result$voxel_error, "\n")
cat("  Screen error:", result$screen_error, "\n")

cat("\nConclusion: The validation expects non-identity affines but gets identity.\n")
cat("We need to either fix the NIfTI generation or adjust the validation logic.\n")
