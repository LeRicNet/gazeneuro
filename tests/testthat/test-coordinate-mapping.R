# Tests for coordinate-mapping.R

# Helper to create mock nifti_data object for testing
create_mock_nifti_data <- function(dims = c(182, 218, 182)) {
  # Create a minimal NVImage-like object for testing
  mock_nvimage <- list(
    dimsRAS = c(3, dims[1], dims[2], dims[3]),
    pixDimsRAS = c(1, 1, 1, 1),
    convertFrac2MM = function(frac, isForceSliceMM = FALSE) {
      # Simple identity-like transformation for testing
      # frac [0,1] maps to mm space with origin at (-90, -126, -72)
      mm <- c(
        frac[1] * dims[1] - 90,
        frac[2] * dims[2] - 126,
        frac[3] * dims[3] - 72
      )
      return(c(mm, 1))
    }
  )

  return(list(
    nvimage = mock_nvimage,
    dims = dims
  ))
}

# ============================================================================
# Tests for get_slice_number
# ============================================================================

test_that("get_slice_number handles basic cases correctly", {
  # First slice
  expect_equal(get_slice_number(0, 100), 1L)

  # Last slice
  expect_equal(get_slice_number(1, 100), 100L)

  # Middle slice (50th of 100)
  expect_equal(get_slice_number(0.5, 99), 50L)
})

test_that("get_slice_number handles edge cases", {
  # Single slice
  expect_equal(get_slice_number(0, 1), 1L)
  expect_equal(get_slice_number(0.5, 1), 1L)
  expect_equal(get_slice_number(1, 1), 1L)

  # Two slices
  expect_equal(get_slice_number(0, 2), 1L)
  expect_equal(get_slice_number(0.49, 2), 1L)
  expect_equal(get_slice_number(0.51, 2), 2L)
  expect_equal(get_slice_number(1, 2), 2L)
})

test_that("get_slice_number respects rounding methods", {
  # With n_slices = 25, index 0.64 maps to 0.64 * 24 = 15.36 (0-based)
  # round: 15 + 1 = 16
  # floor: 15 + 1 = 16
  # ceiling: 16 + 1 = 17

  n <- 25
  idx <- 0.64

  expect_equal(get_slice_number(idx, n, method = "round"), 16L)
  expect_equal(get_slice_number(idx, n, method = "floor"), 16L)
  expect_equal(get_slice_number(idx, n, method = "ceiling"), 17L)
})

test_that("get_slice_number clamps out-of-range values", {
  # Negative index
  expect_equal(get_slice_number(-0.5, 100), 1L)

  # Index > 1
  expect_equal(get_slice_number(1.5, 100), 100L)
})

test_that("get_slice_number is consistent with 25-slice images", {
  # Common case: 25 slices in z-direction
  n <- 25

  # Index 0 -> slice 1
  expect_equal(get_slice_number(0, n), 1L)

  # Index 1 -> slice 25
  expect_equal(get_slice_number(1, n), 25L)

  # Index 0.5 -> slice 13 (middle)
  expect_equal(get_slice_number(0.5, n), 13L)

  # Test specific case that was failing:
  # If normalized_index maps to slice 64 but expected 64, that's correct
  # With 128 slices: index 0.5 -> 0.5 * 127 = 63.5 -> round to 64 -> slice 64
  expect_equal(get_slice_number(0.5, 128), 64L)
})

test_that("get_slice_number returns integer type", {
  result <- get_slice_number(0.5, 100)
  expect_type(result, "integer")
})

# ============================================================================
# Tests for get_normalized_index
# ============================================================================

test_that("get_normalized_index is inverse of get_slice_number", {
  n <- 25

  # Round-trip: slice -> index -> slice
  for (slice in 1:n) {
    idx <- get_normalized_index(slice, n)
    recovered_slice <- get_slice_number(idx, n)
    expect_equal(recovered_slice, slice)
  }
})

test_that("get_normalized_index handles boundaries", {
  expect_equal(get_normalized_index(1, 100), 0)
  expect_equal(get_normalized_index(100, 100), 1)
  expect_equal(get_normalized_index(50, 99), 0.5)
})

test_that("get_normalized_index handles single slice",
          expect_equal(get_normalized_index(1, 1), 0.5)
)

# ============================================================================
# Tests for locate_single_point
# ============================================================================

test_that("locate_single_point validates inputs", {
  expect_error(locate_single_point(0.5, 0.5, 10, NULL))
  expect_error(locate_single_point(0.5, 0.5, 10, list(nvimage = NULL)))
})

test_that("locate_single_point returns correct structure", {
  mock_data <- create_mock_nifti_data()

  result <- locate_single_point(0.5, 0.5, 10, mock_data, plane = "AXIAL")

  expect_type(result, "list")
  expect_named(result, c("voxel", "mm", "frac", "valid"))
  expect_length(result$voxel, 3)
  expect_length(result$mm, 3)
  expect_length(result$frac, 3)
  expect_type(result$valid, "logical")
})

test_that("locate_single_point handles center of image", {
  mock_data <- create_mock_nifti_data(dims = c(100, 100, 50))

  result <- locate_single_point(0.5, 0.5, 25, mock_data, plane = "AXIAL")

  # Center should be at voxel (50, 50, 25) approximately
  expect_equal(result$voxel[1], 50, tolerance = 1)
  expect_equal(result$voxel[2], 50, tolerance = 1)
  expect_equal(result$voxel[3], 25)  # Slice 25 (1-based)
  expect_true(result$valid)
})

test_that("locate_single_point handles different planes", {
  mock_data <- create_mock_nifti_data(dims = c(100, 100, 50))

  # AXIAL: gaze maps to X-Y, slice is Z
  axial <- locate_single_point(0.5, 0.5, 25, mock_data, plane = "AXIAL")
  expect_equal(axial$voxel[3], 25)

  # SAGITTAL: gaze maps to Y-Z, slice is X
  sagittal <- locate_single_point(0.5, 0.5, 50, mock_data, plane = "SAGITTAL")
  expect_equal(sagittal$voxel[1], 50)

  # CORONAL: gaze maps to X-Z, slice is Y
  coronal <- locate_single_point(0.5, 0.5, 50, mock_data, plane = "CORONAL")
  expect_equal(coronal$voxel[2], 50)
})

test_that("locate_single_point clamps gaze coordinates", {
  mock_data <- create_mock_nifti_data()

  # Out of bounds gaze should be clamped
  result_low <- locate_single_point(-0.5, -0.5, 10, mock_data)
  result_high <- locate_single_point(1.5, 1.5, 10, mock_data)

  expect_true(all(result_low$frac >= 0))
  expect_true(all(result_high$frac <= 1))
})

test_that("locate_single_point handles invalid plane", {
  mock_data <- create_mock_nifti_data()
  expect_error(locate_single_point(0.5, 0.5, 10, mock_data, plane = "INVALID"))
})

# ============================================================================
# Tests for map_gaze_to_voxels
# ============================================================================

test_that("map_gaze_to_voxels handles batch processing", {
  mock_data <- create_mock_nifti_data(dims = c(100, 100, 25))

  gaze_df <- data.frame(
    gaze_x = c(0.25, 0.5, 0.75),
    gaze_y = c(0.25, 0.5, 0.75)
  )

  slice_df <- data.frame(
    slice_num = c(5, 13, 20)
  )

  result <- map_gaze_to_voxels(gaze_df, slice_df, mock_data)

  expect_equal(nrow(result), 3)
  expect_true("voxel_i" %in% names(result))
  expect_true("voxel_j" %in% names(result))
  expect_true("voxel_k" %in% names(result))
  expect_true("mm_x" %in% names(result))
  expect_true("mm_y" %in% names(result))
  expect_true("mm_z" %in% names(result))
  expect_true("valid" %in% names(result))
})

test_that("map_gaze_to_voxels validates input lengths", {
  mock_data <- create_mock_nifti_data()

  gaze_df <- data.frame(gaze_x = 1:5, gaze_y = 1:5)
  slice_df <- data.frame(slice_num = 1:3)  # Wrong length

  expect_error(map_gaze_to_voxels(gaze_df, slice_df, mock_data))
})

# ============================================================================
# Tests for get_slice_transitions
# ============================================================================

test_that("get_slice_transitions identifies transitions", {
  z_axis <- data.frame(
    time_sec = c(0, 1, 2, 3, 4),
    slice_num = c(1, 1, 2, 2, 3)
  )

  transitions <- get_slice_transitions(z_axis)

  expect_equal(nrow(transitions), 2)  # Two transitions: 1->2 and 2->3
  expect_equal(transitions$transition_time, c(2, 4))
  expect_equal(transitions$slice_before, c(1, 2))
  expect_equal(transitions$slice_after, c(2, 3))
})

test_that("get_slice_transitions handles no transitions", {
  z_axis <- data.frame(
    time_sec = c(0, 1, 2),
    slice_num = c(5, 5, 5)
  )

  transitions <- get_slice_transitions(z_axis)

  expect_equal(nrow(transitions), 0)
})

test_that("get_slice_transitions handles minimal data", {
  # Single row
  z_axis_single <- data.frame(time_sec = 0, slice_num = 1)
  expect_equal(nrow(get_slice_transitions(z_axis_single)), 0)

  # Empty
  z_axis_empty <- data.frame(time_sec = numeric(0), slice_num = integer(0))
  expect_equal(nrow(get_slice_transitions(z_axis_empty)), 0)
})

test_that("get_slice_transitions handles unsorted input", {
  z_axis <- data.frame(
    time_sec = c(3, 1, 4, 2, 0),  # Out of order
    slice_num = c(2, 1, 3, 2, 1)
  )

  transitions <- get_slice_transitions(z_axis)

  # After sorting: t=0(1), t=1(1), t=2(2), t=3(2), t=4(3)
  # Transitions at t=2 (1->2) and t=4 (2->3)
  expect_equal(transitions$transition_time, c(2, 4))
})
