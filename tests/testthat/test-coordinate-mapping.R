test_that("coordinate mapping works correctly", {
  # Create test NIfTI data
  test_dims <- c(64, 64, 32)
  test_mat <- diag(4)
  test_mat[1:3, 4] <- c(-32, -32, -16)  # Center at origin

  test_nifti_data <- list(
    nvimage = NVImage$new(
      matRAS = test_mat,
      dimsRAS = c(3, test_dims[1], test_dims[2], test_dims[3]),
      pixDimsRAS = c(1, 1, 1, 1)
    ),
    data = array(rnorm(prod(test_dims)), dim = test_dims),
    dims = test_dims
  )
  test_nifti_data$nvimage$calculateOblique()

  # Test single point location
  location <- locate_single_point(
    x = 0.5, y = 0.5, z = 0.5,
    nifti_data = test_nifti_data,
    plane = "AXIAL"
  )

  expect_type(location, "list")
  expect_equal(length(location$mm), 3)
  expect_equal(length(location$vox), 3)
  expect_equal(location$frac, c(0.5, 0.5, 0.5))

  # Test that center voxel maps to near origin
  expect_lt(abs(location$mm[1]), 2)
  expect_lt(abs(location$mm[2]), 2)
  expect_lt(abs(location$mm[3]), 2)
})

test_that("slice number conversion works", {
  dims <- c(256, 256, 128)

  # Test AXIAL
  expect_equal(get_slice_number(0, "AXIAL", dims), 1)
  expect_equal(get_slice_number(0.5, "AXIAL", dims), 64)
  expect_equal(get_slice_number(1, "AXIAL", dims), 128)

  # Test SAGITTAL
  expect_equal(get_slice_number(0, "SAGITTAL", dims), 1)
  expect_equal(get_slice_number(0.5, "SAGITTAL", dims), 128)
  expect_equal(get_slice_number(1, "SAGITTAL", dims), 256)
})

test_that("safe value extraction works", {
  test_data <- array(1:27, dim = c(3, 3, 3))

  # Valid coordinates
  expect_equal(safe_get_value(test_data, c(0, 0, 0), c(3, 3, 3)), 1)
  expect_equal(safe_get_value(test_data, c(2, 2, 2), c(3, 3, 3)), 27)

  # Out of bounds
  expect_true(is.na(safe_get_value(test_data, c(3, 0, 0), c(3, 3, 3))))
  expect_true(is.na(safe_get_value(test_data, c(-1, 0, 0), c(3, 3, 3))))
})

test_that("batch coordinate mapping works", {
  # Create test integrated data
  test_integrated <- data.frame(
    gaze_id = 1:10,
    time_sec = seq(0, 9),
    time_aligned = seq(0, 9),
    gaze_x = runif(10, 0.3, 0.7),
    gaze_y = runif(10, 0.3, 0.7),
    slice_index = rep(0.5, 10),
    plane = rep("AXIAL", 10),
    image_id = rep("test", 10),
    stringsAsFactors = FALSE
  )

  # Create test NIfTI data
  test_dims <- c(64, 64, 32)
  test_nifti_data <- list(
    nvimage = NVImage$new(
      matRAS = diag(4),
      dimsRAS = c(3, test_dims[1], test_dims[2], test_dims[3]),
      pixDimsRAS = c(1, 1, 1, 1)
    ),
    data = array(rnorm(prod(test_dims)), dim = test_dims),
    dims = test_dims
  )
  test_nifti_data$nvimage$calculateOblique()

  # Test mapping
  locations <- map_gaze_to_anatomy(test_integrated, test_nifti_data)

  expect_equal(nrow(locations), 10)
  expect_true(all(c("mm_x", "mm_y", "mm_z", "vox_x", "vox_y", "vox_z") %in% names(locations)))

  # Test summarization
  summary <- batch_locate_gaze(test_integrated, test_nifti_data, summarize = TRUE)
  expect_true(nrow(summary) <= nrow(locations))
})

test_that("resolve_coordinates works with display frame", {
  # Test center point
  coords <- resolve_coordinates(0.5, 0.5)
  expect_true(coords$in_bounds)
  # Center of frame should map near center of canvas
  expect_gt(coords$x, 0.4)
  expect_lt(coords$x, 0.6)
  expect_gt(coords$y, 0.4)
  expect_lt(coords$y, 0.6)

  # Test top-left corner (should be out of bounds)
  coords_tl <- resolve_coordinates(0, 0)
  expect_false(coords_tl$in_bounds)

  # Test bottom-right corner
  coords_br <- resolve_coordinates(1, 1)
  expect_false(coords_br$in_bounds)  # Right border cuts off

  # Test exact canvas boundaries
  # Left edge of canvas: 350/2624 ≈ 0.1334
  left_edge <- 350 / 2624
  coords_left <- resolve_coordinates(left_edge, 0.5)
  expect_equal(coords_left$x, 0, tolerance = 0.001)

  # Right edge of canvas: (350+1924)/2624 ≈ 0.8668
  right_edge <- (350 + 1924) / 2624
  coords_right <- resolve_coordinates(right_edge, 0.5)
  expect_equal(coords_right$x, 1, tolerance = 0.001)

  # Top edge of canvas: 80/1640 ≈ 0.0488
  top_edge <- 80 / 1640
  coords_top <- resolve_coordinates(0.5, top_edge)
  expect_equal(coords_top$y, 0, tolerance = 0.001)

  # Bottom edge of canvas: canvas extends to bottom (y=1)
  coords_bottom <- resolve_coordinates(0.5, 1.0)
  expect_equal(coords_bottom$y, 1, tolerance = 0.001)
  expect_true(coords_bottom$in_bounds)

  # Test with default DPR
  coords_default <- resolve_coordinates(0.5, 0.5)
  expect_equal(coords_default$x, coords$x)  # Should be same with DPR 1.25
})
