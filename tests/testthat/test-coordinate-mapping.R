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

test_that("resolve_coordinates works", {
  coords <- resolve_coordinates(100, 200, dpr = 2)
  expect_equal(coords$x, 200)
  expect_equal(coords$y, 400)

  coords_default <- resolve_coordinates(100, 200)
  expect_equal(coords_default$x, 100)
  expect_equal(coords_default$y, 200)
})
