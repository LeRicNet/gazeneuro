test_that("NVImage class works correctly", {
  # Test data
  test_matRAS <- matrix(c(
    0.4295729696750641, 1.4573031670295222e-10, 0.11543580144643784, -113.94583892822266,
    -0.0032469534780830145, 0.40601974725723267, 1.6360894441604614, -120.05899047851562,
    -0.009373842738568783, -0.14063891768455505, 4.723334312438965, 77.72416687011719,
    0, 0, 0, 1
  ), nrow = 4, ncol = 4, byrow = FALSE)

  test_dimsRAS <- c(3, 512, 512, 25)
  test_pixDimsRAS <- c(1, 0.4296875, 0.4296875, 5)

  # Create NVImage instance
  nvimg <- NVImage$new(
    matRAS = test_matRAS,
    dimsRAS = test_dimsRAS,
    pixDimsRAS = test_pixDimsRAS
  )

  # Test initialization
  expect_equal(nvimg$matRAS, test_matRAS)
  expect_equal(nvimg$dimsRAS, test_dimsRAS)
  expect_equal(nvimg$pixDimsRAS, test_pixDimsRAS)

  # Test calculateOblique
  expect_no_error(nvimg$calculateOblique())
  expect_false(is.null(nvimg$frac2mm))
  expect_false(is.null(nvimg$frac2mmOrtho))

  # Test matrix dimensions
  expect_equal(dim(nvimg$frac2mm), c(4, 4))
  expect_equal(dim(nvimg$frac2mmOrtho), c(4, 4))

  # Test mm2vox
  test_mm <- c(0, 0, 0)
  voxel <- nvimg$mm2vox(test_mm)
  expect_length(voxel, 3)
  expect_true(all(is.finite(voxel)))

  # Test convertFrac2MM
  test_frac <- c(0.5, 0.5, 0.5)
  mm_coords <- nvimg$convertFrac2MM(test_frac)
  expect_length(mm_coords, 4)
  expect_equal(mm_coords[4], 1)  # Homogeneous coordinate
})

test_that("GL-matrix functions work correctly", {
  # Test mat4_create
  identity <- mat4_create()
  expect_equal(dim(identity), c(4, 4))
  expect_equal(diag(identity), rep(1, 4))

  # Test mat4_clone
  mat <- matrix(1:16, 4, 4)
  cloned <- mat4_clone(mat)
  expect_equal(cloned, mat)
  expect_false(identical(cloned, mat))  # Different objects

  # Test mat4_transpose
  mat <- matrix(1:16, 4, 4)
  transposed <- mat4_transpose(mat)
  expect_equal(transposed, t(mat))

  # Test mat4_invert on identity
  identity <- mat4_create()
  inv_identity <- mat4_invert(identity)
  expect_equal(inv_identity, identity)

  # Test vec4_transformMat4
  vec <- c(1, 2, 3, 1)
  identity <- mat4_create()
  transformed <- vec4_transformMat4(vec, identity)
  expect_equal(transformed, vec)
})

test_that("Integration with invalid data fails gracefully", {
  # Test with missing columns
  bad_gaze <- data.frame(
    time = 1:10,
    x = runif(10)
  )

  bad_z <- data.frame(
    time = 1:5,
    slice = 1:5
  )

  expect_error(integrate_all_gaze_points(bad_gaze, bad_z))
})

test_that("Coordinate transformations are reversible", {
  # Create test NVImage
  test_mat <- diag(4)
  test_mat[1:3, 4] <- c(-90, -126, -72)

  nvimg <- NVImage$new(
    matRAS = test_mat,
    dimsRAS = c(3, 182, 218, 182),
    pixDimsRAS = c(1, 1, 1, 1)
  )
  nvimg$calculateOblique()

  # Test round-trip conversion
  original_mm <- c(10, 20, 30)
  voxel <- nvimg$mm2vox(original_mm, frac = TRUE)

  # Convert voxel to fractional
  frac <- voxel / nvimg$dimsRAS[2:4]

  # Convert back to mm
  recovered_mm <- nvimg$convertFrac2MM(frac, isForceSliceMM = TRUE)

  # Check if close (within numerical precision)
  expect_equal(recovered_mm[1:3], original_mm, tolerance = 1e-5)
})
