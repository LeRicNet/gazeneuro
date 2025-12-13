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

test_that("Coordinate transformations are consistent with voxel center convention", {
  # Create test NVImage with simple transformation
  test_mat <- diag(4)
  test_mat[1:3, 4] <- c(-90, -126, -72)  # Translation only

  nvimg <- NVImage$new(
    matRAS = test_mat,
    dimsRAS = c(3, 182, 218, 182),
    pixDimsRAS = c(1, 1, 1, 1)
  )
  nvimg$calculateOblique()

  # Test that frac2mm and mm2frac are inverses
  # The frac2mm matrix includes a -0.5 offset for voxel center convention
  # To test round-trip, we must use the same transformation in both directions

  # Forward: mm -> frac using inverse of frac2mm
  original_mm <- c(10, 20, 30)
  frac <- nvimg$convertMM2Frac(original_mm)

  # Backward: frac -> mm using frac2mm
  recovered_mm <- nvimg$convertFrac2MM(frac, isForceSliceMM = TRUE)

  # Check if close (within numerical precision)
  expect_equal(recovered_mm[1:3], original_mm, tolerance = 1e-5)
})

test_that("vox2mm and mm2vox are inverses", {
  # Create test NVImage
  test_mat <- diag(4)
  test_mat[1:3, 4] <- c(-90, -126, -72)

  nvimg <- NVImage$new(
    matRAS = test_mat,
    dimsRAS = c(3, 182, 218, 182),
    pixDimsRAS = c(1, 1, 1, 1)
  )

  # Test round-trip: mm -> vox -> mm (using vox2mm which applies matRAS directly)
  original_mm <- c(10, 20, 30)
  voxel <- nvimg$mm2vox(original_mm, frac = TRUE)
  recovered_mm <- nvimg$vox2mm(voxel)

  expect_equal(recovered_mm, original_mm, tolerance = 1e-10)
})

test_that("Fractional coordinates are bounded correctly", {
  test_mat <- diag(4)
  test_mat[1:3, 4] <- c(-90, -126, -72)

  nvimg <- NVImage$new(
    matRAS = test_mat,
    dimsRAS = c(3, 182, 218, 182),
    pixDimsRAS = c(1, 1, 1, 1)
  )
  nvimg$calculateOblique()

  # Center of volume should be frac = (0.5, 0.5, 0.5)
  center_frac <- c(0.5, 0.5, 0.5)
  center_mm <- nvimg$convertFrac2MM(center_frac, isForceSliceMM = TRUE)

  # For identity rotation with translation (-90, -126, -72):
  # frac=0.5 should map to voxel center of (91, 109, 91)
  # Then mm = voxel + origin = (91-90, 109-126, 91-72) = (1, -17, 19)
  # But with -0.5 shim: voxel_effective = (91-0.5, 109-0.5, 91-0.5) = (90.5, 108.5, 90.5)
  # mm = (90.5-90, 108.5-126, 90.5-72) = (0.5, -17.5, 18.5)

  # The center of the image in mm coordinates
  expected_center_mm <- c(0.5, -17.5, 18.5)
  expect_equal(center_mm[1:3], expected_center_mm, tolerance = 1e-5)
})

test_that("Edge voxels map correctly", {
  test_mat <- diag(4)
  test_mat[1:3, 4] <- c(-90, -126, -72)

  nvimg <- NVImage$new(
    matRAS = test_mat,
    dimsRAS = c(3, 182, 218, 182),
    pixDimsRAS = c(1, 1, 1, 1)
  )
  nvimg$calculateOblique()

  # frac = (0, 0, 0) should map to center of first voxel
  # With -0.5 shim: voxel = (-0.5, -0.5, -0.5)
  # mm = voxel + origin = (-0.5-90, -0.5-126, -0.5-72) = (-90.5, -126.5, -72.5)
  first_voxel_mm <- nvimg$convertFrac2MM(c(0, 0, 0), isForceSliceMM = TRUE)
  expected_first <- c(-90.5, -126.5, -72.5)
  expect_equal(first_voxel_mm[1:3], expected_first, tolerance = 1e-5)

  # frac = (1, 1, 1) should map to center of last voxel
  # With dims (182, 218, 182) and -0.5 shim:
  # voxel = (182-0.5, 218-0.5, 182-0.5) = (181.5, 217.5, 181.5)
  # mm = voxel + origin = (181.5-90, 217.5-126, 181.5-72) = (91.5, 91.5, 109.5)
  last_voxel_mm <- nvimg$convertFrac2MM(c(1, 1, 1), isForceSliceMM = TRUE)
  expected_last <- c(91.5, 91.5, 109.5)
  expect_equal(last_voxel_mm[1:3], expected_last, tolerance = 1e-5)
})

test_that("Pixel dimensions affect coordinate transformations", {
  test_mat <- diag(4)
  test_mat[1:3, 4] <- c(0, 0, 0)  # No translation for simplicity

  # Create with anisotropic voxels
  nvimg <- NVImage$new(
    matRAS = test_mat,
    dimsRAS = c(3, 100, 100, 50),  # Different z dimension
    pixDimsRAS = c(1, 1, 1, 2)     # z voxels are 2mm thick
  )
  nvimg$calculateOblique()

  # Check that extents are calculated correctly
  # X extent: 100 * 1mm = 100mm
  # Y extent: 100 * 1mm = 100mm
  # Z extent: 50 * 2mm = 100mm

  x_extent <- nvimg$extentsMaxOrtho[1] - nvimg$extentsMinOrtho[1]
  y_extent <- nvimg$extentsMaxOrtho[2] - nvimg$extentsMinOrtho[2]
  z_extent <- nvimg$extentsMaxOrtho[3] - nvimg$extentsMinOrtho[3]

  expect_equal(x_extent, 100, tolerance = 1e-5)
  expect_equal(y_extent, 100, tolerance = 1e-5)
  expect_equal(z_extent, 100, tolerance = 1e-5)
})
