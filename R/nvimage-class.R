#' NVImage Class for NIfTI Image Processing
#'
#' @description
#' R6 class for handling NIfTI image transformations and coordinate conversions.
#' This class provides methods for converting between voxel and world coordinates
#' using gl-matrix style transformations.
#'
#' @export
#' @importFrom R6 R6Class
NVImage <- R6::R6Class(
  "NVImage",
  public = list(
    #' @field matRAS RAS transformation matrix
    matRAS = NULL,

    #' @field dimsRAS Dimensions in RAS space
    dimsRAS = NULL,

    #' @field pixDimsRAS Pixel dimensions in RAS space
    pixDimsRAS = NULL,

    #' @field frac2mm Fractional to mm transformation matrix
    frac2mm = NULL,

    #' @field frac2mmOrtho Orthographic fractional to mm transformation matrix
    frac2mmOrtho = NULL,

    #' @field extentsMinOrtho Minimum extents in orthographic space
    extentsMinOrtho = NULL,

    #' @field extentsMaxOrtho Maximum extents in orthographic space
    extentsMaxOrtho = NULL,

    #' @description
    #' Create a new NVImage object
    #' @param matRAS RAS transformation matrix (4x4)
    #' @param dimsRAS Dimensions vector (length 4: ndims, x, y, z)
    #' @param pixDimsRAS Pixel dimensions vector (length 4: 1, dx, dy, dz)
    #' @return A new NVImage object
    initialize = function(matRAS = NULL, dimsRAS = NULL, pixDimsRAS = NULL) {
      self$matRAS <- matRAS
      self$dimsRAS <- dimsRAS
      self$pixDimsRAS <- pixDimsRAS
    },

    #' @description
    #' Convert mm coordinates to voxel coordinates
    #' @param mm Numeric vector of mm coordinates (x, y, z)
    #' @param frac Logical, whether to return fractional voxel coordinates
    #' @return Numeric vector of voxel coordinates
    mm2vox = function(mm, frac = FALSE) {
      if (is.null(self$matRAS)) {
        stop("matRAS undefined")
      }

      sform <- mat4_clone(self$matRAS)
      out <- mat4_transpose(sform)
      out <- mat4_invert(out)

      pos <- c(mm[1], mm[2], mm[3], 1)
      pos <- vec4_transformMat4(pos, out)
      pos3 <- pos[1:3]

      if (frac) {
        return(pos3)
      }
      return(round(pos3))
    },

    #' @description
    #' Convert voxel coordinates to mm coordinates
    #' @param vox Numeric vector of voxel coordinates (x, y, z)
    #' @return Numeric vector of mm coordinates
    vox2mm = function(vox) {
      if (is.null(self$matRAS)) {
        stop("matRAS undefined")
      }

      sform <- mat4_clone(self$matRAS)
      sform <- mat4_transpose(sform)

      pos <- c(vox[1], vox[2], vox[3], 1)
      pos <- vec4_transformMat4(pos, sform)

      return(pos[1:3])
    },

    #' @description
    #' Calculate oblique transformation matrices
    #' Sets up frac2mm and frac2mmOrtho matrices for coordinate transformations
    calculateOblique = function() {
      if (is.null(self$matRAS) || is.null(self$dimsRAS) || is.null(self$pixDimsRAS)) {
        stop("matRAS, dimsRAS, and pixDimsRAS must be set")
      }

      dim <- c(self$dimsRAS[2], self$dimsRAS[3], self$dimsRAS[4], 1)

      sform <- mat4_clone(self$matRAS)
      sform <- mat4_transpose(sform)

      # The -0.5 shim accounts for voxel center convention:
      # frac=0 maps to center of first voxel, not its corner
      shim <- c(-0.5, -0.5, -0.5)
      sform <- mat4_translate(sform, shim)

      sform <- mat4_scale_columns_directly(sform, dim[1:3])

      self$frac2mm <- sform

      # Orthographic version
      pixdimX <- self$pixDimsRAS[2]
      pixdimY <- self$pixDimsRAS[3]
      pixdimZ <- self$pixDimsRAS[4]

      oform <- mat4_clone(sform)

      # Set orthographic elements
      oform[1, 1] <- pixdimX * dim[1]
      oform[2, 1] <- 0
      oform[3, 1] <- 0
      oform[1, 2] <- 0
      oform[2, 2] <- pixdimY * dim[2]
      oform[3, 2] <- 0
      oform[1, 3] <- 0
      oform[2, 3] <- 0
      oform[3, 3] <- pixdimZ * dim[3]

      originVoxel <- self$mm2vox(c(0, 0, 0), frac = TRUE)

      oform[1, 4] <- (-originVoxel[1] - 0.5) * pixdimX
      oform[2, 4] <- (-originVoxel[2] - 0.5) * pixdimY
      oform[3, 4] <- (-originVoxel[3] - 0.5) * pixdimZ

      self$frac2mmOrtho <- oform

      self$extentsMinOrtho <- c(oform[1, 4], oform[2, 4], oform[3, 4])
      self$extentsMaxOrtho <- c(
        oform[1, 1] + oform[1, 4],
        oform[2, 2] + oform[2, 4],
        oform[3, 3] + oform[3, 4]
      )
    },

    #' @description
    #' Convert fractional coordinates to mm
    #' @param frac Numeric vector of fractional coordinates (0-1)
    #' @param isForceSliceMM Logical, whether to use slice transformation
    #' @return Numeric vector of mm coordinates
    convertFrac2MM = function(frac, isForceSliceMM = FALSE) {
      if (is.null(self$frac2mm) || is.null(self$frac2mmOrtho)) {
        stop("Must call calculateOblique() first")
      }

      pos <- c(frac[1], frac[2], frac[3], 1)

      if (isForceSliceMM) {
        transformation_matrix <- self$frac2mm
      } else {
        transformation_matrix <- self$frac2mmOrtho
      }

      pos_transformed <- vec4_transformMat4(pos, transformation_matrix)
      return(pos_transformed)
    },

    #' @description
    #' Convert mm coordinates to fractional coordinates
    #' This is the inverse of convertFrac2MM with isForceSliceMM=TRUE
    #' @param mm Numeric vector of mm coordinates (x, y, z)
    #' @return Numeric vector of fractional coordinates (0-1)
    convertMM2Frac = function(mm) {
      if (is.null(self$frac2mm)) {
        stop("Must call calculateOblique() first")
      }

      # Invert frac2mm to get mm2frac
      mm2frac <- mat4_invert(self$frac2mm)

      pos <- c(mm[1], mm[2], mm[3], 1)
      pos_transformed <- vec4_transformMat4(pos, mm2frac)

      return(pos_transformed[1:3])
    }
  )
)
