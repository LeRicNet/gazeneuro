# Standalone GL-Matrix Solution for NVImage
# This version doesn't require creating a package - just source this file

library(Rcpp)
library(R6)

# Define the C++ code inline
glmatrix_cpp <- '
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix mat4_create() {
  NumericMatrix out(4, 4);
  out(0,0) = 1.0; out(1,0) = 0.0; out(2,0) = 0.0; out(3,0) = 0.0;
  out(0,1) = 0.0; out(1,1) = 1.0; out(2,1) = 0.0; out(3,1) = 0.0;
  out(0,2) = 0.0; out(1,2) = 0.0; out(2,2) = 1.0; out(3,2) = 0.0;
  out(0,3) = 0.0; out(1,3) = 0.0; out(2,3) = 0.0; out(3,3) = 1.0;
  return out;
}

// [[Rcpp::export]]
NumericMatrix mat4_clone(NumericMatrix a) {
  return Rcpp::clone(a);
}

// [[Rcpp::export]]
NumericMatrix mat4_transpose(NumericMatrix a) {
  NumericMatrix out(4, 4);
  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 4; j++) {
      out(i,j) = a(j,i);
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix mat4_invert(NumericMatrix a) {
  NumericMatrix out(4, 4);
  
  double a00 = a(0,0), a01 = a(0,1), a02 = a(0,2), a03 = a(0,3);
  double a10 = a(1,0), a11 = a(1,1), a12 = a(1,2), a13 = a(1,3);
  double a20 = a(2,0), a21 = a(2,1), a22 = a(2,2), a23 = a(2,3);
  double a30 = a(3,0), a31 = a(3,1), a32 = a(3,2), a33 = a(3,3);
  
  double b00 = a00 * a11 - a01 * a10;
  double b01 = a00 * a12 - a02 * a10;
  double b02 = a00 * a13 - a03 * a10;
  double b03 = a01 * a12 - a02 * a11;
  double b04 = a01 * a13 - a03 * a11;
  double b05 = a02 * a13 - a03 * a12;
  double b06 = a20 * a31 - a21 * a30;
  double b07 = a20 * a32 - a22 * a30;
  double b08 = a20 * a33 - a23 * a30;
  double b09 = a21 * a32 - a22 * a31;
  double b10 = a21 * a33 - a23 * a31;
  double b11 = a22 * a33 - a23 * a32;
  
  double det = b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;
  
  if (std::abs(det) < 1e-10) {
    Rcpp::stop("Matrix is not invertible");
  }
  
  det = 1.0 / det;
  
  out(0,0) = (a11 * b11 - a12 * b10 + a13 * b09) * det;
  out(0,1) = (a02 * b10 - a01 * b11 - a03 * b09) * det;
  out(0,2) = (a31 * b05 - a32 * b04 + a33 * b03) * det;
  out(0,3) = (a22 * b04 - a21 * b05 - a23 * b03) * det;
  out(1,0) = (a12 * b08 - a10 * b11 - a13 * b07) * det;
  out(1,1) = (a00 * b11 - a02 * b08 + a03 * b07) * det;
  out(1,2) = (a32 * b02 - a30 * b05 - a33 * b01) * det;
  out(1,3) = (a20 * b05 - a22 * b02 + a23 * b01) * det;
  out(2,0) = (a10 * b10 - a11 * b08 + a13 * b06) * det;
  out(2,1) = (a01 * b08 - a00 * b10 - a03 * b06) * det;
  out(2,2) = (a30 * b04 - a31 * b02 + a33 * b00) * det;
  out(2,3) = (a21 * b02 - a20 * b04 - a23 * b00) * det;
  out(3,0) = (a11 * b07 - a10 * b09 - a12 * b06) * det;
  out(3,1) = (a00 * b09 - a01 * b07 + a02 * b06) * det;
  out(3,2) = (a31 * b01 - a30 * b03 - a32 * b00) * det;
  out(3,3) = (a20 * b03 - a21 * b01 + a22 * b00) * det;
  
  return out;
}

// [[Rcpp::export]]
NumericMatrix mat4_translate(NumericMatrix a, NumericVector v) {
  NumericMatrix out = Rcpp::clone(a);
  double x = v[0], y = v[1], z = v[2];
  
  out(0,3) = a(0,0) * x + a(0,1) * y + a(0,2) * z + a(0,3);
  out(1,3) = a(1,0) * x + a(1,1) * y + a(1,2) * z + a(1,3);
  out(2,3) = a(2,0) * x + a(2,1) * y + a(2,2) * z + a(2,3);
  out(3,3) = a(3,0) * x + a(3,1) * y + a(3,2) * z + a(3,3);
  
  return out;
}

// [[Rcpp::export]]
NumericMatrix mat4_scale_columns_directly(NumericMatrix mat, NumericVector dim) {
  NumericMatrix out = Rcpp::clone(mat);
  
  // Scale first column by dim[0]
  out(0,0) *= dim[0];
  out(1,0) *= dim[0];
  out(2,0) *= dim[0];
  
  // Scale second column by dim[1]
  out(0,1) *= dim[1];
  out(1,1) *= dim[1];
  out(2,1) *= dim[1];
  
  // Scale third column by dim[2]
  out(0,2) *= dim[2];
  out(1,2) *= dim[2];
  out(2,2) *= dim[2];
  
  return out;
}

// [[Rcpp::export]]
NumericVector vec4_transformMat4(NumericVector a, NumericMatrix m) {
  NumericVector out(4);
  double x = a[0], y = a[1], z = a[2], w = a[3];
  
  out[0] = m(0,0) * x + m(0,1) * y + m(0,2) * z + m(0,3) * w;
  out[1] = m(1,0) * x + m(1,1) * y + m(1,2) * z + m(1,3) * w;
  out[2] = m(2,0) * x + m(2,1) * y + m(2,2) * z + m(2,3) * w;
  out[3] = m(3,0) * x + m(3,1) * y + m(3,2) * z + m(3,3) * w;
  
  return out;
}

// [[Rcpp::export]]
void print_glmatrix(NumericMatrix mat) {
  Rcout << "gl-matrix format (column-major):" << std::endl;
  for(int j = 0; j < 4; j++) {
    Rcout << "[";
    for(int i = 0; i < 4; i++) {
      Rcout << mat(i,j);
      if(i < 3) Rcout << ", ";
    }
    Rcout << "]";
    if(j < 3) Rcout << ",";
    Rcout << std::endl;
  }
}
'

# Compile the C++ code
sourceCpp(code = glmatrix_cpp)

# Now define the NVImage class
NVImage <- R6Class(
  "NVImage",
  public = list(
    matRAS = NULL,
    dimsRAS = NULL,
    pixDimsRAS = NULL,
    frac2mm = NULL,
    frac2mmOrtho = NULL,
    extentsMinOrtho = NULL,
    extentsMaxOrtho = NULL,
    
    initialize = function(matRAS = NULL, dimsRAS = NULL, pixDimsRAS = NULL) {
      self$matRAS <- matRAS
      self$dimsRAS <- dimsRAS
      self$pixDimsRAS <- pixDimsRAS
    },
    
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
    
    calculateOblique = function() {
      if (is.null(self$matRAS) || is.null(self$dimsRAS) || is.null(self$pixDimsRAS)) {
        stop("matRAS, dimsRAS, and pixDimsRAS must be set")
      }
      
      dim <- c(self$dimsRAS[2], self$dimsRAS[3], self$dimsRAS[4], 1)
      
      sform <- mat4_clone(self$matRAS)
      sform <- mat4_transpose(sform)
      
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
    }
  )
)

# Helper function to debug matrix differences
debug_matrix_differences <- function(mat1, mat2, name = "Matrix") {
  cat("\n", name, "differences:\n")
  diff_mat <- abs(mat1 - mat2)
  max_diff_idx <- which(diff_mat == max(diff_mat), arr.ind = TRUE)[1,]
  
  cat("  Maximum difference:", max(diff_mat), "at position [", 
      max_diff_idx[1], ",", max_diff_idx[2], "]\n")
  cat("  R value:", mat1[max_diff_idx[1], max_diff_idx[2]], "\n")
  cat("  JS value:", mat2[max_diff_idx[1], max_diff_idx[2]], "\n")
  
  # Show all differences > 1e-7
  significant_diffs <- which(diff_mat > 1e-7, arr.ind = TRUE)
  if (nrow(significant_diffs) > 0 && nrow(significant_diffs) <= 5) {
    cat("  All differences > 1e-7:\n")
    for (i in 1:nrow(significant_diffs)) {
      r <- significant_diffs[i, 1]
      c <- significant_diffs[i, 2]
      cat("    [", r, ",", c, "]: R=", mat1[r,c], ", JS=", mat2[r,c], 
          ", diff=", diff_mat[r,c], "\n")
    }
  }
}

# Test function
test_nvimage_glmatrix <- function() {
  cat("=== Testing NVImage with GL-Matrix Implementation ===\n\n")
  
  # Real test data
  real_matRAS <- matrix(c(
    0.4295729696750641, 1.4573031670295222e-10, 0.11543580144643784, -113.94583892822266,
    -0.0032469534780830145, 0.40601974725723267, 1.6360894441604614, -120.05899047851562,
    -0.009373842738568783, -0.14063891768455505, 4.723334312438965, 77.72416687011719,
    0, 0, 0, 1
  ), nrow = 4, ncol = 4, byrow = FALSE)
  
  real_dimsRAS <- c(3, 512, 512, 25)
  real_pixDimsRAS <- c(1, 0.4296875, 0.4296875, 5)
  
  # Expected results
  expected_frac2mm <- matrix(c(
    219.9413604736328, -1.6624401807785034, -4.799407482147217, 0,
    7.461392215191154e-8, 207.88211059570312, -72.00712585449219, 0,
    2.885895013809204, 40.90223693847656, 118.08335876464844, 0,
    -114.21834564208984, -121.07842254638672, 75.43750762939453, 1
  ), nrow = 4, ncol = 4, byrow = FALSE)
  
  expected_frac2mmOrtho <- matrix(c(
    220, 0, 0, 0,
    0, 220, 0, 0,
    0, 0, 125, 0,
    -114.91867065429688, -139.10035705566406, 29.007308959960938, 1
  ), nrow = 4, ncol = 4, byrow = FALSE)
  
  # Create instance
  nvimg <- NVImage$new(
    matRAS = real_matRAS,
    dimsRAS = real_dimsRAS,
    pixDimsRAS = real_pixDimsRAS
  )
  
  # Calculate matrices
  nvimg$calculateOblique()
  
  # Compare results
  cat("\nfrac2mm comparison:\n")
  cat("Max difference:", max(abs(nvimg$frac2mm - expected_frac2mm)), "\n")
  
  cat("\nfrac2mmOrtho comparison:\n")
  cat("Max difference:", max(abs(nvimg$frac2mmOrtho - expected_frac2mmOrtho)), "\n")
  
  # Detailed debugging if needed
  if (max(abs(nvimg$frac2mm - expected_frac2mm)) > 1e-7) {
    debug_matrix_differences(nvimg$frac2mm, expected_frac2mm, "frac2mm")
  }
  
  if (max(abs(nvimg$frac2mmOrtho - expected_frac2mmOrtho)) > 1e-7) {
    debug_matrix_differences(nvimg$frac2mmOrtho, expected_frac2mmOrtho, "frac2mmOrtho")
  }
  
  # Test result
  tolerance <- 1e-5  # Reasonable tolerance for floating-point comparison
  frac2mm_diff <- max(abs(nvimg$frac2mm - expected_frac2mm))
  frac2mmOrtho_diff <- max(abs(nvimg$frac2mmOrtho - expected_frac2mmOrtho))
  
  cat("\nDetailed comparison:\n")
  cat("frac2mm max difference:", frac2mm_diff, 
      ifelse(frac2mm_diff < tolerance, "✓ PASS", "✗ FAIL"), "\n")
  cat("frac2mmOrtho max difference:", frac2mmOrtho_diff, 
      ifelse(frac2mmOrtho_diff < tolerance, "✓ PASS", "✗ FAIL"), "\n")
  
  if (frac2mm_diff < tolerance && frac2mmOrtho_diff < tolerance) {
    cat("\n✅ SUCCESS! GL-Matrix implementation matches JavaScript within floating-point precision!\n")
    cat("Note: Differences of", max(frac2mm_diff, frac2mmOrtho_diff), 
        "are negligible for medical imaging applications.\n")
    return(TRUE)
  } else {
    cat("\n❌ FAILED! Results don't match within tolerance.\n")
    return(FALSE)
  }
}

# Run the test
test_nvimage_glmatrix()

# Example usage
cat("\n=== Example Usage ===\n")
# Your actual NIfTI data would go here
example_affine <- diag(4)
example_affine[1:3, 4] <- c(-90, -126, -72)  # Set origin

nvimg <- NVImage$new(
  matRAS = example_affine,
  dimsRAS = c(3, 182, 218, 182),
  pixDimsRAS = c(1, 1, 1, 1)
)

nvimg$calculateOblique()

# Convert a voxel coordinate to world coordinate
voxel <- c(91, 109, 91)  # Center voxel
frac <- voxel / nvimg$dimsRAS[2:4]
world <- nvimg$convertFrac2MM(frac)
cat("Voxel", voxel, "-> World", world[1:3], "\n")