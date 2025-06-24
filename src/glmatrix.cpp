#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' Create a 4x4 identity matrix
 //' @return A 4x4 identity matrix
 //' @export
 // [[Rcpp::export]]
 NumericMatrix mat4_create() {
   NumericMatrix out(4, 4);
   out(0,0) = 1.0; out(1,0) = 0.0; out(2,0) = 0.0; out(3,0) = 0.0;
   out(0,1) = 0.0; out(1,1) = 1.0; out(2,1) = 0.0; out(3,1) = 0.0;
   out(0,2) = 0.0; out(1,2) = 0.0; out(2,2) = 1.0; out(3,2) = 0.0;
   out(0,3) = 0.0; out(1,3) = 0.0; out(2,3) = 0.0; out(3,3) = 1.0;
   return out;
 }

 //' Clone a matrix
 //' @param a Matrix to clone
 //' @return A copy of the input matrix
 //' @export
 // [[Rcpp::export]]
 NumericMatrix mat4_clone(NumericMatrix a) {
   return Rcpp::clone(a);
 }

 //' Transpose a 4x4 matrix
 //' @param a Matrix to transpose
 //' @return Transposed matrix
 //' @export
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

 //' Invert a 4x4 matrix
 //' @param a Matrix to invert
 //' @return Inverted matrix
 //' @export
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

 //' Translate a matrix by a vector
 //' @param a Matrix to translate
 //' @param v Translation vector (x, y, z)
 //' @return Translated matrix
 //' @export
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

 //' Scale matrix columns directly
 //' @param mat Matrix to scale
 //' @param dim Scaling factors for each column
 //' @return Scaled matrix
 //' @export
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

 //' Transform a 4D vector by a 4x4 matrix
 //' @param a Vector to transform (x, y, z, w)
 //' @param m Transformation matrix
 //' @return Transformed vector
 //' @export
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

 //' Print matrix in gl-matrix format
 //' @param mat Matrix to print
 //' @export
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
