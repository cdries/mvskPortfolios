/*
 * Functions copied from PerformanceAnalytics in order to pass CRAN checks
 */

#include <R.h>
#include <Rinternals.h>


SEXP  M3mat2vec(SEXP XX, SEXP PP){
  /*
   arguments
   XX        : numeric vector with full vectorized coskewness matrix
   PP        : integer, number of assets

   Written by Dries Cornilly
   */

  // // declare pointers for the vectors
  double *X;

  // use REAL() to access the C array inside the numeric vector passed in from R
  X = REAL(XX);
  int P = asInteger(PP);

  // allocate and construct the vector with unique coskewness elements
  SEXP M3vec = PROTECT(allocVector(REALSXP, P * (P + 1) * (P + 2) / 6));
  double *rM3vec = REAL(M3vec);

  int iter = 0;
  // construct the vector with unique coskewness elements
  for (int ii = 0; ii < P; ii++) {
    for (int jj = ii; jj < P; jj++) {
      for (int kk = jj; kk < P; kk++) {
        // extract unique element
        rM3vec[iter] = X[(ii * P + jj) * P + kk];
        iter++;
      } // loop kk
    } // loop jj
  } // loop ii

  UNPROTECT(1);
  return M3vec;
}

SEXP  M4mat2vec(SEXP XX, SEXP PP){
  /*
   arguments
   XX        : numeric vector with full vectorized cokurtosis matrix
   PP        : integer, number of assets

   Written by Dries Cornilly
   */

  // // declare pointers for the vectors
  double *X;

  // use REAL() to access the C array inside the numeric vector passed in from R
  X = REAL(XX);
  int P = asInteger(PP);

  // allocate and construct the vector with unique cokurtosis elements
  SEXP M4vec = PROTECT(allocVector(REALSXP, P * (P + 1) * (P + 2) * (P + 3) / 24));
  double *rM4vec = REAL(M4vec);

  int iter = 0;
  // construct the vector with unique cokurtosis elements
  for (int ii = 0; ii < P; ii++) {
    for (int jj = ii; jj < P; jj++) {
      for (int kk = jj; kk < P; kk++) {
        for (int ll = kk; ll < P; ll++) {
          // extract unique element
          rM4vec[iter] = X[(ii * P * P + jj * P + kk) * P + ll];
          iter++;
        } // loop ll
      } // loop kk
    } // loop jj
  } // loop ii

  UNPROTECT(1);
  return M4vec;
}


// // //
// // Compute portfolio skewness and kurtosis and the derivatives with respect to w
// // //

SEXP  M3port(SEXP WW, SEXP XX, SEXP PP){
  /*
   arguments
   WW        : numeric vector with the portfolio weights
   XX        : numeric vector with unique elements of a coskewness matrix
   PP        : integer, number of assets

   Written by Dries Cornilly
   */

  // // declare pointers for the vectors
  double *X, *W;

  // use REAL() to access the C array inside the numeric vector passed in from R
  X = REAL(XX);
  W = REAL(WW);
  int P = asInteger(PP);

  // allocate and compute the portfolio skewness
  SEXP skew = PROTECT(allocVector(REALSXP, 1));
  double *rskew = REAL(skew);
  rskew[0] = 0.0;

  int iter = 0;
  // compute the portfolio skewness
  for (int ii = 0; ii < P; ii++) {
    for (int jj = ii; jj < P; jj++) {
      for (int kk = jj; kk < P; kk++) {
        // add to the matrix product
        if (ii == jj) {
          if (jj == kk) {
            // element phi_iii
            rskew[0] += X[iter] * W[ii] * W[ii] * W[ii];
          } else {
            // element phi_iik
            rskew[0] += 3.0 * X[iter] * W[ii] * W[ii] * W[kk];
          }
        } else {
          if (jj == kk) {
            // element phi_ijj
            rskew[0] += 3.0 * X[iter] * W[ii] * W[jj] * W[jj];
          } else {
            // element phi_ijk
            rskew[0] += 6.0 * X[iter] * W[ii] * W[jj] * W[kk];
          }
        }
        iter++;
      } // loop kk
    } // loop jj
  } // loop ii

  UNPROTECT(1);
  return skew;
}

SEXP  M3port_grad(SEXP WW, SEXP XX, SEXP PP){
  /*
   arguments
   WW        : numeric vector with the portfolio weights
   XX        : numeric vector with unique elements of a coskewness matrix
   PP        : integer, number of assets

   Written by Dries Cornilly
   */

  // // declare pointers for the vectors
  double *X, *W;

  // use REAL() to access the C array inside the numeric vector passed in from R
  X = REAL(XX);
  W = REAL(WW);
  int P = asInteger(PP);

  // allocate and compute gradient of the portfolio skewness
  SEXP skew_grad = PROTECT(allocVector(REALSXP, P));
  double *rskew_grad = REAL(skew_grad);
  for (int ii = 0; ii < P; ii++) {
    rskew_grad[ii] = 0.0;
  }

  int iter = 0;
  // compute the gradient of the portfolio skewness
  for (int ii = 0; ii < P; ii++) {
    for (int jj = ii; jj < P; jj++) {
      for (int kk = jj; kk < P; kk++) {
        // add to the inner product
        if (ii == jj) {
          if (jj == kk) {
            // element phi_iii
            rskew_grad[ii] += 3.0 * X[iter] * W[ii] * W[ii];
          } else {
            // element phi_iik
            rskew_grad[ii] += 6.0 * X[iter] * W[ii] * W[kk];
            rskew_grad[kk] += 3.0 * X[iter] * W[ii] * W[ii];
          }
        } else {
          if (jj == kk) {
            // element phi_ijj
            rskew_grad[ii] += 3.0 * X[iter] * W[jj] * W[jj];
            rskew_grad[jj] += 6.0 * X[iter] * W[ii] * W[jj];
          } else {
            // element phi_ijk
            rskew_grad[ii] += 6.0 * X[iter] * W[jj] * W[kk];
            rskew_grad[jj] += 6.0 * X[iter] * W[ii] * W[kk];
            rskew_grad[kk] += 6.0 * X[iter] * W[ii] * W[jj];
          }
        }
        iter++;
      } // loop kk
    } // loop jj
  } // loop ii

  UNPROTECT(1);
  return skew_grad;
}

SEXP  M4port(SEXP WW, SEXP XX, SEXP PP){
  /*
   arguments
   WW        : numeric vector with the portfolio weights
   XX        : numeric vector with unique elements of a cokurtosis matrix
   PP        : integer, number of assets

   Written by Dries Cornilly
   */

  // // declare pointers for the vectors
  double *X, *W;

  // use REAL() to access the C array inside the numeric vector passed in from R
  X = REAL(XX);
  W = REAL(WW);
  int P = asInteger(PP);

  // allocate and compute the portfolio kurtosis
  SEXP kurt = PROTECT(allocVector(REALSXP, 1));
  double *rkurt = REAL(kurt);
  rkurt[0] = 0.0;

  int iter = 0;
  // compute the portfolio central fourth moment
  for (int ii = 0; ii < P; ii++) {
    for (int jj = ii; jj < P; jj++) {
      for (int kk = jj; kk < P; kk++) {
        for (int ll = kk; ll < P; ll++) {
          // add to the matrix product
          if (ii == jj) {
            if (jj == kk) {
              if (kk == ll) {
                // element psi_iiii
                rkurt[0] += X[iter] * W[ii] * W[ii] * W[ii] * W[ii];
              } else {
                // element psi_iiil
                rkurt[0] += 4.0 * X[iter] * W[ii] * W[ii] * W[ii] * W[ll];
              }
            } else {
              if (kk == ll) {
                // element psi_iikk
                rkurt[0] += 6.0 * X[iter] * W[ii] * W[ii] * W[kk] * W[kk];
              } else {
                // element psi_iikl
                rkurt[0] += 12.0 * X[iter] * W[ii] * W[ii] * W[kk] * W[ll];
              }
            }
          } else {
            if (jj == kk) {
              if (kk == ll) {
                // element psi_ijjj
                rkurt[0] += 4.0 * X[iter] * W[ii] * W[jj] * W[jj] * W[jj];
              } else {
                // element psi_ijjl
                rkurt[0] += 12.0 * X[iter] * W[ii] * W[jj] * W[jj] * W[ll];
              }
            } else {
              if (kk == ll) {
                // element psi_ijkk
                rkurt[0] += 12.0 * X[iter] * W[ii] * W[jj] * W[kk] * W[kk];
              } else {
                // element psi_ijkl
                rkurt[0] += 24.0 * X[iter] * W[ii] * W[jj] * W[kk] * W[ll];
              }
            }
          }
          iter++;
        } // loop ll
      } // loop kk
    } // loop jj
  } // loop ii

  UNPROTECT(1);
  return kurt;
}

SEXP  M4port_grad(SEXP WW, SEXP XX, SEXP PP){
  /*
   arguments
   WW        : numeric vector with the portfolio weights
   XX        : numeric vector with unique elements of a cokurtosis matrix
   PP        : integer, number of assets

   Written by Dries Cornilly
   */

  // // declare pointers for the vectors
  double *X, *W;

  // use REAL() to access the C array inside the numeric vector passed in from R
  X = REAL(XX);
  W = REAL(WW);
  int P = asInteger(PP);

  // allocate and compute the gradient of the portfolio kurtosis
  SEXP kurt_grad = PROTECT(allocVector(REALSXP, P));
  double *rkurt_grad = REAL(kurt_grad);
  for (int ii = 0; ii < P; ii++) {
    rkurt_grad[ii] = 0.0;
  }

  int iter = 0;
  // compute the gradient portfolio central fourth moment
  for (int ii = 0; ii < P; ii++) {
    for (int jj = ii; jj < P; jj++) {
      for (int kk = jj; kk < P; kk++) {
        for (int ll = kk; ll < P; ll++) {
          // add to the matrix product
          if (ii == jj) {
            if (jj == kk) {
              if (kk == ll) {
                // element psi_iiii
                rkurt_grad[ii] += 4.0 * X[iter] * W[ii] * W[ii] * W[ii];
              } else {
                // element psi_iiil
                rkurt_grad[ii] += 12.0 * X[iter] * W[ii] * W[ii] * W[ll];
                rkurt_grad[ll] += 4.0 * X[iter] * W[ii] * W[ii] * W[ii];
              }
            } else {
              if (kk == ll) {
                // element psi_iikk
                rkurt_grad[ii] += 12.0 * X[iter] * W[ii] * W[kk] * W[kk];
                rkurt_grad[kk] += 12.0 * X[iter] * W[ii] * W[ii] * W[kk];
              } else {
                // element psi_iikl
                rkurt_grad[ii] += 24.0 * X[iter] * W[ii] * W[kk] * W[ll];
                rkurt_grad[kk] += 12.0 * X[iter] * W[ii] * W[ii] * W[ll];
                rkurt_grad[ll] += 12.0 * X[iter] * W[ii] * W[ii] * W[kk];
              }
            }
          } else {
            if (jj == kk) {
              if (kk == ll) {
                // element psi_ijjj
                rkurt_grad[ii] += 4.0 * X[iter] * W[jj] * W[jj] * W[jj];
                rkurt_grad[jj] += 12.0 * X[iter] * W[ii] * W[jj] * W[jj];
              } else {
                // element psi_ijjl
                rkurt_grad[ii] += 12.0 * X[iter] * W[jj] * W[jj] * W[ll];
                rkurt_grad[jj] += 24.0 * X[iter] * W[ii] * W[jj] * W[ll];
                rkurt_grad[ll] += 12.0 * X[iter] * W[ii] * W[jj] * W[jj];
              }
            } else {
              if (kk == ll) {
                // element psi_ijkk
                rkurt_grad[ii] += 12.0 * X[iter] * W[jj] * W[kk] * W[kk];
                rkurt_grad[jj] += 12.0 * X[iter] * W[ii] * W[kk] * W[kk];
                rkurt_grad[kk] += 24.0 * X[iter] * W[ii] * W[jj] * W[kk];
              } else {
                // element psi_ijkl
                rkurt_grad[ii] += 24.0 * X[iter] * W[jj] * W[kk] * W[ll];
                rkurt_grad[jj] += 24.0 * X[iter] * W[ii] * W[kk] * W[ll];
                rkurt_grad[kk] += 24.0 * X[iter] * W[ii] * W[jj] * W[ll];
                rkurt_grad[ll] += 24.0 * X[iter] * W[ii] * W[jj] * W[kk];
              }
            }
          }
          iter++;
        } // loop ll
      } // loop kk
    } // loop jj
  } // loop ii

  UNPROTECT(1);
  return kurt_grad;
}
