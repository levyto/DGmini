// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Definitions of linear algebra operations (BLAS-like)
//
// -----------------------------------------------------------------------------

#ifndef LINALG_H
#define LINALG_H

#include "Algebra/Vec.h"
#include "Algebra/Mat.h"

double nrm2(const Vec& x);
double dot(const Vec& x, const Vec& y);
void scal(double alpha, Vec& x);
void axpy(double alpha, const Vec& x, Vec& y);
void gemv(const Mat& A, const Vec& x, Vec& y);
Vec  gemv(const Mat& A, const Vec& x);
void gemm(const Mat& A, const Mat& B, Mat& C);
Mat  gemm(const Mat& A, const Mat& B);

#endif