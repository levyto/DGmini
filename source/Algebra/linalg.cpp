// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Definitions of linear algebra operations (BLAS-like)
//
// -----------------------------------------------------------------------------

#include <cassert>
#include <cmath>

#include "Algebra/linalg.h"

// -----------------------------------------------------------------------------
// Description: L2 norm   ||x||_2
// -----------------------------------------------------------------------------
double nrm2(const Vec& x)
{
  return std::sqrt(dot(x, x));
}

// -----------------------------------------------------------------------------
// Description: Dot product   x^T y
// -----------------------------------------------------------------------------
double dot(const Vec& x, const Vec& y)
{
  assert(x.size() == y.size());

  double result = 0.0;

  for (int i = 0; i < x.size(); ++i)
  {
    result += x[i] * y[i];
  }

  return result;
}

// -----------------------------------------------------------------------------
// Description: Mutiplication by scalar   x <- alpha * x
// -----------------------------------------------------------------------------
void scal(double alpha, Vec& x)
{
  for (int i = 0; i < x.size(); ++i)
  {
    x[i] *= alpha;
  }
}

// -----------------------------------------------------------------------------
// Description: axpy   y <- alpha * x + y
// -----------------------------------------------------------------------------
void axpy(double alpha, const Vec& x, Vec& y)
{
  assert(x.size() == y.size());

  for (int i = 0; i < x.size(); ++i)
  {
    y[i] += alpha * x[i];
  }
}

// -----------------------------------------------------------------------------
// Description: Matrix-vector product y <- A x
// -----------------------------------------------------------------------------
void gemv(const Mat& A, const Vec& x, Vec& y)
{
  assert(A.cols() == x.size());
  assert(y.size() == A.rows());

  for (int i = 0; i < A.rows(); ++i)
  {
    double sum = 0.0;

    for (int j = 0; j < A.cols(); ++j)
    {
      sum += A(i, j) * x[j];
    }

    y[i] = sum;
  }
}

// -----------------------------------------------------------------------------
// Description: Matrix-vector product   A x
// -----------------------------------------------------------------------------
Vec gemv(const Mat& A, const Vec& x)
{
  assert(A.cols() == x.size());

  Vec result(A.rows());

  for (int i = 0; i < A.rows(); ++i)
  {
    double sum = 0.0;

    for (int j = 0; j < A.cols(); ++j)
    {
      sum += A(i, j) * x[j];
    }

  result[i] = sum;
  }

  return result;
}

// -----------------------------------------------------------------------------
// Description: Matrix-matrix product   C <- A B
// -----------------------------------------------------------------------------
void gemm(const Mat& A, const Mat& B, Mat& C)
{
  assert(A.cols() == B.rows());
  assert(C.cols() == B.cols());
  assert(C.rows() == A.rows());

  for (int i = 0; i < A.rows(); ++i)
  {
    for (int j = 0; j < B.cols(); ++j)
    {
      double sum = 0.0;

      for (int k = 0; k < A.cols(); ++k)
      {
        sum += A(i, k) * B(k, j);
      }

      C(i, j) = sum;
    }
  }
}

// -----------------------------------------------------------------------------
// Description: Matrix-matrix product   A B
// -----------------------------------------------------------------------------
Mat gemm(const Mat& A, const Mat& B)
{
  assert(A.cols() == B.rows());

  Mat result(A.rows(), B.cols());

  for (int i = 0; i < A.rows(); ++i)
  {
    for (int j = 0; j < B.cols(); ++j)
    {
      double sum = 0.0;

      for (int k = 0; k < A.cols(); ++k)
      {
        sum += A(i, k) * B(k, j);
      }

      result(i, j) = sum;
    }
  }

  return result;
}
