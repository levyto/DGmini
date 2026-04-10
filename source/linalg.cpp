// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Definitions of linear algebra operations (BLAS-like)
//
// -----------------------------------------------------------------------------

#include <cassert>
#include <cmath>

#include "linalg.h"

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
