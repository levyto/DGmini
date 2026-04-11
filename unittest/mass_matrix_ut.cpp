// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for mass matrix on one element
//
// -----------------------------------------------------------------------------

#include <sstream>

#include "mass_matrix_ut.h"
#include "element1d.h"
#include "quadrature1d.h"
#include "linalg.h"

// -----------------------------------------------------------------------------
// Description: Check if numerically integrated mass matrix using quadrature 
//              matches the analytical one computed using the orthogonality of 
//              Legendre polynomials
// -----------------------------------------------------------------------------
void Test_massMatrix1D_numVsAnalytic()
{
  const int p = 4;
  const double tol = 1e-12;

  Element1D element(-2.0, 4.0);
  Quadrature1D quadrature(2 * p);

  Mat M_num = buildMassMatrix1D(element, quadrature, p);
  Mat M_ref = buildMassMatrix1D(element, p);

  for (int i = 0; i <= p; ++i)
  {
    for (int j = 0; j <= p; ++j)
    {
      std::ostringstream oss;
      oss << "Mass matrix mismatch at (" << i << "," << j << ")";

      CheckEqual(M_num(i, j), M_ref(i, j), tol, oss.str());
    }
  }
}

// -----------------------------------------------------------------------------
// Description: Check if inv(M)*M == I, where I is the identity matrix
// -----------------------------------------------------------------------------
void Test_massMatrix1D_inverse()
{
  const int p = 5;
  const double tol = 1e-12;

  Element1D element(1.1, 3.6);
  Mat M     = buildMassMatrix1D(element, p);
  Mat M_inv = buildMassMatrix1DInverse(element, p);
  Mat I     = gemm(M_inv, M);

  for (int i = 0; i <= p; ++i)
  {
    for (int j = 0; j <= p; ++j)
    {
      double expected = (i == j) ? 1.0 : 0.0;

      std::ostringstream oss;
      oss << "M_inv * M mismatch at (" << i << "," << j << ")";

      CheckEqual(I(i, j), expected, tol, oss.str());
    }
  }
}
