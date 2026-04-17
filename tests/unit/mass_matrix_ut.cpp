// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for mass matrix on one element
//
// -----------------------------------------------------------------------------

#include <sstream>

#include "mass_matrix_ut.h"
#include "FEM/quadrature1d.h"
#include "Algebra/linalg.h"

// -----------------------------------------------------------------------------
// Description: Check if numerically integrated mass matrix using quadrature 
//              matches the analytical one computed using the orthogonality of 
//              Legendre polynomials
// -----------------------------------------------------------------------------
void Test_massMatrix1D_numVsAnalyticLowLvl()
{
  const int p = 4;
  const double tol = 1e-12;

  Quadrature1D quadrature(2 * p);

  Mat M_num(p + 1, p + 1);
  buildMassMatrix1D(quadrature, p, M_num);
  Mat M_ref(p + 1, p + 1);
  buildMassMatrix1D(p, M_ref);

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
// Description: Check if numerically integrated mass matrix using quadrature 
//              matches the analytical one computed using the orthogonality of 
//              Legendre polynomials
// -----------------------------------------------------------------------------
void Test_massMatrix1D_numVsAnalyticHighLvl()
{
  const int p = 4;
  const double tol = 1e-12;

  Quadrature1D quadrature(2 * p);

  Mat M_num = buildMassMatrix1D(quadrature, p);
  Mat M_ref = buildMassMatrix1D(p);

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
void Test_massMatrix1D_inverseLowLvl()
{
  const int p = 5;
  const double tol = 1e-12;

  Mat M(p + 1, p + 1);
  buildMassMatrix1D(p, M);
  Mat M_inv(p + 1, p + 1);
  buildMassMatrix1DInverse(p, M_inv);

  Mat I(p + 1, p + 1);
  gemm(M_inv, M, I);

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

// -----------------------------------------------------------------------------
// Description: Check if inv(M)*M == I, where I is the identity matrix
// -----------------------------------------------------------------------------
void Test_massMatrix1D_inverseHighLvl()
{
  const int p = 5;
  const double tol = 1e-12;

  Mat M     = buildMassMatrix1D(p);
  Mat M_inv = buildMassMatrix1DInverse(p);
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
