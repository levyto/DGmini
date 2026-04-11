// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for stiffness  matrix on one element
//
// -----------------------------------------------------------------------------

#include <sstream>

#include "stiffness_matrix_ut.h"
#include "quadrature1d.h"

// -----------------------------------------------------------------------------
// Description: Check if numerically integrated stiffness matrix using quadrature 
//              matches the analytical one computed using the orthogonality of 
//              Legendre polynomials
// -----------------------------------------------------------------------------
void Test_stiffnessMatrix1D_numVsAnalytic()
{
  const int p = 5;
  const double tol = 1e-12;

  Quadrature1D quadrature(2 * p);

  Mat K_num = buildStiffnessMatrix1D(quadrature, p);
  Mat K_ref = buildStiffnessMatrix1D(p);

  for (int i = 0; i <= p; ++i)
  {
    for (int j = 0; j <= p; ++j)
    {
      std::ostringstream oss;
      oss << "Stiffness matrix mismatch at (" << i << "," << j << ")";

      CheckEqual(K_num(i, j), K_ref(i, j), tol, oss.str());
    }
  }
}
