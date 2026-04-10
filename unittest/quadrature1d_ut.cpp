// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for Quadrature1D class
//
// -----------------------------------------------------------------------------

#include <cmath>

#include "quadrature1d_ut.h"

// -----------------------------------------------------------------------------
// Description: Sum of quadrature weights must be 2 on [-1,1]
// -----------------------------------------------------------------------------
void Test_Quadrature1D_weightsSum()
{
  for (int Q = 1; Q <= 39; ++Q)
  {
    Quadrature1D quadrature(Q);

    double sum = 0.0;
    for (int i = 0; i < quadrature.nip(); ++i)
    {
      sum += quadrature.weight(i);
    }

    CheckEqual(sum, 2.0, 1e-14, "Quadrature weights should sum to 2 on [-1,1]");
  }
}

// -----------------------------------------------------------------------------
// Description: Gauss-Legendre points and weights are symmetric
// -----------------------------------------------------------------------------
void Test_Quadrature1D_symmetry()
{
  for (int Q = 1; Q <= 39; ++Q)
  {
    Quadrature1D quadrature(Q);
    const int nip = quadrature.nip();

    for (int i = 0; i < nip; ++i)
    {
      CheckEqual(quadrature.point(i), -quadrature.point(nip - 1 - i), 1e-14,
                 "Quadrature points are not symmetric");

      CheckEqual(quadrature.weight(i), quadrature.weight(nip - 1 - i), 1e-14,
                 "Quadrature weights are not symmetric");
    }
  }
}

// -----------------------------------------------------------------------------
// Description: Exactness test on monomials x^k up to order Q,
//              integral is equal to 0 for odd k and equal to 2/(k+1) for even k
// -----------------------------------------------------------------------------
void Test_Quadrature1D_exactness()
{
  for (int Q = 1; Q <= 39; ++Q)
  {
    Quadrature1D quadrature(Q);

    for (int k = 0; k <= Q; ++k)
    {
      double numeric = 0.0;
      for (int i = 0; i < quadrature.nip(); ++i)
      {
        numeric += quadrature.weight(i) * std::pow(quadrature.point(i), k);
      }

      double exact = 0.0;
      if (k % 2 == 0)
      {
        exact = 2.0 / static_cast<double>(k + 1);
      }

      CheckEqual(numeric, exact, 1e-13, "Quadrature exactness failed for monomial");
    }
  }
}