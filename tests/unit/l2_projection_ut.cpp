// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for L2 projection
//
// -----------------------------------------------------------------------------

#include <cmath>
#include <sstream>

#include "l2_projection_ut.h"
#include "FEM/basis1d.h"
#include "FEM/quadrature1d.h"

// -----------------------------------------------------------------------------
// Description: L2 projection of a constant function should yield the constant 
//              value in the first mode and zero in all higher modes
// -----------------------------------------------------------------------------
void Test_L2ProjectionOnElement_constant()
{
  const int p = 5;
  const double C = 3.5;

  FESpace1D fe(p);
  Element1D element(2.0, 5.0);

  Vec u_e(fe.DoFs());

  auto u0 = [C](double x) { return C; };

  L2ProjectionOnElement(fe, element, u0, u_e);

  CheckEqual(u_e[0], C, 1.e-12, "Constant projection failed at mode 0");

  for (int i = 1; i <= p; ++i)
  {
    std::ostringstream oss;
    oss << "Constant projection failed at mode " << i;
    CheckEqual(u_e[i], 0.0, 1.e-12, oss.str());
  }
}

// -----------------------------------------------------------------------------
// Description: L2 projection of a linear function should yield nonzero value
//              in the first two modes and zero in all higher modes
// -----------------------------------------------------------------------------
void Test_L2ProjectionOnElement_linear()
{
  const int p = 5;

  FESpace1D fe(p);
  Element1D element(0.0, 2.0);

  Vec u_e(fe.DoFs());

  auto u0 = [](double x) { return x; };

  L2ProjectionOnElement(fe, element, u0, u_e);

  Check(std::abs(u_e[0]) > 0.0, "Linear projection failed: mode 0 is zero");
  Check(std::abs(u_e[1]) > 0.0, "Linear projection failed: mode 1 is zero");

  for (int i = 2; i <= p; ++i)
  {
    std::ostringstream oss;
    oss << "Linear projection failed at higher mode " << i;
    CheckEqual(u_e[i], 0.0, 1.e-12, oss.str());
  }
}

// -----------------------------------------------------------------------------
// Description: L2 projection of a Legendre mode should yield 1.0 in the 
//              corresponding mode and zero in all other modes
// -----------------------------------------------------------------------------
void Test_L2ProjectionOnElement_legendreMode()
{
  const int p = 5;
  const int k = 2;

  FESpace1D fe(p);
  Element1D element(2.0, 5.0);

  Vec u_e(p + 1);

  auto u0 = [&element, k](double x)
  {
    const double xi = element.mapToReference(x);

    Vec  phi(k + 1);
    Vec dphi(k + 1);

    evaluateLegendreBasis(k, xi, phi, dphi);

    return phi[k];
  };

  L2ProjectionOnElement(fe, element, u0, u_e);

  for (int i = 0; i <= p; ++i)
  {
    const double expected = (i == k) ? 1.0 : 0.0;

    std::ostringstream oss;
    oss << "Legendre mode projection failed at mode " << i;
    CheckEqual(u_e[i], expected, 1.e-12, oss.str());
  }
}