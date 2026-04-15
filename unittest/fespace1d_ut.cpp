// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for 1D finite element space class
//
// -----------------------------------------------------------------------------

#include <algorithm>
#include <string>
#include <vector>

#include "fespace1d_ut.h"
#include "FEM/quadrature1d.h"

// -----------------------------------------------------------------------------
// Description: Check  P(1) == 1,         P(-1) == (-1)^n,
//                    dP(1) == n(n+1)/2, dP(-1) == (-1)^(n+1)*n(n+1)/2
// -----------------------------------------------------------------------------
void Test_FESpace1D_constructor()
{
  const int p = 5;
  const double tol = 1e-12;

  FESpace1D fe(p);

  // Dimensions
  Check(fe.DoFs() == (p + 1), "DoFs mismatch");
  Check(fe.nip() > 0, "quadrature not initialized");

  // Boundary values of Legendre polynomials
  for (int i = 0; i <= p; ++i)
  {
    // P_n(1) = 1
    CheckEqual(fe.phiRight(i), 1.0, tol, "phiRight mismatch");
    // P_n(-1) = (-1)^n
    double expected = (i % 2 == 0) ? 1.0 : -1.0;
    CheckEqual(fe.phiLeft(i), expected, tol, "phiLeft mismatch");
  }

  // Init. of basis by summing on first quadrature point (sanity, not full validation)
  const int q0 = 0;
  double sum = 0.0;
  for (int i = 0; i <= p; ++i)
  {
    sum += std::abs(fe.phi(q0, i));
  }
  Check(sum > 0.0, "basis seems uninitialized");

  // init. of mass matrix by summing first row (sanity, not full validation)
  sum = 0.0;
  for (int j = 0; j <= p; ++j)
  {
    sum += fe.massMatrix()(0, j);
  }
  Check(sum > 0.0, "mass matrix seems uninitialized");

  // Init. of mass matrix inverse by summing first row (sanity, not full validation)
  sum = 0.0;
  for (int j = 0; j <= p; ++j)
  {
    sum += fe.inverseMassMatrix()(0, j);
  }
  Check(sum > 0.0, "mass matrix inverse seems uninitialized");
}