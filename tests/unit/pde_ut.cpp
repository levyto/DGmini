// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for Quadrature1D class
//
// -----------------------------------------------------------------------------

#include <cmath>

#include "pde_ut.h"
#include "PDE/linear_advection1d.h"
#include "PDE/burgers1d.h"

// -----------------------------------------------------------------------------
// Description: Test LinearAdvection1D flux and eigenvalue
// -----------------------------------------------------------------------------
void Test_LinearAdvection1D()
{
  const double tol = 1.e-14;

  LinearAdvection1D pde(2.5);

  CheckEqual(pde.convectiveFlux(3.0), 7.5, tol, "Linear advection flux failed");
  CheckEqual(pde.convectiveFlux(-2.0), -5.0, tol, "Linear advection flux failed");
  CheckEqual(pde.maxConvectiveEigenvalue(0.0), 2.5, tol, "Linear advection eigenvalue failed");
  CheckEqual(pde.maxConvectiveEigenvalue(10.0), 2.5, tol, "Linear advection eigenvalue failed");
}

// -----------------------------------------------------------------------------
// Description: Test Burgers1D flux and eigenvalue
// -----------------------------------------------------------------------------
void Test_Burgers1D()
{
  const double tol = 1.e-14;

  Burgers1D pde;

  CheckEqual(pde.convectiveFlux(2.0), 2.0, tol, "Burgers flux failed");
  CheckEqual(pde.convectiveFlux(-2.0), 2.0, tol, "Burgers flux failed");
  CheckEqual(pde.maxConvectiveEigenvalue(3.0), 3.0, tol, "Burgers eigenvalue failed");
  CheckEqual(pde.maxConvectiveEigenvalue(-4.0), 4.0, tol, "Burgers eigenvalue failed");
}
