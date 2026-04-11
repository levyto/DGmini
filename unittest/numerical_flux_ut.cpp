// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for Quadrature1D class
//
// -----------------------------------------------------------------------------

#include <cmath>

#include "pde_ut.h"
#include "Spatial/lax_friedrichs_flux.h"
#include "Spatial/rusanov_flux.h"
#include "PDE/linear_advection1d.h"

// -----------------------------------------------------------------------------
// Description: Test Lax-Friedrichs flux for linear advection, 
//              should reduce to central flux with alpha = max eigenvalue
// -----------------------------------------------------------------------------
void Test_LaxFriedrichsFlux()
{
  const double tol = 1e-14;

  LinearAdvection1D pde(2.0);
  LaxFriedrichsFlux flux(1.0);

  const double uL = 3.0;
  const double uR = -1.0;

  const double fhat = flux.evaluate(pde, uL, uR);

  CheckEqual(fhat, 4.0, tol, "Lax-Friedrichs flux failed");
}

// -----------------------------------------------------------------------------
// Description: Test Rusanov flux for positive advection speed, 
//              should reduce to upwind flux
// -----------------------------------------------------------------------------
void Test_RusanovFlux_advectionPositive()
{
  LinearAdvection1D pde(2.0);
  RusanovFlux flux;

  const double uL = 3.0;
  const double uR = -1.0;

  const double fhat = flux.evaluate(pde, uL, uR);

  CheckEqual(fhat, 2.0 * uL, 1.e-14, "Rusanov flux failed for positive advection speed");
}

// -----------------------------------------------------------------------------
// Description: Test Rusanov flux for negative advection speed, 
//              should reduce to upwind flux
// -----------------------------------------------------------------------------
void Test_RusanovFlux_advectionNegative()
{
  LinearAdvection1D pde(-2.0);
  RusanovFlux flux;

  const double uL = 3.0;
  const double uR = -1.0;

  const double fhat = flux.evaluate(pde, uL, uR);

  CheckEqual(fhat, -2.0 * uR, 1.e-14, "Rusanov flux failed for negative advection speed");
}
