// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for TimeStepController class
//
// -----------------------------------------------------------------------------

#include <cmath>

#include "cfl_number_ut.h"
#include "FEM/fespace1d.h"
#include "Spatial/modal_vector.h"
#include "Mesh/mesh1d.h"
#include "PDE/linear_advection1d.h"
#include "PDE/burgers1d.h"

// -----------------------------------------------------------------------------
// Description: Test that evalMaxConvectiveEigenvalue returns the correct 
//              maximum eigenvalue for linear advection with a constant solution, 
//              which should be the absolute value of the advection velocity
// -----------------------------------------------------------------------------
void Test_evalMaxConvectiveEigenvalue_CFLTimeStepLinearAdvectionConstantSolution()
{
  const double a = 2.5;

  FESpace1D fe(2);
  Mesh1D mesh(0.0, 1.0, 4);
  LinearAdvection1D pde(a);

  ModalVector u(mesh.Ne(), fe.DoFs());

  for (int e = 0; e < mesh.Ne(); ++e)
  {
    u(e,0) = 1.0;
    for (int i = 1; i < fe.DoFs(); ++i)
      u(e,i) = 0.5;
  }

  const double lambda = evalMaxConvectiveEigenvalue(fe, mesh, pde, u);

  CheckEqual(lambda, std::abs(a), 1e-14, "Max eigenvalue for linear advection should be |a|");
}

// -----------------------------------------------------------------------------
// Description: Test that evalMaxConvectiveEigenvalue correctly accounts for 
//              Burgers' equation wave speed, which depends on the solution
// -----------------------------------------------------------------------------
void Test_evalMaxConvectiveEigenvalue_BurgersP0UsesMaxSolutionValue()
{
  const int p  = 0;
  const int Ne = 3;

  FESpace1D fe(p);
  Mesh1D mesh(0.0, 3.0, Ne);
  Burgers1D pde;

  ModalVector u(mesh.Ne(), fe.DoFs());

  u(0,0) = 1.0;
  u(1,0) = -3.0;
  u(2,0) = 2.0;

  const double lambda = evalMaxConvectiveEigenvalue(fe, mesh, pde, u);

  CheckEqual(lambda, 3.0, 1e-14, "Burgers p=0 should use max |u| over elements");
}

// -----------------------------------------------------------------------------
// Description: Test that evalMaxConvectiveEigenvalue correctly accounts for the 
//              fact that the maximum wave speed may occur at the element 
//              boundaries
// -----------------------------------------------------------------------------
void Test_evalMaxConvectiveEigenvalue_ReflectsElementBoundaries()
{
  const int p  = 1;
  const int Ne = 1;

  FESpace1D fe(p);
  Mesh1D mesh(0.0, 1.0, Ne);
  Burgers1D pde;

  ModalVector u(mesh.Ne(), fe.DoFs());
  u.zero();

  // Linear solution on element 0 (from -1 to 1)
  u(0,0) = 0.0;
  u(0,1) = 1.0;

  const double lambda = evalMaxConvectiveEigenvalue(fe, mesh, pde, u);

  CheckEqual(lambda, 1.0, 1e-14, "Max eigenvalue should include boundary values");
}

// -----------------------------------------------------------------------------
// Description: Test that computeEffectiveCFL returns the value predicted by the 
//              analytical formula based on the maximum wave speed, mesh size, 
//              approximation order and time step size
// -----------------------------------------------------------------------------
void Test_computeEffectiveCFL_matchesFormula()
{
  const double dt = 0.01;
  const double a  = 2.0;
  const int p     = 2;
  const int Ne    = 4;

  FESpace1D fe(p);
  Mesh1D mesh(0.0, 1.0, Ne);
  LinearAdvection1D pde(a);

  ModalVector u(mesh.Ne(), fe.DoFs());

  for (int e = 0; e < mesh.Ne(); ++e)
  {
    u(e,0) = 1.0;
    for (int i = 1; i < fe.DoFs(); ++i)
      u(e,i) = 0.0;
  }

  const double cfl = computeEffectiveCFL(fe, mesh, pde, u, dt);

  const double h = mesh.element(0).right() - mesh.element(0).left();
  const double expected = dt * ((2.0 * fe.order() + 1.0) * std::abs(a)) / h;

  CheckEqual(cfl, expected, 1e-14, "Effective CFL does not match analytical formula");
}