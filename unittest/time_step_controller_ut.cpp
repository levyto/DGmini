// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for TimeStepController class
//
// -----------------------------------------------------------------------------

#include <cmath>

#include "time_step_controller_ut.h"
#include "FEM/fespace1d.h"
#include "Spatial/modal_vector.h"
#include "Mesh/mesh1d.h"
#include "Temporal/fixed_time_step.h"
#include "Temporal/cfl_time_step.h"
#include "PDE/linear_advection1d.h"
#include "PDE/burgers1d.h"

// -----------------------------------------------------------------------------
// Description: Test that FixedTimeStep controller returns the prescribed time 
//              step value
// -----------------------------------------------------------------------------
void Test_TimeStepController_FixedTimeStepReturnsGivenValue()
{
  const double dt = 1.25e-3;

  FESpace1D fe(2);
  Mesh1D mesh(0.0, 1.0, 4);
  LinearAdvection1D pde(2.0);
  ModalVector u(mesh.Ne(), fe.DoFs());
  u.zero();

  FixedTimeStep controller(dt);

  const double computed_dt = controller.computeTimeStep(fe, mesh, pde, u);

  CheckEqual(computed_dt, dt, 1e-15, 
             "FixedTimeStep did not return the prescribed dt");
}

// -----------------------------------------------------------------------------
// Description: Test that CFLTimeStep controller returns the correct time step 
//              for linear advection with a constant solution, which should be 
//              the same as the CFL time step based on the maximum eigenvalue of 
//              the flux Jacobian and the mesh size
// -----------------------------------------------------------------------------
void Test_TimeStepController_CFLTimeStepLinearAdvectionConstantSolution()
{
  const double cfl = 0.5;
  const double a   = 2.0;
  const int p      = 2;
  const int Ne     = 4;

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

  CFLTimeStep controller(cfl);

  const double h = mesh.element(0).right() - mesh.element(0).left();
  const double expected_dt = cfl * h / ((2.0 * fe.order() + 1.0) * std::abs(a));

  const double computed_dt = controller.computeTimeStep(fe, mesh, pde, u);

  CheckEqual(computed_dt, expected_dt, 1e-14, 
             "CFLTimeStep returned wrong dt for linear advection");
}

// -----------------------------------------------------------------------------
// Description: Test that CFLTimeStep controller returns the correct time step
//              for Burgers' equation with a piecewise constant solution, which
//              should be the same as the CFL time step based on the maximum
//              eigenvalue of the flux Jacobian (which is the solution value for
//              Burgers) and the mesh size
// -----------------------------------------------------------------------------
void Test_TimeStepController_CFLTimeStepBurgersP0UsesMaxSolutionValue()
{
  const double cfl = 0.5;
  const int p = 0;
  const int Ne = 3;

  FESpace1D fe(p);
  Mesh1D mesh(0.0, 3.0, Ne); // h = 1
  Burgers1D pde;
  ModalVector u(mesh.Ne(), fe.DoFs());

  u(0,0) = 1.0;
  u(1,0) = -3.0;
  u(2,0) = 2.0;

  CFLTimeStep controller(cfl);

  const double expected_dt = cfl * 1.0 / ((2.0 * fe.order() + 1.0) * 3.0);
  const double computed_dt = controller.computeTimeStep(fe, mesh, pde, u);

  CheckEqual(computed_dt, expected_dt, 1e-14,
             "CFLTimeStep returned wrong dt for Burgers p=0");
}

// -----------------------------------------------------------------------------
// Description: Test that CFLTimeStep controller correctly accounts for the
//              solution values at the element boundaries when computing the
//              maximum eigenvalue
// -----------------------------------------------------------------------------
void Test_TimeStepController_CFLTimeStepReflectsElementBoundaries()
{
  const double cfl = 0.5;
  const int      p = 1;
  const int     Ne = 1;

  FESpace1D fe(p);
  Mesh1D mesh(0.0, 1.0, Ne);
  Burgers1D pde;
  ModalVector u(mesh.Ne(), fe.DoFs());
  u.zero();
 
  // Linear solution on element 0 (from -1 to 1)
  u(0,0) = 0.0;
  u(0,1) = 1.0;

  CFLTimeStep controller(cfl);

  const double h = 1.0;
  const double expected_dt = cfl * h / ((2.0 * fe.order() + 1.0) * 1.0);

  const double computed_dt = controller.computeTimeStep(fe, mesh, pde, u);

  CheckEqual(computed_dt, expected_dt, 1e-14,
             "CFLTimeStep should use boundary values when computing max eigenvalue");
}