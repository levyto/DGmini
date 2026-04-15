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
