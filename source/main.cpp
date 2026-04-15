// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Light-weight 1D DG reference implementation intended for
//              experimentation and understanding of DG building blocks
//
// -----------------------------------------------------------------------------

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "FEM/fespace1d.h"
#include "IO/input.h"
#include "IO/output.h"
#include "Mesh/mesh1d.h"
#include "PDE/pde.h"
#include "Spatial/l2_projection.h"
#include "Spatial/modal_vector.h"
#include "Spatial/numerical_flux.h"
#include "Spatial/residual.h"
#include "Temporal/time_integrator.h"
#include "Temporal/time_step_controller.h"
#include "Temporal/cfl_number.h"

using std::cout;

int main()
{
  // ---------------------------------------------------------------------------
  // Problem setup
  // ---------------------------------------------------------------------------
  const int p = 2;
  const int n_elements = 20;
  const double x_left = 0.0;
  const double x_right = 1.0;

  const double final_time = 1.0;
  const double dt = 1.0e-3;

  // ---------------------------------------------------------------------------
  // Discretization objects
  // ---------------------------------------------------------------------------
  FESpace1D fe(p);
  Mesh1D mesh(x_left, x_right, n_elements);

  auto pde = createPDE("linear_advection1d");
  auto flux = createNumericalFlux("rusanov");

  auto integrator = createTimeIntegrator("forward_euler");
  // auto integrator = createTimeIntegrator("runge_kutta_2");
  // auto integrator = createTimeIntegrator("runge_kutta_3_ssp");
  // auto integrator = createTimeIntegrator("runge_kutta_2");
  integrator->initialize(mesh, fe);

  auto dt_controller = createTimeStepController("cfl_time_step", 0.53);
  // auto dt_controller = createTimeStepController("cfl_time_step", 0.45);

  // ---------------------------------------------------------------------------
  // Solution vectors
  // ---------------------------------------------------------------------------
  ModalVector sol(mesh.Ne(), fe.DoFs());
  ModalVector rhs(mesh.Ne(), fe.DoFs());

  // ---------------------------------------------------------------------------
  // Output
  // ---------------------------------------------------------------------------
  TimeSeriesWriter output("output", "solution", 0.01);

  // ---------------------------------------------------------------------------
  // Initial condition
  // ---------------------------------------------------------------------------
  auto u0 = [](double x)
  {
    return std::sin(2.0 * M_PI * x);
  };

  for (int e = 0; e < mesh.Ne(); ++e)
  {
    L2ProjectionOnElement(fe, mesh.element(e), u0, sol.elementPtr(e));
  }

  output.write(mesh, sol, 0.0);

  // ---------------------------------------------------------------------------
  // Time integration: explicit Euler
  // ---------------------------------------------------------------------------
  double time = 0.0;
  int step = 0;

  while (time < final_time)
  {
    double dt = dt_controller->computeTimeStep(fe, mesh, *pde, sol);
    
    const double effective_cfl = computeEffectiveCFL(fe, mesh, *pde, sol, dt);
    if (effective_cfl > integrator->recommendedCFL())
    {
      std::cout << "\nWarning: current CFL = " << effective_cfl
                << " exceeds recommended value "
                << integrator->recommendedCFL()
                << " for selected time integrator.";
    }

    integrator->doTimeStep(fe, mesh, *pde, *flux, dt, sol);

    time += dt;
    step++;

    output.write(mesh, sol, time);
  }

  output.writeFinal(mesh, sol, time);

  std::cout << "\n\nFinished at t = " << time
            << " after " << step << " steps.\n";

  return 0;
}
