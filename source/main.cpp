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
    residual(fe, mesh, *pde, *flux, sol, rhs);

    sol.axpy(dt, rhs);    // u^{n+1} = u^n + dt * rhs

    time += dt;
    step++;

    output.write(mesh, sol, time);
  }

  output.writeFinal(mesh, sol, time);

  std::cout << "\n\nFinished at t = " << time
            << " after " << step << " steps.\n";

  return 0;
}
