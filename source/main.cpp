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
  ModalVector u(mesh.Ne(), fe.DoFs());
  ModalVector rhs(mesh.Ne(), fe.DoFs());

  // ---------------------------------------------------------------------------
  // Initial condition
  // ---------------------------------------------------------------------------
  auto u0 = [](double x)
  {
    return std::sin(2.0 * M_PI * x);
  };

  for (int e = 0; e < mesh.Ne(); ++e)
  {
    L2ProjectionOnElement(fe, mesh.element(e), u0, u.elementPtr(e));
  }

  // Write initial modal solution
  writeModalSolution1D("solution_0000.dat", mesh, u);
  int output_id = 1;

  // ---------------------------------------------------------------------------
  // Time integration: explicit Euler
  // ---------------------------------------------------------------------------
  double time = 0.0;
  int step = 0;

  while (time < final_time)
  {
    residual(fe, mesh, *pde, *flux, u, rhs);

    // u^{n+1} = u^n + dt * rhs
    u.axpy(dt, rhs);

    if (step % 5 == 0)
    {
      std::ostringstream name;
      name << "solution_" << std::setw(4) << std::setfill('0') << output_id << ".dat";
      writeModalSolution1D(name.str(), mesh, u);
      ++output_id;
    }

    time += dt;
    ++step;
  }

  // ---------------------------------------------------------------------------
  // Output
  // ---------------------------------------------------------------------------
  writeModalSolution1D("solution_final.dat", mesh, u);

  std::cout << "Finished at t = " << time
            << " after " << step << " steps.\n";

  return 0;
}
