// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Light-weight 1D DG reference implementation intended for
//              experimentation and understanding of DG building blocks
//
// -----------------------------------------------------------------------------

#include <cmath>
#include <iostream>

#include "FEM/fespace1d.h"
#include "IO/input.h"
#include "IO/output.h"
#include "Mesh/mesh1d.h"
#include "Spatial/l2_projection.h"
#include "Spatial/modal_vector.h"

using std::cout;

int main()
{
  std::cout << "DGmini: startup OK\n";

  auto pde = createPDE("linear_advection1d");

  auto numFlux = createNumericalFlux("rusanov");

  const int Ne = 16;
  const double x0 = 0.0;
  const double x1 = 6.283185307179586; // 2*pi
  const int p = 2;

  FESpace1D   fe(p);

  Mesh1D      mesh(x0, x1, Ne);
  
  ModalVector sol(mesh.Ne(), fe.DoFs());

  for (int e = 0; e < mesh.Ne(); ++e)
  {
    cout << "Element " << e
         << ": ["      << mesh.element(e).left()
         << ", "       << mesh.element(e).right()
         << "]\n";
  }

  auto u0 = [](double x) { return std::sin(x); };

  for (int e = 0; e < mesh.Ne(); e++)
  {
    L2ProjectionOnElement(fe, mesh.element(e), u0, sol.elementPtr(e));
  }
  

  // sol[0][0] = 1.0; sol[0][1] = 0.1; sol[0][2] = 0.0;
  // sol[1][0] = 0.8; sol[1][1] = 0.2; sol[1][2] = 0.0;
  // sol[2][0] = 0.6; sol[2][1] = 0.3; sol[2][2] = 0.1;
  // sol[3][0] = 0.4; sol[3][1] = 0.1; sol[3][2] = 0.0;

  writeModalSolution1D("solution.dat", mesh, sol);

  return 0;
}
