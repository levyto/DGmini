// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Light-weight 1D DG reference implementation intended for
//              experimentation and understanding of DG building blocks
//
// -----------------------------------------------------------------------------

#include <iostream>

#include "FEM/fespace1d.h"
#include "IO/input.h"
#include "IO/output.h"
#include "Mesh/mesh1d.h"

using std::cout;
using ModalSolution1D = std::vector<Vec>;

int main()
{
  std::cout << "DGmini: startup OK\n";

  auto pde = createPDE("linear_advection1d");

  const int Ne = 4;
  const double x0 = 0.0;
  const double x1 = 1.0;
  const int p = 2;

  FESpace1D fe(p);
  Mesh1D mesh(x0, x1, Ne);

  ModalSolution1D sol(mesh.Ne(), Vec(fe.DoFs()));

  for (int i = 0; i < mesh.Ne(); ++i)
  {
    cout << "Element " << i
         << ": ["      << mesh.element(i).left()
         << ", "       << mesh.element(i).right()
         << "]\n";
  }

  sol[0][0] = 1.0; sol[0][1] = 0.1; sol[0][2] = 0.0;
  sol[1][0] = 0.8; sol[1][1] = 0.2; sol[1][2] = 0.0;
  sol[2][0] = 0.6; sol[2][1] = 0.3; sol[2][2] = 0.1;
  sol[3][0] = 0.4; sol[3][1] = 0.1; sol[3][2] = 0.0;

  writeModalSolution1D("solution.dat", mesh, sol);

  return 0;
}
