// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Light-weight 1D DG reference implementation intended for
//              experimentation and understanding of DG building blocks
//
// -----------------------------------------------------------------------------

#include <iostream>

#include "mesh1d.h"
#include "output.h"

using std::cout;

int main()
{
  std::cout << "DGmini: startup OK\n";

  Mesh1D mesh(0.0, 1.0, 4);

  for (int i = 0; i < mesh.Ne(); ++i)
  {
    cout << "Element " << i
         << ": ["      << mesh.element(i).left()
         << ", "       << mesh.element(i).right()
         << "]\n";
  }

  std::vector<Vec> coeffs(mesh.Ne(), Vec(3));

  coeffs[0][0] = 1.0; coeffs[0][1] = 0.1; coeffs[0][2] = 0.0;
  coeffs[1][0] = 0.8; coeffs[1][1] = 0.2; coeffs[1][2] = 0.0;
  coeffs[2][0] = 0.6; coeffs[2][1] = 0.3; coeffs[2][2] = 0.1;
  coeffs[3][0] = 0.4; coeffs[3][1] = 0.1; coeffs[3][2] = 0.0;

  writeModalSolution1D("solution.dat", mesh, coeffs);

  return 0;
}
