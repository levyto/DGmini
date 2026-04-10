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

  return 0;
}
