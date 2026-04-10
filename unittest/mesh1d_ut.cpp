// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for Mesh1D class
//
// -----------------------------------------------------------------------------

#include "mesh1d_ut.h"

// -----------------------------------------------------------------------------
// Description: Basic mesh construction test
// -----------------------------------------------------------------------------
void Test_Mesh1D_basic()
{
  Mesh1D mesh(0.0, 1.0, 4);

  Check(mesh.Ne() == 4, "Mesh should have 4 elements");

  CheckEqual(mesh.element(0).left(),  0.0,  1e-14, "Wrong left boundary of element 0");
  CheckEqual(mesh.element(0).right(), 0.25, 1e-14, "Wrong right boundary of element 0");

  CheckEqual(mesh.element(3).left(),  0.75, 1e-14, "Wrong left boundary of element 3");
  CheckEqual(mesh.element(3).right(), 1.0,  1e-14, "Wrong right boundary of element 3");
}