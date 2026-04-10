// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for Element1D class
//
// -----------------------------------------------------------------------------

#include "element1d_ut.h"

// -----------------------------------------------------------------------------
// Description: Reference-to-physical mapping test
// -----------------------------------------------------------------------------
void Test_Element1D_basic()
{
  Element1D e(2.0, 4.0);

  CheckEqual(e.mapToPhysical(-1.0), 2.0, 1e-14, "Mapping failed at xi = -1");
  CheckEqual(e.mapToPhysical( 0.0), 3.0, 1e-14, "Mapping failed at xi = 0");
  CheckEqual(e.mapToPhysical( 1.0), 4.0, 1e-14, "Mapping failed at xi = 1");

  CheckEqual(e.length(),   2.0, 1e-14, "Element length is wrong");
  CheckEqual(e.jacobian(), 1.0, 1e-14, "Element Jacobian is wrong");
}