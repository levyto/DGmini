// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for Vec class
//
// -----------------------------------------------------------------------------

#include "Vec_ut.h"

// -----------------------------------------------------------------------------
// Description: Basic Vec operations
// -----------------------------------------------------------------------------
void Test_Vec_basic()
{
  Vec v(3);

  Check(v.size() == 3, "Vec size should be 3");

  v[0] = 1.0;
  v[1] = 2.0;
  v[2] = 3.0;

  CheckEqual(v[0], 1.0, 1e-14, "Vec[0] mismatch");
  CheckEqual(v[1], 2.0, 1e-14, "Vec[1] mismatch");
  CheckEqual(v[2], 3.0, 1e-14, "Vec[2] mismatch");

  v.fill(5.0);

  CheckEqual(v[0], 5.0, 1e-14, "Vec fill failed at index 0");
  CheckEqual(v[1], 5.0, 1e-14, "Vec fill failed at index 1");
  CheckEqual(v[2], 5.0, 1e-14, "Vec fill failed at index 2");
}

