// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for Mat class
//
// -----------------------------------------------------------------------------

#include "Mat_ut.h"

// -----------------------------------------------------------------------------
// Description: Basic Mat operations
// -----------------------------------------------------------------------------
void Test_Mat_basic()
{
  Mat A(2, 3);

  Check(A.rows() == 2, "Mat rows mismatch");
  Check(A.cols() == 3, "Mat cols mismatch");

  A(1,2) = 7.5;
  CheckEqual(A(1,2), 7.5, 1e-14, "Mat entry access failed");

  A.fill(3.0);

  CheckEqual(A(0,0), 3.0, 1e-14, "Mat fill failed at (0,0)");
  CheckEqual(A(1,2), 3.0, 1e-14, "Mat fill failed at (1,2)");
}
