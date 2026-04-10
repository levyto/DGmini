// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for linear algebra operations
//
// -----------------------------------------------------------------------------

#include "linalg_ut.h"
#include "Vec.h"

// -----------------------------------------------------------------------------
// Description: (Vec) L2 norm test
// -----------------------------------------------------------------------------
void Test_linalg_nrm2()
{
  Vec x(2);
  x[0] = 3.0;
  x[1] = 4.0;

  CheckEqual(nrm2(x), 5.0, 1e-14, "nrm2 failed");
}

// -----------------------------------------------------------------------------
// Description: (Vec) Dot product test
// -----------------------------------------------------------------------------
void Test_linalg_dot()
{
  Vec x(3);
  x[0] = 1.0; 
  x[1] = 2.0; 
  x[2] = 3.0;

  Vec y(3);
  y[0] = 4.0; 
  y[1] = 5.0; 
  y[2] = 6.0;

  CheckEqual(dot(x, y), 32.0, 1e-14, "dot product failed");
}

// -----------------------------------------------------------------------------
// Description: (Vec) Multiplication by scalar test
// -----------------------------------------------------------------------------
void Test_linalg_scal()
{
  Vec x(3);
  x[0] = 1.0;
  x[1] = -2.0;
  x[2] = 3.0;

  scal(2.0, x);

  CheckEqual(x[0],  2.0, 1e-14, "scal failed at index 0");
  CheckEqual(x[1], -4.0, 1e-14, "scal failed at index 1");
  CheckEqual(x[2],  6.0, 1e-14, "scal failed at index 2");
}

// -----------------------------------------------------------------------------
// Description: (Vec) alpha*x + y test
// -----------------------------------------------------------------------------
void Test_linalg_axpy()
{
  Vec x(3);
  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;

  Vec y(3);
  y[0] = 4.0;
  y[1] = 5.0;
  y[2] = 6.0;

  axpy(0.5, x, y);

  CheckEqual(y[0], 4.5, 1e-14, "axpy failed at index 0");
  CheckEqual(y[1], 6.0, 1e-14, "axpy failed at index 1");
  CheckEqual(y[2], 7.5, 1e-14, "axpy failed at index 2");
}
