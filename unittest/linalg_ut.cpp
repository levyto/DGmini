// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for linear algebra operations
//
// -----------------------------------------------------------------------------

#include "linalg_ut.h"
#include "Vec.h"
#include "Mat.h"

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
// Description: (Vec) y <- alpha*x + y test
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

// -----------------------------------------------------------------------------
// Description: (Mat-Vec) y <- A*x test
// -----------------------------------------------------------------------------
void Test_linalg_gemv()
{
  Mat A(2, 3);
  Vec x(3);
  Vec y(2);

  A(0,0) = 1.0; A(0,1) = 2.0; A(0,2) = 3.0;
  A(1,0) = 4.0; A(1,1) = 5.0; A(1,2) = 6.0;

  x[0] = 1.0;
  x[1] = 0.0;
  x[2] = -1.0;

  gemv(A, x, y);

  CheckEqual(y[0], -2.0, 1e-14, "gemv failed at row 0");
  CheckEqual(y[1], -2.0, 1e-14, "gemv failed at row 1");
}

// -----------------------------------------------------------------------------
// Description: (Mat-Vec) y = A*x test
// -----------------------------------------------------------------------------
void Test_linalg_gemvHighLvl()
{
  Mat A(2, 3);
  Vec x(3);

  A(0,0) = 1.0; A(0,1) = 2.0; A(0,2) = 3.0;
  A(1,0) = 4.0; A(1,1) = 5.0; A(1,2) = 6.0;

  x[0] = 1.0;
  x[1] = 0.0;
  x[2] = -1.0;

  Vec y = gemv(A, x);

  CheckEqual(y[0], -2.0, 1e-14, "high-level gemv failed at row 0");
  CheckEqual(y[1], -2.0, 1e-14, "high-level gemv failed at row 1");
}

// -----------------------------------------------------------------------------
// Description: (Mat-Mat) C <- A*B test
// -----------------------------------------------------------------------------
void Test_linalg_gemm()
{
  Mat A(2, 2);
  Mat B(2, 2);
  Mat C(2, 2);

  A(0,0) = 1.0; A(0,1) = 2.0;
  A(1,0) = 3.0; A(1,1) = 4.0;

  B(0,0) = 5.0; B(0,1) = 6.0;
  B(1,0) = 7.0; B(1,1) = 8.0;

  gemm(A, B, C);

  CheckEqual(C(0,0), 19.0, 1e-14, "gemm failed at (0,0)");
  CheckEqual(C(0,1), 22.0, 1e-14, "gemm failed at (0,1)");
  CheckEqual(C(1,0), 43.0, 1e-14, "gemm failed at (1,0)");
  CheckEqual(C(1,1), 50.0, 1e-14, "gemm failed at (1,1)");
}

// -----------------------------------------------------------------------------
// Description: (Mat-Mat) C = A*B test
// -----------------------------------------------------------------------------
void Test_linalg_gemmHighLvl()
{
  Mat A(2, 2);
  Mat B(2, 2);

  A(0,0) = 1.0; A(0,1) = 2.0;
  A(1,0) = 3.0; A(1,1) = 4.0;

  B(0,0) = 5.0; B(0,1) = 6.0;
  B(1,0) = 7.0; B(1,1) = 8.0;

  Mat C = gemm(A, B);

  CheckEqual(C(0,0), 19.0, 1e-14, "high-level gemm failed at (0,0)");
  CheckEqual(C(0,1), 22.0, 1e-14, "high-level gemm failed at (0,1)");
  CheckEqual(C(1,0), 43.0, 1e-14, "high-level gemm failed at (1,0)");
  CheckEqual(C(1,1), 50.0, 1e-14, "high-level gemm failed at (1,1)");
}
