// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Light-weight 1D DG reference implementation intended for
//              experimentation and understanding of DG building blocks
//
// -----------------------------------------------------------------------------

#include <iostream>

#include "linalg.h"

using std::cout;

int main()
{

  Vec x(3);
  Vec y(3);

  x[0] = 1.0; x[1] = 2.0; x[2] = 3.0;
  y[0] = 4.0; y[1] = 5.0; y[2] = 6.0;

  double d = dot(x, y);   // 32

  cout << d << "\n";

  scal(2.0, x);           // x = [2,4,6]

  axpy(0.5, x, y);        // y = y + 0.5*x

  std::cout << "DGmini: startup OK\n";
  return 0;
}
