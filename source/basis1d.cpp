// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: 1D basis functions and their derivatives
//
// -----------------------------------------------------------------------------

#include <cassert>
#include <cmath>

#include "basis1d.h"

// -----------------------------------------------------------------------------
// Description: Evaluate the Legendre basis functions and their derivatives
//              at a given point xi from [-1,1] where the basis functions are 
//              defined. Note this basis is modal, i.e., the basis functions 
//              are defined on the reference element and are not associated with 
//              any particular element. 
//
// Here: p  ... order of the basis functions (number of basis functions = p + 1)
//       xi ... point in [-1,1] where the basis functions are evaluated
//       P  ... vector of size p + 1 to store the basis function values
//       dP ... vector of size p + 1 to store the basis function derivatives
// -----------------------------------------------------------------------------
void evaluateLegendreBasis(int p, double xi, Vec& P, Vec& dP)
{
  assert(P.size()  == p + 1);
  assert(dP.size() == p + 1);

   P[0] = 1.0;
  dP[0] = 0.0;

  if (p == 0)
    return;

   P[1] = xi;
  dP[1] = 1.0;

  for (int n = 1; n < p; ++n)
  {
     P[n + 1] = ((2.0 * n + 1.0) * xi * P[n] - n * P[n - 1]) / (n + 1.0);
    dP[n + 1] = (n + 1.0) * P[n] + xi * dP[n];
  }
}
