// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Uniform 1D mesh class
//
// -----------------------------------------------------------------------------

#include <cassert>

#include "mesh1d.h"

// -----------------------------------------------------------------------------
// Description: Mesh1D constructor
// -----------------------------------------------------------------------------
Mesh1D::Mesh1D(double x0, double x1, int Ne)
{
  assert(x0 < x1);
  assert(Ne > 0);

  const double delta_x = (x1 - x0) / Ne;

  for (int i = 0; i < Ne; ++i)
  {
    const double xl = x0 +    i    * delta_x;
    const double xr = x0 + (i + 1) * delta_x;
    elements_.push_back(Element1D(xl, xr));
  }
}


// -----------------------------------------------------------------------------
// Description: Access element by index
// -----------------------------------------------------------------------------
const Element1D& Mesh1D::element(int i) const
{
  assert(i >= 0 && i < Ne());
  return elements_[i];
}
