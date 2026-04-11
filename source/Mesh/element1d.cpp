// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Element geometry class
//
// -----------------------------------------------------------------------------

#include <cassert>

#include "Mesh/element1d.h"

// -----------------------------------------------------------------------------
// Description: Element1D constructor
// -----------------------------------------------------------------------------
Element1D::Element1D(double x_left, double x_right)
  : x_left_(x_left), x_right_(x_right)
{
  assert(x_left < x_right);
}

// -----------------------------------------------------------------------------
// Description: Map xi in [-1,1] to physical coordinate x
// -----------------------------------------------------------------------------
double Element1D::mapToPhysical(double xi) const
{
  return 0.5 * (x_left_ + x_right_) + 0.5 * (x_right_ - x_left_) * xi;
}

// -----------------------------------------------------------------------------
// Description: Map x in [x_left_, x_right_] to reference coordinate xi in [-1,1]
// -----------------------------------------------------------------------------
double Element1D::mapToReference(double x) const
{
  return (2.0 * x - (x_left_ + x_right_)) / (x_right_ - x_left_);
}