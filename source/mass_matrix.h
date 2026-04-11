// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Mass matrix on one element
//
// -----------------------------------------------------------------------------

#ifndef MASSMATRIX_H
#define MASSMATRIX_H

#include "element1d.h"
#include "quadrature1d.h"
#include "Mat.h"

Mat buildMassMatrix1D(const Element1D& element, const Quadrature1D& quadrature, int p);
Mat buildMassMatrix1D(const Element1D& element, int p);
Mat buildMassMatrix1DInverse(const Element1D& element, int p);

#endif