// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Mass matrix on reference element
//
// -----------------------------------------------------------------------------

#ifndef MASSMATRIX_H
#define MASSMATRIX_H

#include "element1d.h"
#include "quadrature1d.h"
#include "Mat.h"

Mat buildMassMatrix1D(const Quadrature1D& quadrature, int p);
Mat buildMassMatrix1D(int p);
Mat buildMassMatrix1DInverse(int p);

#endif