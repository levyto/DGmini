// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Mass matrix on reference element
//
// -----------------------------------------------------------------------------

#ifndef MASSMATRIX_H
#define MASSMATRIX_H

#include "Algebra/Mat.h"
#include "FEM/quadrature1d.h"
#include "Mesh/element1d.h"

Mat buildMassMatrix1D(const Quadrature1D& quadrature, int p);
Mat buildMassMatrix1D(int p);
Mat buildMassMatrix1DInverse(int p);

#endif