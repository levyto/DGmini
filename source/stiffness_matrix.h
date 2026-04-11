// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Stiffness matrix on one element
//
// -----------------------------------------------------------------------------

#ifndef STIFFNESSMATRIX_H
#define STIFFNESSMATRIX_H

#include "quadrature1d.h"
#include "Mat.h"

Mat buildStiffnessMatrix1D(const Quadrature1D& quadrature, int p);
Mat buildStiffnessMatrix1D(int p);

#endif