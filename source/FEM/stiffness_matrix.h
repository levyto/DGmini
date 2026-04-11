// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Stiffness matrix on one element
//
// -----------------------------------------------------------------------------

#ifndef STIFFNESSMATRIX_H
#define STIFFNESSMATRIX_H

#include "Algebra/Mat.h"
#include "FEM/quadrature1d.h"

Mat buildStiffnessMatrix1D(const Quadrature1D& quadrature, int p);
Mat buildStiffnessMatrix1D(int p);

#endif