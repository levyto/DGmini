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

void buildMassMatrix1D(const Quadrature1D& quadrature, int p, Mat& M_e);
Mat  buildMassMatrix1D(const Quadrature1D& quadrature, int p);
void buildMassMatrix1D(int p, Mat& M_e);
Mat  buildMassMatrix1D(int p);
void buildMassMatrix1DInverse(int p, Mat& M_e);
Mat  buildMassMatrix1DInverse(int p);

#endif