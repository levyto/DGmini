// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: 1D basis functions and their derivatives
//
// -----------------------------------------------------------------------------

#ifndef BASIS1D_H
#define BASIS1D_H

#include "Algebra/Vec.h"

void evaluateLegendreBasis(int p, double xi, Vec& P, Vec& dP);

#endif