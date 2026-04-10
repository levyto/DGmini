// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Definitions of linear algebra operations (BLAS-like)
//
// -----------------------------------------------------------------------------

#ifndef LINALG_H
#define LINALG_H

#include "Vec.h"

double dot(const Vec& x, const Vec& y);
void scal(double alpha, Vec& x);
void axpy(double alpha, const Vec& x, Vec& y);

#endif