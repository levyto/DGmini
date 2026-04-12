// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: L2 projection of a prescribed function onto the space 
//              of piecewise polynomials
//
// -----------------------------------------------------------------------------

#ifndef L2_PROJECTION_H
#define L2_PROJECTION_H

#include <cassert>
#include <functional>

#include "FEM/fespace1d.h"
#include "Mesh/element1d.h"

void L2ProjectionOnElement(const FESpace1D& fe, 
                           const Element1D& element,
                           std::function<double(double)> u0,
                           double* u_e);

void L2ProjectionOnElement(const FESpace1D& fe, 
                           const Element1D& element,
                           std::function<double(double)> u0,
                           Vec& u_e);

#endif