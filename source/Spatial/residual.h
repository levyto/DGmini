// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: du/dt = R residual assembly
//
// -----------------------------------------------------------------------------

#ifndef RESIDUAL_H 
#define RESIDUAL_H

#include "FEM/fespace1d.h"
#include "Mesh/mesh1d.h"
#include "PDE/pde.h"
#include "Spatial/modal_vector.h"
#include "Spatial/numerical_flux.h"

void residual(const FESpace1D& fe,
              const Mesh1D& mesh,
              const PDE& pde,
              const NumericalFlux& flux,
              const ModalVector& u,
              ModalVector& rhs);

#endif