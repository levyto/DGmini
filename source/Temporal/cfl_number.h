// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Compute effective CFL number for a given time step
//
// -----------------------------------------------------------------------------

#ifndef CFL_NUMBER_H
#define CFL_NUMBER_H

#include "FEM/fespace1d.h"
#include "Mesh/mesh1d.h"
#include "PDE/pde.h"
#include "Spatial/modal_vector.h"

double computeEffectiveCFL
(
  const FESpace1D& fe,
  const Mesh1D& mesh,
  const PDE& pde,
  const ModalVector& solution,
  double dt
);

double evalMaxConvectiveEigenvalue
(
  const FESpace1D& fe,
  const Mesh1D& mesh,
  const PDE& pde,
  const ModalVector& solution
);

#endif