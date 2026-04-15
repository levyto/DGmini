// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Abstract TimeStepController class, intended to be used as a base 
//              class for specific time step controllers (e.g. CFL-based, 
//              fixed time step, etc.) and to be extended
//
// -----------------------------------------------------------------------------

#ifndef TIME_STEP_CONTROLLER_H
#define TIME_STEP_CONTROLLER_H

#include "FEM/fespace1d.h"
#include "Mesh/mesh1d.h"
#include "PDE/pde.h"
#include "Spatial/modal_vector.h"

class TimeStepController
{
  public:
    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------
    virtual ~TimeStepController() = default;

    // -------------------------------------------------------------------------
    // Access
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Modification
    // -------------------------------------------------------------------------
    virtual double computeTimeStep
    (
      const FESpace1D& fe,
      const Mesh1D& mesh,
      const PDE& pde,
      const ModalVector& solution
    ) const = 0;

    // -------------------------------------------------------------------------
    // Operations
    // -------------------------------------------------------------------------

  private:
    // -------------------------------------------------------------------------
    // Data
    // -------------------------------------------------------------------------
};

#endif
