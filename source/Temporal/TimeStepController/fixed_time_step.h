// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Class for fixed time step controller, which always returns 
//              a user-specified time step size
//
// -----------------------------------------------------------------------------

#ifndef FIXED_TIME_STEP_H
#define FIXED_TIME_STEP_H

#include "Temporal/TimeStepController/time_step_controller.h"

class FixedTimeStep : public TimeStepController
{
  public:
    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------
    explicit FixedTimeStep(double dt) : dt_(dt) {}
    virtual ~FixedTimeStep() = default;

    // -------------------------------------------------------------------------
    // Access
    // -------------------------------------------------------------------------
    double getTimeStep() const { return dt_; }

    // -------------------------------------------------------------------------
    // Modification
    // -------------------------------------------------------------------------
    virtual double computeTimeStep
    (
      const FESpace1D& fe,
      const Mesh1D& mesh,
      const PDE& pde,
      const ModalVector& solution
    ) const override
    {
      return dt_;
    }

    // -------------------------------------------------------------------------
    // Operations
    // -------------------------------------------------------------------------

  private:
    // -------------------------------------------------------------------------
    // Data
    // -------------------------------------------------------------------------
    double dt_ = 0.0;
};

#endif
