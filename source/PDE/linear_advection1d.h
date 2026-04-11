// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Class for 1D linear advection PDE, i.e. u_t + a * u_x = 0, 
//              where a is the advection speed.
//
// -----------------------------------------------------------------------------

#ifndef LINEAR_ADVECTION1D_H
#define LINEAR_ADVECTION1D_H

#include <cmath>

#include "pde.h"

class LinearAdvection1D : public PDE
{
  public:
    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------
    explicit LinearAdvection1D(double velocity) : velocity_(velocity) {}

    // -------------------------------------------------------------------------
    // Access
    // -------------------------------------------------------------------------
    double velocity() const { return velocity_; }

    double convectiveFlux(double state) const override
    {
      return velocity_ * state;
    }

    // -------------------------------------------------------------------------
    // Modification
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Operations
    // -------------------------------------------------------------------------
    double maxConvectiveEigenvalue(double state) const override
    {
      (void) state;
      return std::abs(velocity_);
    }

  private:
    // -------------------------------------------------------------------------
    // Data
    // -------------------------------------------------------------------------
    double velocity_ = 0.0;
};
#endif