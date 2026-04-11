// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Class for Lax-Friedrichs numerical flux
//
// -----------------------------------------------------------------------------

#ifndef LAX_FRIEDRICHS_FLUX_H
#define LAX_FRIEDRICHS_FLUX_H

#include <algorithm>

#include "Spatial/numerical_flux.h"

class LaxFriedrichsFlux : public NumericalFlux
{
  public:
    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------
    explicit LaxFriedrichsFlux(double alpha) : alpha_(alpha) {}

    // -------------------------------------------------------------------------
    // Access
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Modification
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Operations
    // -------------------------------------------------------------------------
    double evaluate(const PDE& pde, double uL, double uR) const override
    {
      const double fL = pde.convectiveFlux(uL);
      const double fR = pde.convectiveFlux(uR);

      return 0.5 * (fL + fR) - 0.5 * alpha_ * (uR - uL);
    }

  private:
    // -------------------------------------------------------------------------
    // Data
    // -------------------------------------------------------------------------
    double alpha_ = 0.0;
};

#endif
