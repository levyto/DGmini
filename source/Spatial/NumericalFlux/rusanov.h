// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Class for Rusanov numerical flux
//
// -----------------------------------------------------------------------------

#ifndef RUSANOV_FLUX_H
#define RUSANOV_FLUX_H

#include <algorithm>

#include "Spatial/NumericalFlux/numerical_flux.h"

class RusanovFlux : public NumericalFlux
{
  public:
    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------
    RusanovFlux() = default;

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

      const double lambda =
        std::max(pde.maxConvectiveEigenvalue(uL),
                 pde.maxConvectiveEigenvalue(uR));

      return 0.5 * (fL + fR) - 0.5 * lambda * (uR - uL);
    }

  private:
    // -------------------------------------------------------------------------
    // Data
    // -------------------------------------------------------------------------
};

#endif