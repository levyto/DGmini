// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Class for 1D Buckley-Leverett PDE, i.e. 
//              
//              u_t + (f(u))_x = 0  
//              with  
//              f(u) = u^2 / (u^2 + M*(1-u)^2)
//
// -----------------------------------------------------------------------------

#ifndef BUCKLEY_LEVERETT1D_H
#define BUCKLEY_LEVERETT1D_H

#include <cmath>

#include "pde.h"

class BuckleyLeverett1D : public PDE
{
  public:
    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------
    BuckleyLeverett1D(double M) : M_(M) {}

    // -------------------------------------------------------------------------
    // Access
    // -------------------------------------------------------------------------
    double convectiveFlux(double state) const override
    {
      return state * state / (state * state + M_ * (1 - state) * (1 - state))  ;
    }

    // -------------------------------------------------------------------------
    // Modification
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Operations
    // -------------------------------------------------------------------------
    double maxConvectiveEigenvalue(double state) const override
    {
      double f_prime = 2 * state * M_ * (1 - state) * (1 - state) / 
                       std::pow(state * state + M_ * (1 - state) * (1 - state), 2);
      return std::abs(f_prime);
    }

  private:
    // -------------------------------------------------------------------------
    // Data
    // -------------------------------------------------------------------------
    double M_; // Viscosity ratio, M = mu_w / mu_o
};

#endif
