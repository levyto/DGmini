// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Class for Forward Euler time integrator
//
//              u^{n+1} = u^n + dt * R(u^n)
//  
//              where R is the residual of the spatial discretization
// -----------------------------------------------------------------------------

#ifndef FORWARD_EULER_H
#define FORWARD_EULER_H

#include <algorithm>

#include "Temporal/time_integrator.h"
#include "Spatial/residual.h"

class ForwardEuler : public TimeIntegrator
{
  public:
    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------
    ForwardEuler() = default;

    // -------------------------------------------------------------------------
    // Access
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Modification
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Operations
    // -------------------------------------------------------------------------
    void doTimeStep
    (
      const FESpace1D& fe,
      const Mesh1D& mesh,
      const PDE& pde,
      const NumericalFlux& flux,
      double dt,
      ModalVector& solution
    ) override
    {
      assert(is_initialized_);

      residual(fe, mesh, pde, flux, solution, rhs_);
      solution.axpy(dt, rhs_);
    }
  
  private:
    // -------------------------------------------------------------------------
    // Data
    // -------------------------------------------------------------------------
};

#endif