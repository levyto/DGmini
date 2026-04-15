// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Class for Explicit Runge-Kutta 2 (RK2) time integrator
//              (2nd order, 2 stages, written in modified Euler form)
//
//              Stage 1: u^{1} = u^{n}
//              Stage 2: u^{2} = u^{n} + 0.5 * dt * R(u^{1})
//              Final: u^{n+1} = u^{n} + dt * R(u^{2})     
//
//              where R is the residual of the spatial discretization
//     
//              Butcher tableau:     0   | 0    0
//                                   1/2 | 1/2  0
//                                   ------------
//                                       | 0    1
//
// -----------------------------------------------------------------------------

#ifndef RUNGE_KUTTA_2_H
#define RUNGE_KUTTA_2_H

#include <algorithm>

#include "Temporal/TimeIntegrator/time_integrator.h"
#include "Spatial/residual.h"

class RungeKutta2 : public TimeIntegrator
{
  public:
    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------
    RungeKutta2() = default;

    // -------------------------------------------------------------------------
    // Access
    // -------------------------------------------------------------------------
    double recommendedCFL() const override { return 0.5; }
    
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

      solution_1_ = solution;

      // R(u^{1})
      residual(fe, mesh, pde, flux, solution_1_, rhs_);

      // u^{2} = u^{n} + 0.5 * dt * R(u^{1})
      solution_2_ = solution;
      solution_2_.axpy(0.5 * dt, rhs_);

      // R(u^{2})
      residual(fe, mesh, pde, flux, solution_2_, rhs_);

      // u^{n+1} = u^{n} + dt * R(u^{2})
      solution.axpy(dt, rhs_);
    }
  
  private:
    // -------------------------------------------------------------------------
    // Data
    // -------------------------------------------------------------------------
};

#endif