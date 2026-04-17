// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Class for Explicit Strong Stability Preserving Runge-Kutta 3
//              (SSP RK3) time integrator
//              (3rd order, 3 stages, written in Shu-Osher form)
//
//                u^{1} = u^{n} + dt R(u^{n})
//                u^{2} = (3/4) u^{n} + (1/4) u^{1} + (1/4) dt R(u^{1})
//              u^{n+1} = (1/3) u^{n} + (2/3) u^{2} + (2/3) dt R(u^{2})
//
//              where R is the residual of the spatial discretization
//     
//              Butcher tableau:     0   | 0    0    0
//                                   1   | 1    0    0
//                                   1/2 | 1/4  1/4  0
//                                   -------------------
//                                       | 1/6  1/6  2/3
//
// -----------------------------------------------------------------------------

#ifndef RUNGE_KUTTA_3_SSP_H
#define RUNGE_KUTTA_3_SSP_H

#include <algorithm>

#include "Temporal/TimeIntegrator/time_integrator.h"
#include "Spatial/residual.h"

class RungeKutta3SSP : public TimeIntegrator
{
  public:
    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------
    RungeKutta3SSP() = default;

    // -------------------------------------------------------------------------
    // Access
    // -------------------------------------------------------------------------
    double recommendedCFL() const override { return 1.0; }

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
      const BoundaryConditions1D& bc,
      const double time,
      const double dt,
      ModalVector& solution
    ) override
    {
      assert(is_initialized_);

      solution_1_ = solution;

      // u^{n+1} = 1/3 * u^{n}
      solution.zero();
      solution.axpy(1.0 / 3.0, solution_1_);

      // R(u^{n})
      double t = time;
      residual(fe, mesh, pde, flux, bc, t, solution_1_, rhs_);

      // u^{1} = u^{n} + dt R(u^{n})
      solution_2_ = solution_1_;
      solution_2_.axpy(dt, rhs_);

      // R(u^{1})
      t = time + dt;
      residual(fe, mesh, pde, flux, bc, t, solution_2_, rhs_);

      // u^{2} = (3/4) u^{n} + (1/4) u^{1} + (1/4) dt R(u^{1})
      solution_2_.axpy(-3.0 / 4.0, solution_2_);
      solution_2_.axpy(3.0 / 4.0, solution_1_);
      solution_2_.axpy(1.0 / 4.0 * dt, rhs_);

      // u^{n+1} += 2/3 * u^{2}
      solution.axpy(2.0 / 3.0, solution_2_);

      // R(u^{2})
      t = time + 0.5 * dt;
      residual(fe, mesh, pde, flux, bc, t, solution_2_, rhs_);

      // u^{n+1} += 2/3 * dt R(u^{2})
      solution.axpy(2.0 / 3.0 * dt, rhs_);
    }
  
  private:
    // -------------------------------------------------------------------------
    // Data
    // -------------------------------------------------------------------------
};

#endif