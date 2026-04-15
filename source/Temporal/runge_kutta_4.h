// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Class for Explicit Runge-Kutta 4 (RK4) time integrator
//              (4th order, 4 stages, written in modified Euler form)
//
//              Stage 1: u^{1} = u^{n}
//              Stage 2: u^{2} = u^{n} + 0.5 * dt * R(u^{1})
//              Stage 3: u^{3} = u^{n} + 0.5 * dt * R(u^{2})
//              Stage 4: u^{4} = u^{n} + dt * R(u^{3})
//              Final: u^{n+1} = u^{n} + (dt/6) * (R(u^{1}) + 2*R(u^{2}) + 2*R(u^{3}) + R(u^{4}))
//
//              where R is the residual of the spatial discretization
//     
//              Butcher tableau:     0   | 0    0    0    0
//                                   1/2 | 1/2  0    0    0
//                                   1/2 | 0    1/2  0    0
//                                   1   | 0    0    1    0
//                                   ------------------------
//                                       | 1/6  1/3  1/3  1/6
//
// -----------------------------------------------------------------------------

#ifndef RUNGE_KUTTA_4_H
#define RUNGE_KUTTA_4_H

#include <algorithm>

#include "Temporal/time_integrator.h"
#include "Spatial/residual.h"

class RungeKutta4 : public TimeIntegrator
{
  public:
    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------
    RungeKutta4() = default;

    // -------------------------------------------------------------------------
    // Access
    // -------------------------------------------------------------------------
    double recommendedCFL() const override { return 1.5; }

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

      // -----------------------------------------------------------------------
      // R(u^{1})
      residual(fe, mesh, pde, flux, solution_1_, rhs_);

      // u^{n+1} += 1/6*dt * R(u^{1})
      solution.axpy(1.0 / 6.0 * dt, rhs_); 

      // u^{2}
      solution_2_ = solution_1_;
      solution_2_.axpy(0.5 * dt, rhs_); 
     
      // -----------------------------------------------------------------------
      // R(u^{2})
      residual(fe, mesh, pde, flux, solution_2_, rhs_); 

      // u^{n+1} += 1/3*dt * R(u^{2})
      solution.axpy(1.0 / 3.0 * dt, rhs_); 

      // u^{3}
      solution_2_ = solution_1_;
      solution_2_.axpy(0.5 * dt, rhs_); 

      // -----------------------------------------------------------------------
      // R(u^{3})
      residual(fe, mesh, pde, flux, solution_2_, rhs_); 

      // u^{n+1} += 1/3*dt * R(u^{3})
      solution.axpy(1.0 / 3.0 * dt, rhs_); 
      
      // u^{4}
      solution_2_ = solution_1_;
      solution_2_.axpy(dt, rhs_);

      // ----------------------------------------------------------------------- 
      // R(u^{4})
      residual(fe, mesh, pde, flux, solution_2_, rhs_); 

      // u^{n+1} += 1/6*dt * R(u^{4})
      solution.axpy(1.0 / 6.0 * dt, rhs_); 
    }
  
  private:
    // -------------------------------------------------------------------------
    // Data
    // -------------------------------------------------------------------------
};

#endif