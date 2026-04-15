// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Class for CFL-based time step controller, which computes the 
//              time step size based on the CFL condition
//
//              Assuming for DG of order p, the maximum stable time step size is 
//              given by
//
//              dt = CFL * h / ((2p+1) * max_lambda)
// -----------------------------------------------------------------------------

#ifndef CFL_TIME_STEP_H
#define CFL_TIME_STEP_H

#include "Temporal/TimeStepController/time_step_controller.h"
#include "Temporal/cfl_number.h"

class CFLTimeStep : public TimeStepController
{
  public:
    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------
    explicit CFLTimeStep(double cfl) : cfl_(cfl) {}
    virtual ~CFLTimeStep() = default;

    // -------------------------------------------------------------------------
    // Access
    // -------------------------------------------------------------------------
    double getCFL() const { return cfl_; }

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
      const double max_lambda = evalMaxConvectiveEigenvalue(fe, mesh, pde, solution);

      // Uniform mesh
      const double h = mesh.element(0).right() - mesh.element(0).left();

      return cfl_ * h / ((2.0 * fe.order() + 1.0) * max_lambda);
    }

    // -------------------------------------------------------------------------
    // Operations
    // -------------------------------------------------------------------------

  private:
    // -------------------------------------------------------------------------
    // Data
    // -------------------------------------------------------------------------
    double cfl_ = 0.0;
};

#endif
