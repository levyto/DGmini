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

#include "Temporal/time_step_controller.h"

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
      double max_lambda = 0.0;

      for (int e = 0; e < mesh.Ne(); ++e)
      {
        const double* u_e = solution.elementPtr(e);

        // Max over volume quadrature points
        for (int q = 0; q < fe.nip(); ++q)
        {
          double uq = 0.0;

          for (int j = 0; j < fe.DoFs(); ++j)
          {
            uq += u_e[j] * fe.phi(q, j);
          }

          max_lambda = std::max(max_lambda, pde.maxConvectiveEigenvalue(uq));
        }

        // Max over left and right boundaries
        double uL = 0.0;
        double uR = 0.0;

        for (int j = 0; j < fe.DoFs(); ++j)
        {
          uL += u_e[j] * fe.phiLeft(j);
          uR += u_e[j] * fe.phiRight(j);
        }

        max_lambda = std::max(max_lambda, pde.maxConvectiveEigenvalue(uL));
        max_lambda = std::max(max_lambda, pde.maxConvectiveEigenvalue(uR));
      }

      assert(max_lambda > 0.0);

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
