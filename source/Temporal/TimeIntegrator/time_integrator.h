// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Abstract TimeIntegrator class, intended to be used as a base 
//              class for specific time integrators (e.g. Forward Euler, SSPRK,
//              etc.) and to be extended
//
// -----------------------------------------------------------------------------

#ifndef TIME_INTEGRATOR_H
#define TIME_INTEGRATOR_H

#include "FEM/fespace1d.h"
#include "Mesh/mesh1d.h"
#include "PDE/pde.h"
#include "Spatial/modal_vector.h"
#include "Spatial/NumericalFlux/numerical_flux.h"

class TimeIntegrator
{
  public:
    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------
    TimeIntegrator() = default;
    virtual ~TimeIntegrator() = default;

    // -------------------------------------------------------------------------
    // Access
    // -------------------------------------------------------------------------
    virtual double recommendedCFL() const = 0; // Empirical value for stability

    // -------------------------------------------------------------------------
    // Modification
    // -------------------------------------------------------------------------
    virtual void initialize(const Mesh1D& mesh, const FESpace1D& fe)
    {
      solution_1_ = ModalVector(mesh.Ne(), fe.DoFs());
      solution_2_ = ModalVector(mesh.Ne(), fe.DoFs());
      rhs_        = ModalVector(mesh.Ne(), fe.DoFs());
      is_initialized_ = true;
    }

    // -------------------------------------------------------------------------
    // Operations
    // -------------------------------------------------------------------------
    virtual void doTimeStep
    (
      const FESpace1D& fe,
      const Mesh1D& mesh,
      const PDE& pde,
      const NumericalFlux& flux,
      double dt,
      ModalVector& solution
    ) = 0;

  protected:
    // -------------------------------------------------------------------------
    // Data
    // -------------------------------------------------------------------------
    ModalVector solution_1_;
    ModalVector solution_2_;
    ModalVector rhs_;
    bool is_initialized_ = false;

  private:
    // -------------------------------------------------------------------------
    // Data
    // -------------------------------------------------------------------------
};

#endif