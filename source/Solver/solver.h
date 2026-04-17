// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Solver class definition
//
// -----------------------------------------------------------------------------

#ifndef SOLVER_H
#define SOLVER_H

#include <cassert>
#include <memory>

#include "IO/input_config.h"
#include "IO/output.h"
#include "PDE/pde.h"
#include "Spatial/NumericalFlux/numerical_flux.h"
#include "Temporal/TimeIntegrator/time_integrator.h"
#include "Temporal/TimeStepController/time_step_controller.h"
#include "Mesh/boundary_conditions1d.h"
#include "Mesh/mesh1d.h"
#include "FEM/fespace1d.h"
#include "Spatial/modal_vector.h"

class Solver
{
  public:
    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------
    explicit Solver(const InputConfig& config);

    // -------------------------------------------------------------------------
    // Access
    // -------------------------------------------------------------------------
    void printSettings() const;

    // -------------------------------------------------------------------------
    // Modification
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Operations
    // -------------------------------------------------------------------------
    void run();

  private:
    // -------------------------------------------------------------------------
    // Access
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Modification
    // -------------------------------------------------------------------------
    void initializeObjects();
    void initializeSolution();

    // -------------------------------------------------------------------------
    // Operations
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Data
    // -------------------------------------------------------------------------
    InputConfig config_;

    BoundaryConditions1D bc_;

    Mesh1D mesh_;
    FESpace1D fe_;

    std::unique_ptr<PDE> pde_;
    std::unique_ptr<NumericalFlux> flux_;
    std::unique_ptr<TimeIntegrator> integrator_;
    std::unique_ptr<TimeStepController> dt_controller_;

    TimeSeriesWriter output_;
    ModalVector solution_;
};

#endif
