// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Definitions of inputs and readers
//
// -----------------------------------------------------------------------------

#ifndef INPUT_H
#define INPUT_H

#include <memory>

#include "PDE/pde.h"
#include "Spatial/NumericalFlux/numerical_flux.h"
#include "Temporal/TimeIntegrator/time_integrator.h"
#include "Temporal/TimeStepController/time_step_controller.h"

std::unique_ptr<PDE> createPDE(const std::string& name);
std::unique_ptr<NumericalFlux> createNumericalFlux(const std::string& name);
std::unique_ptr<TimeIntegrator> createTimeIntegrator(const std::string& name);
std::unique_ptr<TimeStepController> createTimeStepController(const std::string& name, double value);

#endif