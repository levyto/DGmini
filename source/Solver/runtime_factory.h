// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Definitions of runtime creation of objects 
//              (PDE, numerical flux, time integrator, time step controller)
//
// -----------------------------------------------------------------------------

#ifndef RUNTIME_FACTORY_H
#define RUNTIME_FACTORY_H

#include <memory>

#include "IO/input_config.h"
#include "PDE/pde.h"
#include "Spatial/NumericalFlux/numerical_flux.h"
#include "Temporal/TimeIntegrator/time_integrator.h"
#include "Temporal/TimeStepController/time_step_controller.h"

std::unique_ptr<PDE> createPDE(const InputConfig& config);
std::unique_ptr<NumericalFlux> createNumericalFlux(const InputConfig& config);
std::unique_ptr<TimeIntegrator> createTimeIntegrator(const InputConfig& config);
std::unique_ptr<TimeStepController> createTimeStepController(const InputConfig& config);

#endif
