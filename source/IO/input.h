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
#include "Spatial/numerical_flux.h"
#include "Temporal/time_integrator.h"

std::unique_ptr<PDE> createPDE(const std::string& name);
std::unique_ptr<NumericalFlux> createNumericalFlux(const std::string& name);
std::unique_ptr<TimeIntegrator> createTimeIntegrator(const std::string& name);

#endif