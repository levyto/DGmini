// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Definitions of inputs and readers
//
// -----------------------------------------------------------------------------

#include <stdexcept>

#include "IO/input.h"
#include "PDE/linear_advection1d.h"
#include "PDE/burgers1d.h"
#include "Spatial/lax_friedrichs_flux.h"
#include "Spatial/rusanov_flux.h"
#include "Temporal/forward_euler.h"
#include "Temporal/runge_kutta_2.h"
#include "Temporal/runge_kutta_3_ssp.h"
#include "Temporal/runge_kutta_4.h"
#include "Temporal/fixed_time_step.h"
#include "Temporal/cfl_time_step.h"
  
// -----------------------------------------------------------------------------
// Description: Create a PDE instance based on its name
// -----------------------------------------------------------------------------
std::unique_ptr<PDE> createPDE(const std::string& name)
{
  if (name == "linear_advection1d")
    return std::make_unique<LinearAdvection1D>(1.0);

  if (name == "burgers")
    return std::make_unique<Burgers1D>();

  throw std::runtime_error("Unknown PDE model: " + name);
}

// -----------------------------------------------------------------------------
// Description: Create a NumericalFlux instance based on its name
// -----------------------------------------------------------------------------
std::unique_ptr<NumericalFlux> createNumericalFlux(const std::string& name)
{
  if (name == "lax_friedrichs")
    return std::make_unique<LaxFriedrichsFlux>(1.0);

  if (name == "rusanov")
    return std::make_unique<RusanovFlux>();

  throw std::runtime_error("Unknown numerical flux: " + name);
}

// -----------------------------------------------------------------------------
// Description: Create a TimeIntegrator instance based on its name
// -----------------------------------------------------------------------------
std::unique_ptr<TimeIntegrator> createTimeIntegrator(const std::string& name)
{
  if (name == "forward_euler")
    return std::make_unique<ForwardEuler>();

  if (name == "runge_kutta_2")
    return std::make_unique<RungeKutta2>();

  if (name == "runge_kutta_3_ssp")
    return std::make_unique<RungeKutta3SSP>();

  if (name == "runge_kutta_4")
    return std::make_unique<RungeKutta4>();

  throw std::runtime_error("Unknown time integrator: " + name);
}

// -----------------------------------------------------------------------------
// Description: Create a TimeStepController instance based on its name
// -----------------------------------------------------------------------------
std::unique_ptr<TimeStepController> createTimeStepController(const std::string& name, double value)
{
  if (name == "fixed_time_step")
    return std::make_unique<FixedTimeStep>(value);

  if (name == "cfl_time_step")
    return std::make_unique<CFLTimeStep>(value);

  throw std::runtime_error("Unknown time step controller: " + name);
}