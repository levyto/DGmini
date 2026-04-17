// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Definitions of runtime creation of objects 
//              (PDE, numerical flux, time integrator, time step controller)
//
// -----------------------------------------------------------------------------

#include <stdexcept>

#include "Solver/runtime_factory.h"

#include "IO/input_config.h"

#include "PDE/linear_advection1d.h"
#include "PDE/burgers1d.h"

#include "Spatial/NumericalFlux/lax_friedrichs.h"
#include "Spatial/NumericalFlux/rusanov.h"

#include "Temporal/TimeIntegrator/forward_euler.h"
#include "Temporal/TimeIntegrator/runge_kutta_2.h"
#include "Temporal/TimeIntegrator/runge_kutta_3_ssp.h"
#include "Temporal/TimeIntegrator/runge_kutta_4.h"

#include "Temporal/TimeStepController/fixed_time_step.h"
#include "Temporal/TimeStepController/cfl_time_step.h"

// -----------------------------------------------------------------------------
// Description: Create a BoundaryCondition1D instance based on the input
// -----------------------------------------------------------------------------
BoundaryCondition createBC(const BCInput& input)
{
  BoundaryCondition bc;

  if (input.type == "periodic")
  {
    bc.type = BoundaryConditionType::Periodic;
    return bc;
  }

  if (input.type == "dirichlet")
  {
    bc.type = BoundaryConditionType::Dirichlet;
    bc.expression = std::make_unique<ExpressionFunction>(input.expression);
    return bc;
  }

  if (input.type == "outflow")
  {
    bc.type = BoundaryConditionType::Outflow;
    return bc;
  }

  throw std::runtime_error("Unknown boundary condition type: " + input.type);
}

// -----------------------------------------------------------------------------
// Description: Create a PDE instance based on its name
// -----------------------------------------------------------------------------
std::unique_ptr<PDE> createPDE(const InputConfig& config)
{
  const auto &pde = config.pde;

  if (pde.type == "linear_advection")
    return std::make_unique<LinearAdvection1D>(pde.advection_speed);

  if (pde.type == "burgers")
    return std::make_unique<Burgers1D>();

  throw std::runtime_error("Unknown PDE model: " + pde.type);
}

// -----------------------------------------------------------------------------
// Description: Create a NumericalFlux instance based on its name
// -----------------------------------------------------------------------------
std::unique_ptr<NumericalFlux> createNumericalFlux(const InputConfig& config)
{
  const auto &flux = config.flux;

  if (flux.type == "lax_friedrichs")
    return std::make_unique<LaxFriedrichsFlux>(flux.alpha);

  if (flux.type == "rusanov")
    return std::make_unique<RusanovFlux>();

  throw std::runtime_error("Unknown numerical flux: " + flux.type);
}

// -----------------------------------------------------------------------------
// Description: Create a TimeIntegrator instance based on its name
// -----------------------------------------------------------------------------
std::unique_ptr<TimeIntegrator> createTimeIntegrator(const InputConfig& config)
{
  const auto &ti = config.time_integrator;

  if (ti.type == "forward_euler")
    return std::make_unique<ForwardEuler>();

  if (ti.type == "rk2")
    return std::make_unique<RungeKutta2>();

  if (ti.type == "rk3_ssp")
    return std::make_unique<RungeKutta3SSP>();

  if (ti.type == "rk4")
    return std::make_unique<RungeKutta4>();

  throw std::runtime_error("Unknown time integrator: " + ti.type);
}

// -----------------------------------------------------------------------------
// Description: Create a TimeStepController instance based on its name
// -----------------------------------------------------------------------------
std::unique_ptr<TimeStepController> createTimeStepController(const InputConfig& config)
{
  const auto &tsc = config.time_step_controller;

  if (tsc.type == "fixed")
    return std::make_unique<FixedTimeStep>(tsc.dt);

  if (tsc.type == "cfl")
    return std::make_unique<CFLTimeStep>(tsc.cfl);

  throw std::runtime_error("Unknown time step controller: " + tsc.type);
}
