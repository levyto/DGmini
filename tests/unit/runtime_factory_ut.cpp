// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Unittests for runtime factory
//
// -----------------------------------------------------------------------------

#include <memory>
#include <stdexcept>

#include "runtime_factory_ut.h"
#include "PDE/burgers1d.h"
#include "PDE/linear_advection1d.h"
#include "Spatial/NumericalFlux/lax_friedrichs.h"
#include "Spatial/NumericalFlux/rusanov.h"
#include "Temporal/TimeIntegrator/forward_euler.h"
#include "Temporal/TimeIntegrator/runge_kutta_4.h"
#include "Temporal/TimeStepController/cfl_time_step.h"
#include "Temporal/TimeStepController/fixed_time_step.h"

// -----------------------------------------------------------------------------
// Description: Helper functions for UTs
// -----------------------------------------------------------------------------
namespace
{
  // ---------------------------------------------------------------------------
  // Description: Create a valid InputConfig with default values for all fields,
  //              which can be modified by tests as needed
  // ---------------------------------------------------------------------------
  InputConfig makeValidConfig()
  {
    InputConfig config;

    config.boundary_conditions.left.type = "periodic";
    config.boundary_conditions.right.type = "periodic";

    config.mesh.x_left = 0.0;
    config.mesh.x_right = 1.0;
    config.mesh.n_elements = 4;

    config.fem.order = 1;

    config.pde.type = "linear_advection";
    config.pde.advection_speed = 2.5;

    config.flux.type = "rusanov";
    config.flux.alpha = 1.7;

    config.time_integrator.type = "rk4";

    config.time_step_controller.type = "fixed";
    config.time_step_controller.dt = 0.125;
    config.time_step_controller.cfl = 0.4;

    config.run.final_time = 1.0;
    config.initial_condition.expression = "sin(2*pi*x)";

    config.output.directory = "output";
    config.output.prefix = "sol";
    config.output.output_dt = 0.1;
    config.output.write_initial = true;
    config.output.write_final = true;

    return config;
  }

  // ---------------------------------------------------------------------------
  // Description: Helper function to check if a factory function throws an 
  //              exception for a given input
  // ---------------------------------------------------------------------------
  template <typename Factory>
  bool factoryThrows(Factory&& factory)
  {
    try
    {
      (void) factory();
    }
    catch (const std::exception&)
    {
      return true;
    }

    return false;
  }
}

// -----------------------------------------------------------------------------
// Description: Test that createBC correctly creates BoundaryCondition objects
//              based on the provided BCInput, and that the created objects
//              have the expected properties
// -----------------------------------------------------------------------------
void Test_RuntimeFactory_createBoundaryCondition()
{
  {
    BCInput input;
    input.type = "periodic";

    const BoundaryCondition bc = createBC(input);
    Check(bc.type == BoundaryConditionType::Periodic,
          "Periodic BC should be created with Periodic type");
    Check(!bc.expression,
          "Periodic BC should not allocate an expression");
  }

  {
    BCInput input;
    input.type = "outflow";

    const BoundaryCondition bc = createBC(input);
    Check(bc.type == BoundaryConditionType::Outflow,
          "Outflow BC should be created with Outflow type");
    Check(!bc.expression,
          "Outflow BC should not allocate an expression");
  }

  {
    BCInput input;
    input.type = "dirichlet";
    input.expression = "x + 2*t";

    const BoundaryCondition bc = createBC(input);
    Check(bc.type == BoundaryConditionType::Dirichlet,
          "Dirichlet BC should be created with Dirichlet type");
    Check(static_cast<bool>(bc.expression),
          "Dirichlet BC should allocate an expression");
    CheckEqual((*(bc.expression))(1.0, 0.5), 2.0, 1e-14,
               "Dirichlet BC expression should be evaluable");
  }
}

// -----------------------------------------------------------------------------
// Description: Test that factory functions createPDE, createNumericalFlux,
//              createTimeIntegrator, and createTimeStepController correctly create
//              objects based on the provided InputConfig, and that the created
//              objects have the expected properties
// -----------------------------------------------------------------------------
void Test_RuntimeFactory_createObjects()
{
  // Test PDE creation
  {
    InputConfig config = makeValidConfig();

    std::unique_ptr<PDE> pde = createPDE(config);
    auto* advection = dynamic_cast<LinearAdvection1D*>(pde.get());

    Check(advection != nullptr,
          "linear_advection should create LinearAdvection1D");
    CheckEqual(advection->velocity(), 2.5, 1e-14,
               "LinearAdvection1D should preserve configured velocity");
  }

  // Test PDE creation with different type
  {
    InputConfig config = makeValidConfig();
    config.pde.type = "burgers";

    std::unique_ptr<PDE> pde = createPDE(config);
    Check(dynamic_cast<Burgers1D*>(pde.get()) != nullptr,
          "burgers should create Burgers1D");
  }

  // Test numerical flux creation
  {
    InputConfig config = makeValidConfig();
    config.flux.type = "lax_friedrichs";
    config.flux.alpha = 3.0;

    std::unique_ptr<NumericalFlux> flux = createNumericalFlux(config);
    auto* lf = dynamic_cast<LaxFriedrichsFlux*>(flux.get());

    Check(lf != nullptr,
          "lax_friedrichs should create LaxFriedrichsFlux");
    CheckEqual(lf->evaluate(LinearAdvection1D(2.0), 1.0, 4.0), 0.5, 1e-14,
               "Lax-Friedrichs flux should preserve configured alpha");
  }

  // Test numerical flux creation with different type
  {
    InputConfig config = makeValidConfig();

    std::unique_ptr<NumericalFlux> flux = createNumericalFlux(config);
    Check(dynamic_cast<RusanovFlux*>(flux.get()) != nullptr,
          "rusanov should create RusanovFlux");
  }

  // Test time integrator creation
  {
    InputConfig config = makeValidConfig();
    config.time_integrator.type = "forward_euler";

    std::unique_ptr<TimeIntegrator> integrator = createTimeIntegrator(config);
    Check(dynamic_cast<ForwardEuler*>(integrator.get()) != nullptr,
          "forward_euler should create ForwardEuler");
  }

  // Test time integrator creation with different type
  {
    InputConfig config = makeValidConfig();

    std::unique_ptr<TimeIntegrator> integrator = createTimeIntegrator(config);
    Check(dynamic_cast<RungeKutta4*>(integrator.get()) != nullptr,
          "rk4 should create RungeKutta4");
  }

  // Test time step controller creation
  {
    InputConfig config = makeValidConfig();

    std::unique_ptr<TimeStepController> controller = createTimeStepController(config);
    auto* fixed = dynamic_cast<FixedTimeStep*>(controller.get());

    Check(fixed != nullptr,
          "fixed should create FixedTimeStep");
    CheckEqual(fixed->getTimeStep(), 0.125, 1e-14,
               "FixedTimeStep should preserve configured dt");
  }

  // Test time step controller creation with different type
  {
    InputConfig config = makeValidConfig();
    config.time_step_controller.type = "cfl";
    config.time_step_controller.cfl = 0.35;

    std::unique_ptr<TimeStepController> controller = createTimeStepController(config);
    auto* cfl = dynamic_cast<CFLTimeStep*>(controller.get());

    Check(cfl != nullptr,
          "cfl should create CFLTimeStep");
    CheckEqual(cfl->getCFL(), 0.35, 1e-14,
               "CFLTimeStep should preserve configured CFL number");
  }
}

// -----------------------------------------------------------------------------
// Description: Test that factory functions throw exceptions when given invalid
//              inputs, such as unknown types in the InputConfig or invalid
//              parameters for specific types
// -----------------------------------------------------------------------------
void Test_RuntimeFactory_invalidInputsThrow()
{
  // Invalid BC type
  {
    BCInput input;
    input.type = "invalid_bc";
    Check(factoryThrows([&]() { return createBC(input); }),
          "Unknown BC type should throw");
  }

  // Invalid PDE type
  {
    InputConfig config = makeValidConfig();
    config.pde.type = "invalid_pde";
    Check(factoryThrows([&]() { return createPDE(config); }),
          "Unknown PDE type should throw");
  }

  // Invalid flux type
  {
    InputConfig config = makeValidConfig();
    config.flux.type = "invalid_flux";
    Check(factoryThrows([&]() { return createNumericalFlux(config); }),
          "Unknown numerical flux type should throw");
  }

  // Invalid time integrator type
  {
    InputConfig config = makeValidConfig();
    config.time_integrator.type = "invalid_integrator";
    Check(factoryThrows([&]() { return createTimeIntegrator(config); }),
          "Unknown time integrator type should throw");
  }

  // Invalid time step controller type
  {
    InputConfig config = makeValidConfig();
    config.time_step_controller.type = "invalid_controller";
    Check(factoryThrows([&]() { return createTimeStepController(config); }),
          "Unknown time step controller type should throw");
  }
}
