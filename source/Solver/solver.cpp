// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Solver class definition
//
// -----------------------------------------------------------------------------

#include <iostream>

#include "Solver/solver.h"
#include "IO/expression_function.h"
#include "Spatial/l2_projection.h"
#include "Temporal/TimeStepController/fixed_time_step.h"
#include "Temporal/cfl_number.h"
#include "Solver/runtime_factory.h"

using std::cout;

// -----------------------------------------------------------------------------
// Description: Solver constructor, initializes all components based on input 
//              config
// -----------------------------------------------------------------------------
Solver::Solver(const InputConfig& config)
  : config_(config),
    mesh_
    (
      config.mesh.x_left,
      config.mesh.x_right,
      config.mesh.n_elements
    ),
    fe_(config.fem.order),
    output_
    (
      config.output.directory,
      config.output.prefix,
      config.output.output_dt
    ),
    solution_(mesh_.Ne(), fe_.DoFs())
{
  initializeObjects();
  initializeSolution();
}

// -----------------------------------------------------------------------------
// Description: Initialize runtime objects 
//              (PDE, numerical flux, time integrator, time step controller)
// -----------------------------------------------------------------------------
void Solver::initializeObjects()
{
  pde_           = createPDE(config_);
  flux_          = createNumericalFlux(config_);
  integrator_    = createTimeIntegrator(config_);
  dt_controller_ = createTimeStepController(config_);

  if (    !(pde_)
       || !(flux_)
       || !(integrator_)
       || !(dt_controller_)
     )
  {
    throw std::runtime_error("Solver initialization failed: one or more runtime objects were not created.");
  }

  integrator_->initialize(mesh_, fe_);
}

// -----------------------------------------------------------------------------
// Description: Initialize solution vector, e.g. by projecting initial condition
// -----------------------------------------------------------------------------
void Solver::initializeSolution()
{
  solution_.zero();

  ExpressionFunction u0(config_.initial_condition.expression);

  for (int e = 0; e < mesh_.Ne(); ++e)
  {
    L2ProjectionOnElement(fe_, mesh_.element(e), u0, solution_.elementPtr(e));
  }
}

// -----------------------------------------------------------------------------
// Description: Print summary of solver settings
// -----------------------------------------------------------------------------
void Solver::printSettings() const
{
  cout << "\n--- DGmini Settings ----------------------------------\n\n";

  /* Mesh ------------------------------------------------------------------- */
  cout << "Mesh:\n";
  cout << "  domain: ["  << config_.mesh.x_left
       << ", "           << config_.mesh.x_right    << "]\n";
  cout << "  elements: " << config_.mesh.n_elements << "\n\n";

  /* FEM -------------------------------------------------------------------- */
  cout << "FEM:\n";
  cout << "  order: " << config_.fem.order 
       << " (total DoFs = " << mesh_.Ne() * fe_.DoFs() << ")\n\n";

  /* PDE -------------------------------------------------------------------- */
  cout << "PDE:\n";
  cout << "  type: " << config_.pde.type << "\n";

  if (config_.pde.type == "linear_advection")
  {
    cout << "  advection_speed: " << config_.pde.advection_speed << "\n";
  }

  cout << "\n";

  /* Flux ------------------------------------------------------------------- */
  cout << "Flux:\n";
  cout << "  type: " << config_.flux.type << "\n";
  if (config_.flux.type == "lax_friedrichs")
  {
    cout << "  alpha: " << config_.flux.alpha << "\n";
  }
  cout << "\n";

  /* Time integration -------------------------------------------------------- */
  cout << "Time integration:\n";
  cout << "  integrator: " << config_.time_integrator.type
       << " (recommended CFL = " << integrator_->recommendedCFL() << ")\n";
       
  cout << "  controller: " << config_.time_step_controller.type;

  if (config_.time_step_controller.type == "fixed")
  {
    cout << " (dt = "  << config_.time_step_controller.dt  << ")";
  }
  else if (config_.time_step_controller.type == "cfl")
  {
    cout << " (cfl = " << config_.time_step_controller.cfl << ")";
  }

  cout << "\n\n";

  /* Run -------------------------------------------------------------------- */
  cout << "Run:\n";
  cout << "  final_time: " << config_.run.final_time << "\n\n";

  /* Initial condition ------------------------------------------------------ */
  cout << "Initial condition:\n";
  cout << "  expression: " << config_.initial_condition.expression << "\n\n";

  /* Output ------------------------------------------------------------------ */
  cout << "Output:\n";
  cout << "  directory: "  << config_.output.directory << "\n";
  cout << "  prefix: "     << config_.output.prefix    << "\n";
  cout << "  output_dt: "  << config_.output.output_dt << "\n";

  if (config_.output.write_initial == true)
  {
    cout << "  write_initial: true \n";
  }
  if (config_.output.write_final == true)
  {
    cout << "  write_final: true \n";
  }

  cout << "\n------------------------------------------------------\n\n";
}

// -----------------------------------------------------------------------------
// Description: Main time-stepping loop
// -----------------------------------------------------------------------------
void Solver::run()
{
  printSettings();

  double time = 0.0;
  double dt = 0.0;
  int step = 0;

  if (config_.output.write_initial)
  {
    output_.write(mesh_, solution_, time);
  }


  while (time < config_.run.final_time)
  {
    dt = dt_controller_->computeTimeStep(fe_, mesh_, *pde_, solution_);
    dt = std::min(dt, config_.run.final_time - time);


    /* Check current CFL number --------------------------------------------- */
    const double effective_cfl = computeEffectiveCFL(fe_, mesh_, *pde_, solution_, dt);

    if (effective_cfl > integrator_->recommendedCFL())
    {
      std::cout << "\n\nWARNING: Current CFL = " << effective_cfl
                << ", recommended is "           << integrator_->recommendedCFL()
                << "\n";
    }
    /* ---------------------------------------------------------------------- */


    integrator_->doTimeStep(fe_, mesh_, *pde_, *flux_, dt, solution_);

    time += dt;
    step++;

    output_.write(mesh_, solution_, time);
  }


  if (config_.output.write_final)
  {
    output_.writeFinal(mesh_, solution_, time);
  }

  cout << "\n\nFinished at t = " << time
       << " after " << step << " steps.\n";
}