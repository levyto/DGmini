// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: ConfigReader class declaration
//
// -----------------------------------------------------------------------------

#include <filesystem>
#include <stdexcept>
#include <string>
#include <yaml-cpp/yaml.h>

#include "IO/config_reader.h"

// -----------------------------------------------------------------------------
// Description: Helper functions for YAML parsing
// -----------------------------------------------------------------------------
namespace 
{
  // ---------------------------------------------------------------------------
  // Description: Helper function to get a required section from the YAML root 
  //              node, throwing an error if the section is missing
  // ---------------------------------------------------------------------------
  const YAML::Node getRequiredSection(const YAML::Node& root, const std::string& key)
  {
    if (!root[key])
      throw std::runtime_error("ERROR: YAML config file: Missing required section: " + key + "\n");
    return root[key];
  }

  // ---------------------------------------------------------------------------
  // Description: Helper function to get a required value from a YAML node, 
  //              throwing an error if the key is missing
  // ---------------------------------------------------------------------------
  template<typename T>
  T getRequired(const YAML::Node& node, const std::string& key, const std::string& path)
  {
    if (!node[key])
      throw std::runtime_error("ERROR: YAML config file: Missing required key: " + path + "\n");
    return node[key].as<T>();
  }

  // ---------------------------------------------------------------------------
  // Description: Helper function to get an optional value from a YAML node, 
  //              returning a default value if the key is missing
  // ---------------------------------------------------------------------------
  template<typename T>
  T getOptional(const YAML::Node& node, const std::string& key, const T& default_value)
  {
    return node[key] ? node[key].as<T>() : default_value;
  }
}

// -----------------------------------------------------------------------------
// Description: Input YAML config parser that reads the config file and 
//              populates the InputConfig struct
// -----------------------------------------------------------------------------
InputConfig ConfigReader::read(const std::string& filename) const
{
  InputConfig config;

  if (!std::filesystem::exists(filename))
  {
    throw std::runtime_error(
      "\n\nERROR: "
      "Config file does not exist: " + filename + "\n");
  }

  YAML::Node root;
  try
  {
    root = YAML::LoadFile(filename);
  }
  catch (const YAML::Exception& e)
  {
    throw std::runtime_error(
      "\n\nERROR: "
      "Error parsing YAML file '" + filename + "': " + e.what() + "\n");
  }

  /* mesh ------------------------------------------------------------------- */
  const YAML::Node mesh = getRequiredSection(root, "mesh");

  config.mesh.x_left     = getRequired<double>(mesh, "x_left",  "mesh.x_left");
  config.mesh.x_right    = getRequired<double>(mesh, "x_right", "mesh.x_right");
  config.mesh.n_elements = getRequired<int>(mesh, "n_elements", "mesh.n_elements");

  /* fem -------------------------------------------------------------------- */
  const YAML::Node fem = getRequiredSection(root, "fem");

  config.fem.order = getRequired<int>(fem, "order", "fem.order");

  /* pde -------------------------------------------------------------------- */
  const YAML::Node pde = getRequiredSection(root, "pde");

  config.pde.type = getRequired<std::string>(pde, "type", "pde.type");
  config.pde.advection_speed = getOptional<double>(pde, "advection_speed", 0.0);

  /* flux ------------------------------------------------------------------- */
  const YAML::Node flux = getRequiredSection(root, "flux");

  config.flux.type  = getRequired<std::string>(flux, "type", "flux.type");
  config.flux.alpha = getOptional<double>(flux, "alpha", 1.0);

  /* time integrator -------------------------------------------------------- */
  const YAML::Node ti = getRequiredSection(root, "time_integrator");

  config.time_integrator.type = getRequired<std::string>(ti, "type", "time_integrator.type");

  /* time step controller --------------------------------------------------- */
  const YAML::Node tsc = getRequiredSection(root, "time_step_controller");

  config.time_step_controller.type = getRequired<std::string>(tsc, "type", "time_step_controller.type");
  config.time_step_controller.dt   = getOptional<double>(tsc, "dt", 0.0);
  config.time_step_controller.cfl  = getOptional<double>(tsc, "cfl", 0.0);

  /* run -------------------------------------------------------------------- */
  const YAML::Node run = getRequiredSection(root, "run");

  config.run.final_time = getRequired<double>(run, "final_time", "run.final_time");

  /* initial condition ------------------------------------------------------ */
  const YAML::Node ic = getRequiredSection(root, "initial_condition");

  config.initial_condition.expression = getRequired<std::string>(ic, "expression", "initial_condition.expression");

  /* output ----------------------------------------------------------------- */
  const YAML::Node output = getRequiredSection(root, "output");

  config.output.directory     = getRequired<std::string>(output, "directory", "output.directory");
  config.output.prefix        = getRequired<std::string>(output, "prefix", "output.prefix");
  config.output.output_dt     = getRequired<double>(output, "output_dt", "output.output_dt");
  config.output.write_initial = getOptional<bool>(output, "write_initial", true);
  config.output.write_final   = getOptional<bool>(output, "write_final", true);

  /* validation ------------------------------------------------------------- */
  this->validate(config);

  return config;
}

// -----------------------------------------------------------------------------
// Description: Validate the input configuration
// -----------------------------------------------------------------------------
void ConfigReader::validate(const InputConfig& config) const
{
  /* mesh ------------------------------------------------------------------- */
  if (config.mesh.x_right <= config.mesh.x_left)
  {
    throw std::runtime_error(
      "\n\nERROR: YAML config file: Invalid 'mesh: x_right' entry, must be greater than 'mesh: x_left'.");
  }

  if (config.mesh.n_elements <= 0)
  {
    throw std::runtime_error(
      "\n\nERROR: YAML config file: Invalid 'mesh: n_elements' entry, must be positive.");
  }

  /* fem -------------------------------------------------------------------- */
  if (config.fem.order < 0)
  {
    throw std::runtime_error(
      "\n\nERROR: YAML config file: Invalid 'fem: order' entry, must be non-negative.");
  }

  /* pde -------------------------------------------------------------------- */
  if (    (config.pde.type != "linear_advection")
       && (config.pde.type != "burgers")
     )
  {
    throw std::runtime_error(
      "\n\nERROR: YAML config file: Invalid 'pde: type' entry, supported values are:"
      "\n\t'linear_advection'"
      "\n\t'burgers'"
      "\n"
    );
  }

  /* flux ------------------------------------------------------------------- */
  if (    (config.flux.type != "rusanov")
       && (config.flux.type != "lax_friedrichs")
     )
  {
    throw std::runtime_error(
      "\n\nERROR: YAML config file: Invalid 'flux: type' entry, supported values are:"
      "\n\t'lax_friedrichs'"
      "\n\t'rusanov'"
      "\n"
    );
  }

  if (    (config.flux.type == "lax_friedrichs")
       && (config.flux.alpha <= 0.0)
     )
  {
    throw std::runtime_error(
      "\n\nERROR: YAML config file: Invalid 'flux: alpha' entry, must be positive for 'lax_friedrichs' flux.");
  }

  /* time integrator -------------------------------------------------------- */
  if (    (config.time_integrator.type != "forward_euler")
       && (config.time_integrator.type != "rk2")
       && (config.time_integrator.type != "rk3_ssp")
       && (config.time_integrator.type != "rk4")
     )
  {
    throw std::runtime_error(
      "\n\nERROR: YAML config file: Invalid 'time_integrator: type' entry, supported values are:"
      "\n\t'forward_euler'"
      "\n\t'rk2'"
      "\n\t'rk3_ssp'"
      "\n\t'rk4'"
      "\n"
  );
  }

  /* time step controller --------------------------------------------------- */
  if (    (config.time_step_controller.type != "fixed")
       && (config.time_step_controller.type != "cfl")
     )
  {
    throw std::runtime_error(
      "\n\nERROR: YAML config file: Invalid 'time_step_controller: type' entry, supported values are:"
      "\n\t'fixed'"
      "\n\t'cfl'"
      "\n"
    );
  }

  if (config.time_step_controller.type == "fixed")
  {
    if (config.time_step_controller.dt <= 0.0)
    {
      throw std::runtime_error(
        "\n\nERROR: YAML config file: Invalid 'time_step_controller: dt' entry, must be positive for 'fixed' type.");
    }
  }

  if (config.time_step_controller.type == "cfl")
  {
    if (config.time_step_controller.cfl <= 0.0)
    {
      throw std::runtime_error(
        "\n\nERROR: YAML config file: Invalid 'time_step_controller: cfl' entry, must be positive for 'cfl' type.\n");
    }
  }

  /* run -------------------------------------------------------------------- */
  if (config.run.final_time <= 0.0)
  {
    throw std::runtime_error(
      "\n\nERROR: YAML config file: Invalid 'run: final_time' entry, must be positive.");
  }

  /* initial condition ------------------------------------------------------ */
  if (config.initial_condition.expression.empty())
  {
    throw std::runtime_error(
      "\n\nERROR: YAML config file: Invalid 'initial_condition: expression' entry, must not be empty.");
  }

  /* output ----------------------------------------------------------------- */
  if (config.output.directory.empty())
  {
    throw std::runtime_error(
      "\n\nERROR: YAML config file: Invalid 'output: directory' entry, must not be empty.");
  }

  if (config.output.prefix.empty())
  {
    throw std::runtime_error(
      "\n\nERROR: YAML config file: Invalid 'output: prefix' entry, must not be empty.");
  }

  if (config.output.output_dt <= 0.0)
  {
    throw std::runtime_error(
      "\n\nERROR: YAML config file: Invalid 'output: output_dt' entry, must be positive.");
  }
}