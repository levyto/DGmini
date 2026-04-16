// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for ConfigReader class
//
// -----------------------------------------------------------------------------

#include <fstream>
#include <filesystem>
#include <string>

#include "config_reader_ut.h"

// -----------------------------------------------------------------------------
// Description: Helper functions for tests
// -----------------------------------------------------------------------------
namespace
{
  // ---------------------------------------------------------------------------
  // Description: Write text content into a temporary YAML file
  // ---------------------------------------------------------------------------
  void writeTextFile(const std::string& filename, const std::string& content)
  {
    std::ofstream out(filename);
    Check(static_cast<bool>(out), "Failed to create temporary YAML file");
    out << content;
  }

  // ---------------------------------------------------------------------------
  // Description: Helper function to check if ConfigReader::read throws 
  //              an exception for a given YAML file 
  //  ---------------------------------------------------------------------------
  bool readThrows(const ConfigReader& reader, const std::string& filename)
  {
    try
    {
      (void) reader.read(filename);
    }
    catch (const std::exception&)
    {
      return true;
    }

    return false;
  }
  
  // ---------------------------------------------------------------------------
  // Description: A valid YAML config string for testing
  // ---------------------------------------------------------------------------
  std::string valid_yaml =
  "mesh:\n"
  "  x_left: 0.0\n"
  "  x_right: 1.0\n"
  "  n_elements: 10\n"
  "fem:\n"
  "  order: 2\n"
  "pde:\n"
  "  type: linear_advection\n"
  "  advection_speed: 1.0\n"
  "flux:\n"
  "  type: rusanov\n"
  "  alpha: 1.0\n"
  "time_integrator:\n"
  "  type: rk3_ssp\n"
  "time_step_controller:\n"
  "  type: fixed\n"
  "  dt: 0.01\n"
  "run:\n"
  "  final_time: 0.1\n"
  "initial_condition:\n"
  "  expression: \"sin(2*pi*x)\"\n"
  "output:\n"
  "  directory: output\n"
  "  prefix: sol\n"
  "  output_dt: 0.01\n"
  "  write_initial: true\n"
  "  write_final: true\n";
}

// -----------------------------------------------------------------------------
// Description: Incrementally build YAML config and verify required sections/keys
// -----------------------------------------------------------------------------
void Test_ConfigReader_incrementalRequiredKeys()
{
  const std::string filename = "output/unittest_invalid_types.yaml";
  std::filesystem::create_directories("output");

  ConfigReader reader;

  // 0) Empty file -> must fail
  writeTextFile(filename, "");
  Check(readThrows(reader, filename), "Empty YAML file should throw");

  // 1) Missing fem, pde, flux, ... -> must fail
  std::string yaml =
    "mesh:\n"
    "  x_left: 0.0\n"
    "  x_right: 1.0\n"
    "  n_elements: 10\n";

  writeTextFile(filename, yaml);
  Check(readThrows(reader, filename), "Config with only mesh section should throw");

  // 2) Add fem, still missing others -> must fail
  yaml +=
    "fem:\n"
    "  order: 2\n";

  writeTextFile(filename, yaml);
  Check(readThrows(reader, filename), "Config missing PDE/flux/time/output sections should throw");

  // 3) Add pde, still missing others -> must fail
  yaml +=
    "pde:\n"
    "  type: linear_advection\n"
    "  advection_speed: 1.0\n";

  writeTextFile(filename, yaml);
  Check(readThrows(reader, filename), "Config missing flux/time/output sections should throw");

  // 4) Add flux -> still must fail
  yaml +=
    "flux:\n"
    "  type: rusanov\n"
    "  alpha: 1.0\n";

  writeTextFile(filename, yaml);
  Check(readThrows(reader, filename), "Config missing time integrator/controller/run/output should throw");

  // 5) Add time integrator -> still must fail
  yaml +=
    "time_integrator:\n"
    "  type: rk3_ssp\n";

  writeTextFile(filename, yaml);
  Check(readThrows(reader, filename), "Config missing time step controller/run/output should throw");

  // 6) Add time step controller without required dt for fixed -> must fail
  yaml +=
    "time_step_controller:\n"
    "  type: fixed\n";

  writeTextFile(filename, yaml);
  Check(readThrows(reader, filename), "Fixed time step controller without dt should throw");

  // 7) Add dt, but still missing run/output/initial_condition -> must fail
  yaml +=
    "  dt: 0.001\n";

  writeTextFile(filename, yaml);
  Check(readThrows(reader, filename), "Config missing run/output/initial_condition should throw");

  // 8) Add run -> still must fail
  yaml +=
    "run:\n"
    "  final_time: 0.1\n";

  writeTextFile(filename, yaml);
  Check(readThrows(reader, filename), "Config missing initial_condition/output should throw");

  // 9) Add initial_condition without expression -> must fail
  yaml +=
    "initial_condition:\n";

  writeTextFile(filename, yaml);
  Check(readThrows(reader, filename), "Initial condition without expression should throw");

  // 10) Add expression, still missing output -> must fail
  yaml +=
    "  expression: sin(2*pi*x)\n";

  writeTextFile(filename, yaml);
  Check(readThrows(reader, filename), "Config missing output section should throw");

  // 11) Add incomplete output -> must fail
  yaml +=
    "output:\n"
    "  directory: output\n";

  writeTextFile(filename, yaml);
  Check(readThrows(reader, filename), "Incomplete output section should throw");

  // 12) Complete output -> config should pass
  yaml +=
    "  prefix: solution\n"
    "  output_dt: 0.01\n"
    "  write_initial: true\n"
    "  write_final: true\n";

  writeTextFile(filename, yaml);
  Check(!readThrows(reader, filename), "Complete valid config should not throw");

  std::filesystem::remove(filename);
}

// -----------------------------------------------------------------------------
// Description: Test that invalid types/values in config sections cause
//              validation to fail
// -----------------------------------------------------------------------------  
void Test_ConfigReader_invalidEnumeratedValues()
{
  const std::string filename = "output/unittest_invalid_types.yaml";
  std::filesystem::create_directories("output");

  ConfigReader reader;

  auto rewrite = [&](const std::string& yaml)
  {
    writeTextFile(filename, yaml);
    return readThrows(reader, filename);
  };

  // Invalid PDE type
  {
    std::string bad = valid_yaml;
    const size_t pos = bad.find("linear_advection");
    bad.replace(pos, std::string("linear_advection").size(), "invalid_pde");

    Check(rewrite(bad), "Invalid PDE type should throw");
  }

  // Invalid flux type
  {
    std::string bad = valid_yaml;
    const size_t pos = bad.find("rusanov");
    bad.replace(pos, std::string("rusanov").size(), "invalid_flux");

    Check(rewrite(bad), "Invalid flux type should throw");
  }

  // Invalid time integrator
  {
    std::string bad = valid_yaml;
    const size_t pos = bad.find("rk3_ssp");
    bad.replace(pos, std::string("rk3_ssp").size(), "invalid_integrator");

    Check(rewrite(bad), "Invalid time integrator should throw");
  }

  // Invalid time step controller
  {
    std::string bad = valid_yaml;
    const size_t pos = bad.find("fixed");
    bad.replace(pos, std::string("fixed").size(), "invalid_controller");

    Check(rewrite(bad), "Invalid time step controller should throw");
  }

  std::filesystem::remove(filename);
}

// -----------------------------------------------------------------------------
// Description: Test that invalid numeric values in config sections cause
//              validation to fail
// -----------------------------------------------------------------------------
void Test_ConfigReader_invalidNumericValues()
{
  const std::string filename = "output/unittest_invalid_numeric.yaml";
  std::filesystem::create_directories("output");

  ConfigReader reader;

  auto rewrite = [&](const std::string& yaml)
  {
    writeTextFile(filename, yaml);
    return readThrows(reader, filename);
  };

  // x_right <= x_left
  {
    std::string bad = valid_yaml;
    const size_t pos = bad.find("x_right: 1.0");
    bad.replace(pos, std::string("x_right: 1.0").size(), "x_right: 0.0");
    Check(rewrite(bad), "mesh.x_right <= mesh.x_left should throw");
  }

  // Invalid n_elements
  {
    std::string bad = valid_yaml;
    const size_t pos = bad.find("n_elements: 10");
    bad.replace(pos, std::string("n_elements: 10").size(), "n_elements: 0");
    Check(rewrite(bad), "Non-positive mesh.n_elements should throw");
  }

  // Invalid FEM order
  {
    std::string bad = valid_yaml;
    const size_t pos = bad.find("order: 2");
    bad.replace(pos, std::string("order: 2").size(), "order: -1");
    Check(rewrite(bad), "Negative FEM order should throw");
  }

  // Invalid dt for fixed time step controller
  {
    std::string bad = valid_yaml;
    const size_t pos = bad.find("dt: 0.01");
    bad.replace(pos, std::string("dt: 0.01").size(), "dt: 0.0");
    Check(rewrite(bad), "Non-positive dt should throw");
  }

  // Invalid final_time
  {
    std::string bad = valid_yaml;
    const size_t pos = bad.find("final_time: 0.1");
    bad.replace(pos, std::string("final_time: 0.1").size(), "final_time: 0.0");
    Check(rewrite(bad), "Non-positive final_time should throw");
  }

  // Invalid output_dt
  {
    std::string bad = valid_yaml;
    const size_t pos = bad.find("output_dt: 0.01");
    bad.replace(pos, std::string("output_dt: 0.01").size(), "output_dt: 0.0");
    Check(rewrite(bad), "Non-positive output_dt should throw");
  }

  // Empty output directory
  {
    std::string bad = valid_yaml;
    const size_t pos = bad.find("directory: output");
    bad.replace(pos, std::string("directory: output").size(), "directory: \"\"");
    Check(rewrite(bad), "Empty output directory should throw");
  }

  // Empty output prefix
  {
    std::string bad = valid_yaml;
    const size_t pos = bad.find("prefix: sol");
    bad.replace(pos, std::string("prefix: sol").size(), "prefix: \"\"");
    Check(rewrite(bad), "Empty output prefix should throw");
  }

  // Empty initial condition expression
  {
    std::string bad = valid_yaml;
    const size_t pos = bad.find("expression: \"sin(2*pi*x)\"");
    bad.replace(pos, std::string("expression: \"sin(2*pi*x)\"").size(), "expression: \"\"");
    Check(rewrite(bad), "Empty initial condition expression should throw");
  }

  // Non-positive alpha for lax_friedrichs
  {
    std::string bad = valid_yaml;
    const size_t pos1 = bad.find("type: rusanov");
    bad.replace(pos1, std::string("type: rusanov").size(), "type: lax_friedrichs");

    const size_t pos2 = bad.find("alpha: 1.0");
    bad.replace(pos2, std::string("alpha: 1.0").size(), "alpha: 0.0");

    Check(rewrite(bad), "Non-positive alpha for lax_friedrichs should throw");
  }

  // Non-positive CFL for cfl time step controller
  {
    std::string bad = valid_yaml;
    const size_t pos1 = bad.find("type: fixed");
    bad.replace(pos1, std::string("type: fixed").size(), "type: cfl");

    const size_t pos2 = bad.find("dt: 0.01");
    bad.replace(pos2, std::string("dt: 0.01").size(), "cfl: 0.0");

    Check(rewrite(bad), "Non-positive CFL should throw for cfl controller");
  }

  std::filesystem::remove(filename);
}

void Test_ConfigReader_optionalDefaults()
{
  const std::string filename = "output/unittest_optional_defaults.yaml";
  std::filesystem::create_directories("output");

  const std::string yaml =
  "mesh:\n"
  "  x_left: 0.0\n"
  "  x_right: 1.0\n"
  "  n_elements: 10\n"
  "fem:\n"
  "  order: 2\n"
  "pde:\n"
  "  type: linear_advection\n"
  "flux:\n"
  "  type: rusanov\n"
  "time_integrator:\n"
  "  type: rk3_ssp\n"
  "time_step_controller:\n"
  "  type: cfl\n"
  "  cfl: 0.35\n"
  "run:\n"
  "  final_time: 0.1\n"
  "initial_condition:\n"
  "  expression: \"sin(2*pi*x)\"\n"
  "output:\n"
  "  directory: output\n"
  "  prefix: sol\n"
  "  output_dt: 0.01\n";


  writeTextFile(filename, yaml);

  ConfigReader reader;
  InputConfig config = reader.read(filename);

  Check(config.pde.advection_speed == 0.0,
        "Missing pde.advection_speed should default to 0.0");

  Check(config.flux.alpha == 1.0,
        "Missing flux.alpha should default to 1.0");

  Check(config.time_step_controller.dt == 0.0,
        "Missing fixed dt should default to 0.0");

  Check(config.time_step_controller.cfl == 0.35,
        "Explicit CFL value should be preserved");

  Check(config.output.write_initial == true,
        "Missing output.write_initial should default to true");

  Check(config.output.write_final == true,
        "Missing output.write_final should default to true");

  std::filesystem::remove(filename);
}
