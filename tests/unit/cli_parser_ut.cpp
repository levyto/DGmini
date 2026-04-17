// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Unittests for CLI parser
//
// -----------------------------------------------------------------------------

#include <filesystem>
#include <fstream>
#include <string>

#include "cli_parser_ut.h"

// -----------------------------------------------------------------------------
// Description: Helper functuons for UTs
// -----------------------------------------------------------------------------
namespace
{
  // ---------------------------------------------------------------------------
  // Description: Helper structure to temporarily change the current working 
  //              directory for a test, and automatically restore it when the 
  //              test is done
  // ---------------------------------------------------------------------------
  struct CurrentPathGuard
  {
    explicit CurrentPathGuard(const std::filesystem::path& path)
      : old_path(std::filesystem::current_path())
    {
      std::filesystem::current_path(path);
    }

    ~CurrentPathGuard()
    {
      std::filesystem::current_path(old_path);
    }

    std::filesystem::path old_path;
  };
}

// -----------------------------------------------------------------------------
// Description: Providing an explicit --config argument should return the 
//              provided path
// -----------------------------------------------------------------------------
void Test_CLIParser_explicitConfigArgument()
{
  char arg0[] = "dgmini";
  char arg1[] = "--config";
  char arg2[] = "examples/advection.yaml";
  char* argv[] = {arg0, arg1, arg2};

  const std::string config_file = parseCmdLineArgs(3, argv);

  Check(config_file == "examples/advection.yaml",
        "Explicit --config argument should be returned unchanged");
}

// -----------------------------------------------------------------------------
// Description: Providing --config without a value should cause the parser to 
//              throw an error
// -----------------------------------------------------------------------------
void Test_CLIParser_missingConfigValueThrows()
{
  char arg0[] = "dgmini";
  char arg1[] = "--config";
  char* argv[] = {arg0, arg1};

  bool did_throw = false;
  try
  {
    (void) parseCmdLineArgs(2, argv);
  }
  catch (const std::exception&)
  {
    did_throw = true;
  }

  Check(did_throw, "Missing value after --config should throw");
}

// -----------------------------------------------------------------------------
// Description: Providing an unknown command-line argument should cause the parser              
//              to throw an error
// -----------------------------------------------------------------------------
void Test_CLIParser_unknownArgumentThrows()
{
  char arg0[] = "dgmini";
  char arg1[] = "--unknown";
  char* argv[] = {arg0, arg1};

  bool did_throw = false;
  try
  {
    (void) parseCmdLineArgs(2, argv);
  }
  catch (const std::exception&)
  {
    did_throw = true;
  }

  Check(did_throw, "Unknown command-line argument should throw");
}

// -----------------------------------------------------------------------------
// Description: If no --config argument is provided, the parser should look for
//              a default config.yaml in the current working directory and return
//              its path if found
// -----------------------------------------------------------------------------
void Test_CLIParser_usesDefaultConfigInWorkingDirectory()
{
  const std::filesystem::path temp_dir = "output/unittest_cli_parser_default";
  std::filesystem::create_directories(temp_dir);

  {
    std::ofstream out(temp_dir / "config.yaml");
    Check(static_cast<bool>(out), "Failed to create temporary default config.yaml");
    out << "dummy: true\n";
  }

  char arg0[] = "dgmini";
  char* argv[] = {arg0};

  std::string config_file;
  {
    CurrentPathGuard guard(temp_dir);
    config_file = parseCmdLineArgs(1, argv);
  }

  Check(config_file == "config.yaml",
        "Parser should fall back to config.yaml in the current working directory");

  std::filesystem::remove_all(temp_dir);
}

// -----------------------------------------------------------------------------
// Description: If no --config argument is provided and no default config.yaml
//              exists in the working directory, the parser should throw an error
// -----------------------------------------------------------------------------
void Test_CLIParser_missingConfigThrowsWhenNoDefaultExists()
{
  const std::filesystem::path temp_dir = "output/unittest_cli_parser_missing";
  std::filesystem::create_directories(temp_dir);

  char arg0[] = "dgmini";
  char* argv[] = {arg0};

  bool did_throw = false;
  {
    CurrentPathGuard guard(temp_dir);

    try
    {
      (void) parseCmdLineArgs(1, argv);
    }
    catch (const std::exception&)
    {
      did_throw = true;
    }
  }

  Check(did_throw, "Missing --config should throw if config.yaml is not present");

  std::filesystem::remove_all(temp_dir);
}
