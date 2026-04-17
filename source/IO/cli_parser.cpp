// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Command-Line Interface (CLI) parser
//
// -----------------------------------------------------------------------------

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

#include "IO/cli_parser.h"

// -----------------------------------------------------------------------------
// Description: Parses command-line arguments to extract the configuration 
//              file path. Expects the format: --config path/to/input.yaml
// -----------------------------------------------------------------------------
std::string parseCmdLineArgs(int argc, char* argv[])
{
  std::string config_file;
  
  for (int i = 1; i < argc; ++i)
  {
    const std::string arg = argv[i];

    if (arg == "--config")
    {
      if (i + 1 >= argc)
      {
        throw std::runtime_error("Missing value for --config");
      }

      config_file = argv[++i];
    }
    else if (arg == "--help")
    {
      printHelp();
      std::exit(0);
    }
    else
    {
      throw std::runtime_error("Unknown argument: " + arg);
    }
  }

  if (config_file.empty())
  {
    if (std::filesystem::exists("./config.yaml"))
    {
      config_file = "config.yaml";
    }
    else
    {
      throw std::runtime_error("Missing required argument: --config");
    }
  }  

  return config_file;
}

// -----------------------------------------------------------------------------
// Description: Wrapper function to handle exceptions from argument parsing and
//              print usage instructions if an error occurs
// -----------------------------------------------------------------------------
std::string getCmdLineArgsOrExit(int argc, char* argv[])
{
  try
  {
    return parseCmdLineArgs(argc, argv);
  }
  catch (const std::exception& e)
  {
    std::cerr << "[ERROR] " << e.what() << "\n\n";
    printUsage();
    std::exit(1);
  }
}

// -----------------------------------------------------------------------------
// Description: Prints detailed help information for the CLI
// -----------------------------------------------------------------------------
void printHelp()
{
  std::cout
    << "DGmini - a minimal 1D discontinuous Galerkin solver\n\n";
  printUsage();
  std::cout    
    << "Arguments:\n"
    << "  --config <file>   Path to YAML input file.\n"
    << "  --help            Print this help message and exit.\n\n"
    << "Default behavior:\n"
    << "  If no command-line argument is given, DGmini looks for\n"
    << "  a file named 'config.yaml' in the current working directory.\n\n"
    << "Examples:\n"
    << "  dgmini --config examples/advection.yaml\n"
    << "  dgmini\n";
}

// -----------------------------------------------------------------------------
// Description: Prints usage instructions for the CLI
// -----------------------------------------------------------------------------
void printUsage()
{
  std::cout 
    << "Usage:\n"
    << "  dgmini --config path/to/input.yaml\n"
    << "  dgmini\n\n";
}