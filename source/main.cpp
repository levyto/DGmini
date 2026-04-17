// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
//              A lightweight, self-contained C++ implementation of              
//                       1D discontinuous Galerkin solver                       
//       for experimentation and understanding of core DG building blocks       
//
// -----------------------------------------------------------------------------

#include "IO/config_reader.h"
#include "IO/input_config.h"
#include "IO/cli_parser.h"
#include "Solver/solver.h"

int main(int argc, char* argv[])
{
  const std::string config_file = getCmdLineArgsOrExit(argc, argv);

  ConfigReader reader;
  InputConfig config = reader.read(config_file);

  Solver solver(config);
  solver.run();

  return 0;
}

// int main()
// {
//   ConfigReader reader;
//   InputConfig config = reader.read("config.yaml");

//   Solver solver(config);

//   solver.run();

//   return 0;
// }
