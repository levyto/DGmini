// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Command-Line Interface (CLI) parser
//
// -----------------------------------------------------------------------------

#ifndef CLI_PARSER_H
#define CLI_PARSER_H

#include <string>

std::string getCmdLineArgsOrExit(int argc, char* argv[]);
std::string parseCmdLineArgs(int argc, char* argv[]);

void printHelp();
void printUsage();

#endif
