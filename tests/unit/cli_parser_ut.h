// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Unittests for CLI parser
//
// -----------------------------------------------------------------------------

#ifndef CLI_PARSER_UT_H
#define CLI_PARSER_UT_H

#include "unittest.h"
#include "IO/cli_parser.h"

// -----------------------------------------------------------------------------
// Description: CLI parser UTs
// -----------------------------------------------------------------------------
void Test_CLIParser_explicitConfigArgument();
void Test_CLIParser_missingConfigValueThrows();
void Test_CLIParser_unknownArgumentThrows();
void Test_CLIParser_usesDefaultConfigInWorkingDirectory();
void Test_CLIParser_missingConfigThrowsWhenNoDefaultExists();

// -----------------------------------------------------------------------------
// Description: CLI parser UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_CLIParser(TestRegistry& registry)
{
  registry.add("Test_CLIParser_explicitConfigArgument",
                Test_CLIParser_explicitConfigArgument);
  registry.add("Test_CLIParser_missingConfigValueThrows",
                Test_CLIParser_missingConfigValueThrows);
  registry.add("Test_CLIParser_unknownArgumentThrows",
                Test_CLIParser_unknownArgumentThrows);
  registry.add("Test_CLIParser_usesDefaultConfigInWorkingDirectory",
                Test_CLIParser_usesDefaultConfigInWorkingDirectory);
  registry.add("Test_CLIParser_missingConfigThrowsWhenNoDefaultExists",
                Test_CLIParser_missingConfigThrowsWhenNoDefaultExists);
}

#endif
