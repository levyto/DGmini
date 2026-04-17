// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for ExpressionFunction class
//
// -----------------------------------------------------------------------------

#ifndef EXPRESSION_FUNCTION_UT_H
#define EXPRESSION_FUNCTION_UT_H

#include "unittest.h"
#include "IO/expression_function.h"

// -----------------------------------------------------------------------------
// Description: ExpressionFunction UTs
// -----------------------------------------------------------------------------
void Test_ExpressionFunction_constant();
void Test_ExpressionFunction_identity();
void Test_ExpressionFunction_sine();
void Test_ExpressionFunction_gaussian();
void Test_ExpressionFunction_invalidExpressionThrows();
void Test_ExpressionFunction_unknownSymbolThrows();

// -----------------------------------------------------------------------------
// Description: ExpressionFunction UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_ExpressionFunction(TestRegistry& registry)
{
  registry.add("Test_ExpressionFunction_constant",
                Test_ExpressionFunction_constant);
  registry.add("Test_ExpressionFunction_identity",
                Test_ExpressionFunction_identity);
  registry.add("Test_ExpressionFunction_sine",
                Test_ExpressionFunction_sine);
  registry.add("Test_ExpressionFunction_gaussian",
                Test_ExpressionFunction_gaussian);
  registry.add("Test_ExpressionFunction_invalidExpressionThrows",
                Test_ExpressionFunction_invalidExpressionThrows);
  registry.add("Test_ExpressionFunction_unknownSymbolThrows",
                Test_ExpressionFunction_unknownSymbolThrows);
}

#endif
