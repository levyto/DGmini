// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for ExpressionFunction class
//
// -----------------------------------------------------------------------------

#include <string>

#include "expression_function_ut.h"

// -----------------------------------------------------------------------------
// Description: Helper functions for ExpressionFunction UTs
// -----------------------------------------------------------------------------
namespace
{
  bool ExpressionFunctionEvalThrows(const std::string& expression)
  {
    try
    {
      ExpressionFunction f(expression);
      (void) f(0.5);
    }
    catch (const std::exception&)
    {
      return true;
    }

    return false;
  }
}

// -----------------------------------------------------------------------------
// Description: Test constant expression
// -----------------------------------------------------------------------------
void Test_ExpressionFunction_constant()
{
  ExpressionFunction f("1.0");

  CheckEqual(f(0.0), 1.0, 1.0e-14,
             "Constant expression evaluated incorrectly at x = 0.0");
  CheckEqual(f(0.5), 1.0, 1.0e-14,
             "Constant expression evaluated incorrectly at x = 0.5");
  CheckEqual(f(1.0), 1.0, 1.0e-14,
             "Constant expression evaluated incorrectly at x = 1.0");
}

// -----------------------------------------------------------------------------
// Description: Test identity expression
// -----------------------------------------------------------------------------
void Test_ExpressionFunction_identity()
{
  ExpressionFunction f("x");

  CheckEqual(f(0.0), 0.0, 1.0e-14,
             "Identity expression evaluated incorrectly at x = 0.0");
  CheckEqual(f(0.25), 0.25, 1.0e-14,
             "Identity expression evaluated incorrectly at x = 0.25");
  CheckEqual(f(0.75), 0.75, 1.0e-14,
             "Identity expression evaluated incorrectly at x = 0.75");
}

// -----------------------------------------------------------------------------
// Description: Test trigonometric expression with pi
// -----------------------------------------------------------------------------
void Test_ExpressionFunction_sine()
{
  ExpressionFunction f("sin(2*pi*x)");

  CheckEqual(f(0.0), 0.0, 1.0e-14,
             "Sine expression evaluated incorrectly at x = 0.0");
  CheckEqual(f(0.25), 1.0, 1.0e-12,
             "Sine expression evaluated incorrectly at x = 0.25");
  CheckEqual(f(0.5), 0.0, 1.0e-12,
             "Sine expression evaluated incorrectly at x = 0.5");
}

// -----------------------------------------------------------------------------
// Description: Test Gaussian-like expression
// -----------------------------------------------------------------------------
void Test_ExpressionFunction_gaussian()
{
  ExpressionFunction f("exp(-100*(x-0.5)^2)");

  CheckEqual(f(0.5), 1.0, 1.0e-14,
             "Gaussian expression evaluated incorrectly at x = 0.5");
  Check(f(0.0) < f(0.25),
        "Gaussian expression should increase towards the center");
  Check(f(0.25) < f(0.5),
        "Gaussian expression should increase towards the center");
}

// -----------------------------------------------------------------------------
// Description: Invalid syntax should throw
// -----------------------------------------------------------------------------
void Test_ExpressionFunction_invalidExpressionThrows()
{
  Check(ExpressionFunctionEvalThrows("sin(2*pi*)"),
        "Invalid expression syntax should throw");
}

// -----------------------------------------------------------------------------
// Description: Unknown symbol should throw
// -----------------------------------------------------------------------------
void Test_ExpressionFunction_unknownSymbolThrows()
{
  Check(ExpressionFunctionEvalThrows("sin(2*pi*y)"),
        "Unknown symbol in expression should throw");
}
