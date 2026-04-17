// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: ExpressionFunction class definition
//
// -----------------------------------------------------------------------------

#ifndef EXPRESSION_FUNCTION_H
#define EXPRESSION_FUNCTION_H

#include <cmath>
#include <stdexcept>
#include <string>
#include <muParser.h>

class ExpressionFunction
{
  public:
    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------
    explicit ExpressionFunction(const std::string& expr)
    {
      parser_.DefineVar("x", &x_);
      parser_.DefineVar("t", &t_);
      parser_.DefineConst("pi", M_PI);

      try
      {
        parser_.SetExpr(expr);
      }
      catch (mu::Parser::exception_type& e)
      {
        throw std::runtime_error(
          "Invalid expression '" + expr + "': " + e.GetMsg());
      }
    }

    // -------------------------------------------------------------------------
    // Access
    // -------------------------------------------------------------------------
    double operator()(double x) const
    {
      return (*this)(x, 0.0);
    }

    double operator()(double x, double t) const
    {
      x_ = x;
      t_ = t;

      try
      {
        return parser_.Eval();
      }
      catch (mu::Parser::exception_type& e)
      {
        throw std::runtime_error("Expression evaluation failed: " + e.GetMsg());
      }
    }

    // -------------------------------------------------------------------------
    // Modification
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Operations
    // -------------------------------------------------------------------------

  private:
    // -------------------------------------------------------------------------
    // Data
    // -------------------------------------------------------------------------
    mutable double x_ = 0.0;
    mutable double t_ = 0.0;
    mu::Parser parser_;
};

#endif
