// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Definition of BoundaryConditions1D struct, which holds the 
//              boundary condition information for the 1D solver
//
// -----------------------------------------------------------------------------

#ifndef BOUNDARY_CONDITIONS_1D_H
#define BOUNDARY_CONDITIONS_1D_H

#include <memory>
#include <string>

#include "IO/expression_function.h"

enum class BoundaryConditionType
{
  Periodic,
  Dirichlet,
  Outflow
};

struct BoundaryCondition
{
  BoundaryConditionType type;
  std::unique_ptr<ExpressionFunction> expression;
};

struct BoundaryConditions1D
{
  BoundaryCondition left;
  BoundaryCondition right;
};

#endif
