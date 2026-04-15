// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for CFL number computation 
//
// -----------------------------------------------------------------------------

#ifndef CFL_NUMBER_UT_H
#define CFL_NUMBER_UT_H

#include "unittest.h"
#include "Temporal/cfl_number.h"

// -----------------------------------------------------------------------------
// Description: CFL number UTs
// -----------------------------------------------------------------------------
void Test_evalMaxConvectiveEigenvalue_CFLTimeStepLinearAdvectionConstantSolution();
void Test_evalMaxConvectiveEigenvalue_ReflectsElementBoundaries();
void Test_evalMaxConvectiveEigenvalue_BurgersP0UsesMaxSolutionValue();
void Test_computeEffectiveCFL_matchesFormula();

// -----------------------------------------------------------------------------
// Description: CFL number UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_CFLNumber(TestRegistry& registry)
{
  registry.add("Test_evalMaxConvectiveEigenvalue_CFLTimeStepLinearAdvectionConstantSolution",
                Test_evalMaxConvectiveEigenvalue_CFLTimeStepLinearAdvectionConstantSolution);
  registry.add("Test_evalMaxConvectiveEigenvalue_ReflectsElementBoundaries",
                Test_evalMaxConvectiveEigenvalue_ReflectsElementBoundaries);
  registry.add("Test_evalMaxConvectiveEigenvalue_BurgersP0UsesMaxSolutionValue",
                Test_evalMaxConvectiveEigenvalue_BurgersP0UsesMaxSolutionValue);
  registry.add("Test_computeEffectiveCFL_matchesFormula",
                Test_computeEffectiveCFL_matchesFormula);
}
 
#endif