// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Unittests for runtime factory
//
// -----------------------------------------------------------------------------

#ifndef RUNTIME_FACTORY_UT_H
#define RUNTIME_FACTORY_UT_H

#include "unittest.h"
#include "Solver/runtime_factory.h"

// -----------------------------------------------------------------------------
// Description: Runtime factory UTs
// -----------------------------------------------------------------------------
void Test_RuntimeFactory_createBoundaryCondition();
void Test_RuntimeFactory_createObjects();
void Test_RuntimeFactory_invalidInputsThrow();

// -----------------------------------------------------------------------------
// Description: Runtime factory UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_RuntimeFactory(TestRegistry& registry)
{
  registry.add("Test_RuntimeFactory_createBoundaryCondition",
                Test_RuntimeFactory_createBoundaryCondition);
  registry.add("Test_RuntimeFactory_createObjects",
                Test_RuntimeFactory_createObjects);
  registry.add("Test_RuntimeFactory_invalidInputsThrow",
                Test_RuntimeFactory_invalidInputsThrow);
}

#endif
