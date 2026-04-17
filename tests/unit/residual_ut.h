// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for residual assembly
//
// -----------------------------------------------------------------------------

#ifndef RESIDUAL_UT_H
#define RESIDUAL_UT_H

#include "unittest.h"
#include "Spatial/residual.h"

// -----------------------------------------------------------------------------
// Description: Residual UTs
// -----------------------------------------------------------------------------
void Test_residual_linearAdvectionZeroSolution();
void Test_residual_linearAdvectionConstant();
void Test_residual_burgersConstant();
void Test_residual_linearAdvectionP0MatchesFV();
void Test_residual_linearAdvectionScalingWithVelocity();

// -----------------------------------------------------------------------------
// Description: Residual UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_residual(TestRegistry& registry)
{
    registry.add("Test_residual_linearAdvectionZeroSolution",        Test_residual_linearAdvectionZeroSolution       );
    registry.add("Test_residual_linearAdvectionConstant",            Test_residual_linearAdvectionConstant           );
    registry.add("Test_residual_burgersConstant",                    Test_residual_burgersConstant                   );
    registry.add("Test_residual_linearAdvectionP0MatchesFV",         Test_residual_linearAdvectionP0MatchesFV        );
    registry.add("Test_residual_linearAdvectionScalingWithVelocity", Test_residual_linearAdvectionScalingWithVelocity);
}

#endif
