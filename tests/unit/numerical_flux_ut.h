// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for NumericalFlux class
//
// -----------------------------------------------------------------------------

#ifndef NUMERICALFLUX_UT_H
#define NUMERICALFLUX_UT_H

#include "unittest.h"
#include "PDE/pde.h"

// -----------------------------------------------------------------------------
// Description: NumericalFlux UTs
// -----------------------------------------------------------------------------
void Test_RusanovFlux_advectionPositive();
void Test_RusanovFlux_advectionNegative();
void Test_LaxFriedrichsFlux();

// -----------------------------------------------------------------------------
// Description: NumericalFlux UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_NumericalFlux(TestRegistry& registry)
{
  registry.add("Test_RusanovFlux_advectionPositive", Test_RusanovFlux_advectionPositive);
  registry.add("Test_RusanovFlux_advectionNegative", Test_RusanovFlux_advectionNegative);
  registry.add("Test_LaxFriedrichsFlux",             Test_LaxFriedrichsFlux            );
}

#endif