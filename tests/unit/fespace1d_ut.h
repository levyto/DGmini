// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for 1D finite element space class
//
// -----------------------------------------------------------------------------

#ifndef FESPACE1D_UT_H
#define FESPACE1D_UT_H

#include "unittest.h"
#include "FEM/fespace1d.h"

// -----------------------------------------------------------------------------
// Description: FESpace1D UTs
// -----------------------------------------------------------------------------
void Test_FESpace1D_constructor();

// -----------------------------------------------------------------------------
// Description: FESpace1D UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_FESpace1D(TestRegistry& registry)
{
  registry.add("Test_FESpace1D_constructor", Test_FESpace1D_constructor);
}

#endif
