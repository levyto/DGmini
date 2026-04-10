// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for Mesh1D class
//
// -----------------------------------------------------------------------------

#ifndef MESH1D_UT_H
#define MESH1D_UT_H

#include "unittest.h"
#include "mesh1d.h"

// -----------------------------------------------------------------------------
// Description: Mesh1D UTs
// -----------------------------------------------------------------------------
void Test_Mesh1D_basic();

// -----------------------------------------------------------------------------
// Description: Mesh1D UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_Mesh1D(TestRegistry& registry)
{
  registry.add("Test_Mesh1D_basic", Test_Mesh1D_basic);
}

#endif
