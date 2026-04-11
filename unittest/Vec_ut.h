// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for Vec class
//
// -----------------------------------------------------------------------------

#ifndef VEC_UT_H
#define VEC_UT_H

#include "unittest.h"
#include "Algebra/Vec.h"

// -----------------------------------------------------------------------------
// Description: Vec UTs
// -----------------------------------------------------------------------------
void Test_Vec_basic();

// -----------------------------------------------------------------------------
// Description: Vec UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_Vec(TestRegistry& registry)
{
  registry.add("Test_Vec_basic", Test_Vec_basic);
}

#endif
