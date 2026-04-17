// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for Mat class
//
// -----------------------------------------------------------------------------

#ifndef MAT_UT_H
#define MAT_UT_H

#include "unittest.h"
#include "Algebra/Mat.h"

// -----------------------------------------------------------------------------
// Description: Mat UTs
// -----------------------------------------------------------------------------
void Test_Mat_basic();

// -----------------------------------------------------------------------------
// Description: Mat UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_Mat(TestRegistry& registry)
{
  registry.add("Test_Mat_basic", Test_Mat_basic);
}

#endif
