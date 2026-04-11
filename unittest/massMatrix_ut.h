// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for mass matrix on one element
//
// -----------------------------------------------------------------------------

#ifndef MASSMATRIX_UT_H
#define MASSMATRIX_UT_H

#include "unittest.h"
#include "massMatrix.h"

// -----------------------------------------------------------------------------
// Description: Mass matrix UTs
// -----------------------------------------------------------------------------
void Test_massMatrix1D_numVsAnalytic();
void Test_massMatrix1D_inverse();

// -----------------------------------------------------------------------------
// Description: Mass matrix UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_massMatrix1D(TestRegistry& registry)
{
  registry.add("Test_massMatrix1D_numVsAnalytic", Test_massMatrix1D_numVsAnalytic);
  registry.add("Test_massMatrix1D_inverse",       Test_massMatrix1D_inverse);
}

#endif
