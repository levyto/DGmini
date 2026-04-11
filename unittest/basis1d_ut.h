// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for 1D basis functions and their derivatives
//
// -----------------------------------------------------------------------------

#ifndef BASIS1D_UT_H
#define BASIS1D_UT_H

#include "unittest.h"
#include "basis1d.h"

// -----------------------------------------------------------------------------
// Description: Basis1D UTs
// -----------------------------------------------------------------------------
void Test_basis1D_boundaryValues();
void Test_basis1D_orthogonality();
void Test_basis1D_lowOrderExactness();

// -----------------------------------------------------------------------------
// Description: Basis1D UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_basis1D(TestRegistry& registry)
{
  registry.add("Test_basis1D_boundaryValues",    Test_basis1D_boundaryValues);
  registry.add("Test_basis1D_orthogonality",     Test_basis1D_orthogonality);
  registry.add("Test_basis1D_lowOrderExactness", Test_basis1D_lowOrderExactness);
}

#endif
