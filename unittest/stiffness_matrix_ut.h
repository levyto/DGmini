// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for stiffness matrix on one element
//
// -----------------------------------------------------------------------------

#ifndef STIFFNESSMATRIX_UT_H
#define STIFFNESSMATRIX_UT_H

#include "unittest.h"
#include "FEM/stiffness_matrix.h"

// -----------------------------------------------------------------------------
// Description: Stiffness matrix UTs
// -----------------------------------------------------------------------------
void Test_stiffnessMatrix1D_numVsAnalytic();

// -----------------------------------------------------------------------------
// Description: Stiffness matrix UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_stiffnessMatrix1D(TestRegistry& registry)
{
  registry.add("Test_stiffnessMatrix1D_numVsAnalytic", Test_stiffnessMatrix1D_numVsAnalytic);
}

#endif
