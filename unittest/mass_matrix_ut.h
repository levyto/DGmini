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
#include "FEM/mass_matrix.h"

// -----------------------------------------------------------------------------
// Description: Mass matrix UTs
// -----------------------------------------------------------------------------
void Test_massMatrix1D_numVsAnalyticLowLvl();
void Test_massMatrix1D_numVsAnalyticHighLvl();
void Test_massMatrix1D_inverseLowLvl();
void Test_massMatrix1D_inverseHighLvl();

// -----------------------------------------------------------------------------
// Description: Mass matrix UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_massMatrix1D(TestRegistry& registry)
{
  registry.add("Test_massMatrix1D_numVsAnalyticLowLvl",  Test_massMatrix1D_numVsAnalyticLowLvl );
  registry.add("Test_massMatrix1D_numVsAnalyticHighLvl", Test_massMatrix1D_numVsAnalyticHighLvl);
  registry.add("Test_massMatrix1D_inverseLowLvl",        Test_massMatrix1D_inverseLowLvl       );
  registry.add("Test_massMatrix1D_inverseHighLvl",       Test_massMatrix1D_inverseHighLvl      );
}

#endif
