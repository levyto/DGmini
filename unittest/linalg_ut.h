// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for linear algebra operations
//
// -----------------------------------------------------------------------------

#ifndef LINALG_UT_H
#define LINALG_UT_H

#include "unittest.h"
#include "linalg.h"

// -----------------------------------------------------------------------------
// Description: linalg UTs
// -----------------------------------------------------------------------------
void Test_linalg_nrm2();
void Test_linalg_dot();
void Test_linalg_scal();
void Test_linalg_axpy();

// -----------------------------------------------------------------------------
// Description: linalg UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_linalg(TestRegistry& registry)
{
  registry.add("Test_linalg_nrm2", Test_linalg_nrm2);
  registry.add("Test_linalg_dot",  Test_linalg_dot );
  registry.add("Test_linalg_scal", Test_linalg_scal);
  registry.add("Test_linalg_axpy", Test_linalg_axpy);
}

#endif
