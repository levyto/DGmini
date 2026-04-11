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
void Test_linalg_gemv();
void Test_linalg_gemvHighLvl();
void Test_linalg_gemm();
void Test_linalg_gemmHighLvl();

// -----------------------------------------------------------------------------
// Description: linalg UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_linalg(TestRegistry& registry)
{       
  registry.add("Test_linalg_nrm2",        Test_linalg_nrm2);
  registry.add("Test_linalg_dot",         Test_linalg_dot );
  registry.add("Test_linalg_scal",        Test_linalg_scal);
  registry.add("Test_linalg_axpy",        Test_linalg_axpy);
  registry.add("Test_linalg_gemv",        Test_linalg_gemv);
  registry.add("Test_linalg_gemvHighLvl", Test_linalg_gemvHighLvl);
  registry.add("Test_linalg_gemm",        Test_linalg_gemm);
  registry.add("Test_linalg_gemmHighLvl", Test_linalg_gemmHighLvl);
}

#endif
