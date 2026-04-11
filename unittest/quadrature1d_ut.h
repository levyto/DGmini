// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for Quadrature1D class
//
// -----------------------------------------------------------------------------

#ifndef QUADRATURE1D_UT_H
#define QUADRATURE1D_UT_H

#include "unittest.h"
#include "FEM/quadrature1d.h"

// -----------------------------------------------------------------------------
// Description: Quadrature1D UTs
// -----------------------------------------------------------------------------
void Test_Quadrature1D_weightsSum();
void Test_Quadrature1D_symmetry();
void Test_Quadrature1D_exactness();
void Test_Quadrature1D_elemIntegration();

// -----------------------------------------------------------------------------
// Description: Quadrature1D UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_Quadrature1D(TestRegistry& registry)
{
  registry.add("Test_Quadrature1D_weightsSum",       Test_Quadrature1D_weightsSum      );
  registry.add("Test_Quadrature1D_symmetry",         Test_Quadrature1D_symmetry        );
  registry.add("Test_Quadrature1D_exactness",        Test_Quadrature1D_exactness       );
  registry.add("Test_Quadrature1D_elemIntegration",  Test_Quadrature1D_elemIntegration );
}

#endif