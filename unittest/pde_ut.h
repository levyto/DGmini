// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for PDE class
//
// -----------------------------------------------------------------------------

#ifndef PDE_UT_H
#define PDE_UT_H

#include "unittest.h"
#include "PDE/pde.h"

// -----------------------------------------------------------------------------
// Description: PDE UTs
// -----------------------------------------------------------------------------
void Test_LinearAdvection1D();
void Test_Burgers1D ();

// -----------------------------------------------------------------------------
// Description: PDE UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_PDE(TestRegistry& registry)
{
  registry.add("Test_LinearAdvection1D", Test_LinearAdvection1D);
  registry.add("Test_Burgers1D",         Test_Burgers1D        );
}

#endif