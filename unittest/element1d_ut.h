// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for Element1D class
//
// -----------------------------------------------------------------------------

#ifndef ELEMENT1D_UT_H
#define ELEMENT1D_UT_H

#include "unittest.h"
#include "Mesh/element1d.h"

// -----------------------------------------------------------------------------
// Description: Element1D UTs
// -----------------------------------------------------------------------------
void Test_Element1D_basic();

// -----------------------------------------------------------------------------
// Description: Element1D UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_Element1D(TestRegistry& registry)
{
  registry.add("Test_Element1D_basic", Test_Element1D_basic);
}

#endif
