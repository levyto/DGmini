// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for L2 projection
//
// -----------------------------------------------------------------------------

#ifndef L2_PROJECTION_UT_H
#define L2_PROJECTION_UT_H

#include "unittest.h"
#include "Spatial/l2_projection.h"

// -----------------------------------------------------------------------------
// Description: L2 projection UTs
// -----------------------------------------------------------------------------
void Test_L2ProjectionOnElement_constant();
void Test_L2ProjectionOnElement_linear();
void Test_L2ProjectionOnElement_legendreMode();

// -----------------------------------------------------------------------------
// Description: L2 projection UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_L2Projection(TestRegistry& registry)
{
  registry.add("Test_L2ProjectionOnElement_constant",     Test_L2ProjectionOnElement_constant    );
  registry.add("Test_L2ProjectionOnElement_linear",       Test_L2ProjectionOnElement_linear      );
  registry.add("Test_L2ProjectionOnElement_legendreMode", Test_L2ProjectionOnElement_legendreMode);
}

#endif
