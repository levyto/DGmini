// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for ModalVector data structure
//
// -----------------------------------------------------------------------------

#ifndef MODAL_VECTOR_UT_H
#define MODAL_VECTOR_UT_H

#include "unittest.h"
#include "Spatial/modal_vector.h"

// -----------------------------------------------------------------------------
// Description: ModalVector UTs
// -----------------------------------------------------------------------------
void Test_ModalVector_construction();
void Test_ModalVector_indexing();
void Test_ModalVector_fillZero();
void Test_ModalVector_elementPtr();
void Test_ModalVector_axpy();

// -----------------------------------------------------------------------------
// Description: ModalVector UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_ModalVector(TestRegistry& registry)
{
  registry.add("Test_ModalVector_construction", Test_ModalVector_construction);
  registry.add("Test_ModalVector_indexing",     Test_ModalVector_indexing    );
  registry.add("Test_ModalVector_fillZero",     Test_ModalVector_fillZero    );
  registry.add("Test_ModalVector_elementPtr",   Test_ModalVector_elementPtr  );
  registry.add("Test_ModalVector_axpy",         Test_ModalVector_axpy        );
}

#endif
