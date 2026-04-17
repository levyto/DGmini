// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for ModalVector data structure
//
// -----------------------------------------------------------------------------

#include "modal_vector_ut.h"
#include "FEM/quadrature1d.h"

// -----------------------------------------------------------------------------
// Description: Check correct construction of ModalVector
// -----------------------------------------------------------------------------
void Test_ModalVector_construction()
{
  ModalVector u(4, 3);

  Check(u.Ne() == 4, "Wrong number of elements");
  Check(u.localDoFs() == 3, "Wrong local DoFs");
  Check(u.DoFs() == 12, "Wrong total DoFs");
}

// -----------------------------------------------------------------------------
// Description: Check correct indexing of ModalVector
// -----------------------------------------------------------------------------
void Test_ModalVector_indexing()
{
  ModalVector u(3, 2);

  u(0,0) = 1.0;
  u(0,1) = 2.0;
  u(1,0) = 3.0;
  u(1,1) = 4.0;
  u(2,0) = 5.0;
  u(2,1) = 6.0;

  CheckEqual(u(0,0), 1.0, 1e-14, "Indexing failed");
  CheckEqual(u(1,1), 4.0, 1e-14, "Indexing failed");
  CheckEqual(u(2,1), 6.0, 1e-14, "Indexing failed");
}

// -----------------------------------------------------------------------------
// Description: Check fill() and zero() methods of ModalVector
// -----------------------------------------------------------------------------
void Test_ModalVector_fillZero()
{
  ModalVector u(3, 4);

  u.fill(2.5);

  for (int e = 0; e < u.Ne(); ++e)
  {
    for (int i = 0; i < u.localDoFs(); ++i)
    {
      CheckEqual(u(e,i), 2.5, 1e-14, "Fill failed");
    }
  }

  u.zero();

  for (int e = 0; e < u.Ne(); ++e)
  {
    for (int i = 0; i < u.localDoFs(); ++i)
    {
      CheckEqual(u(e,i), 0.0, 1e-14, "Zero failed");
    }
  }
}

// -----------------------------------------------------------------------------
// Description: Check elementPtr() method of ModalVector
// -----------------------------------------------------------------------------
void Test_ModalVector_elementPtr()
{
  ModalVector u(3, 3);

  int counter = 0;
  for (int e = 0; e < u.Ne(); ++e)
  {
    for (int i = 0; i < u.localDoFs(); ++i)
    {
      u(e,i) = static_cast<double>(counter++);
    }
  }

  for (int e = 0; e < u.Ne(); ++e)
  {
    const double* ue = u.elementPtr(e);
    int offset = u.elementOffset(e);

    for (int i = 0; i < u.localDoFs(); ++i)
    {
      CheckEqual(ue[i], u(e,i), 1e-14, "elementPtr mismatch");
      CheckEqual(ue[i], u.data()[offset + i], 1e-14, "offset mismatch");
    }
  }
}

// -----------------------------------------------------------------------------
// Description: Check axpy operation of ModalVector
// -----------------------------------------------------------------------------
void Test_ModalVector_axpy()
{
  ModalVector u(2, 2);
  ModalVector x(2, 2);

  u.fill(1.0);
  x.fill(2.0);

  u.axpy(0.5, x); // u = 1 + 0.5*2 = 2

  for (int e = 0; e < u.Ne(); ++e)
  {
    for (int i = 0; i < u.localDoFs(); ++i)
    {
      CheckEqual(u(e,i), 2.0, 1e-14, "axpy failed");
    }
  }
}