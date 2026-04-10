// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittest runner for DGmini building blocks
//
// -----------------------------------------------------------------------------

#include <iostream>

#include "unittest.h"
#include "Vec_ut.h"
#include "linalg_ut.h"
#include "element1d_ut.h"
#include "mesh1d_ut.h"

int main()
{
  TestRegistry registry;

  Register_Test_Vec(registry);
  Register_Test_linalg(registry);
  Register_Test_Element1D(registry);
  Register_Test_Mesh1D(registry);

  int n_failed = 0;

  for (const auto& [name, test] : registry.tests())
  {
    try
    {
      test();
      std::cout << "[PASS] " << name << '\n';
    }
    catch (const std::exception& e)
    {
      n_failed++;
      std::cout << "[FAIL] " << name << " : " << e.what() << '\n';
    }
  }

  if (n_failed == 0)
  {
    std::cout << "\nAll tests passed.\n";
    return 0;
  }
  else
  {
    std::cout << "\n" << n_failed << " test(s) failed.\n";
    return 1;
  }
}
