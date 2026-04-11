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
#include "quadrature1d_ut.h"
#include "basis1d_ut.h"
#include "Mat_ut.h"
#include "mass_matrix_ut.h"
#include "stiffness_matrix_ut.h"
#include "fespace1d_ut.h"
#include "pde_ut.h"
#include "numerical_flux_ut.h"

int main()
{
  auto& registry = GetTestRegistry();

  Register_Test_Vec(registry);
  Register_Test_linalg(registry);
  Register_Test_Element1D(registry);
  Register_Test_Mesh1D(registry);
  Register_Test_Quadrature1D(registry);
  Register_Test_basis1D(registry);
  Register_Test_Mat(registry);
  Register_Test_massMatrix1D(registry);
  Register_Test_stiffnessMatrix1D(registry);
  Register_Test_FESpace1D(registry);
  Register_Test_PDE(registry);
  Register_Test_NumericalFlux(registry);


  registry.run_all();


  std::cout << "\nChecks : " << registry.checks()   << "\n";

  if (registry.failures() == 0)
  {
    std::cout << "\n[OK] All tests passed\n\n";
    return 0;
  }
  else
  {
    std::cout << "\n[ERROR] " << registry.failures() << " check(s) failed\n\n";
    return 1;
  }
}
