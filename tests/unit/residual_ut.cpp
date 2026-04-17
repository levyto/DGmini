// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for residual assembly
//
// -----------------------------------------------------------------------------

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "residual_ut.h"
#include "IO/expression_function.h"
#include "Mesh/boundary_conditions1d.h"
#include "PDE/burgers1d.h"
#include "PDE/linear_advection1d.h"
#include "Spatial/NumericalFlux/rusanov.h"

// -----------------------------------------------------------------------------
// Description: Zero solution should yield zero residual for any PDE
// -----------------------------------------------------------------------------
void Test_residual_linearAdvectionZeroSolution()
{
  const int p = 2;
  const int Ne = 3;
  const double tol = 1e-14;

  const double t = 0.0;

  BoundaryConditions1D bc;
  bc.left.type  = BoundaryConditionType::Periodic;
  bc.right.type = BoundaryConditionType::Periodic;

  FESpace1D fe(p);
  Mesh1D mesh(0.0, 1.0, Ne);
  LinearAdvection1D pde(1.0);
  RusanovFlux flux;

  ModalVector u(Ne, fe.DoFs());
  ModalVector rhs(Ne, fe.DoFs());

  u.zero();

  residual(fe, mesh, pde, flux, bc, t, u, rhs);

  for (int e = 0; e < Ne; ++e)
  {
    for (int i = 0; i < fe.DoFs(); ++i)
    {
      CheckEqual(rhs(e, i), 0.0, tol, "Zero solution residual is not zero");
    }
  }
}

// -----------------------------------------------------------------------------
// Description: For linear advection, test for zero residual with constant solution
// -----------------------------------------------------------------------------
void Test_residual_linearAdvectionConstant()
{
  const int p = 3;
  const int Ne = 4;
  const double tol = 1e-12;

  const double t = 0.0;

  BoundaryConditions1D bc;
  bc.left.type  = BoundaryConditionType::Periodic;
  bc.right.type = BoundaryConditionType::Periodic;

  FESpace1D fe(p);
  Mesh1D mesh(0.0, 1.0, Ne);
  LinearAdvection1D pde(2.0);
  RusanovFlux flux;

  ModalVector u(Ne, fe.DoFs());
  ModalVector rhs(Ne, fe.DoFs());

  const double C = 3.0;

  for (int e = 0; e < Ne; ++e)
  {
    u(e, 0) = C;
    for (int i = 1; i < fe.DoFs(); ++i)
    {
      u(e, i) = 0.0;
    }
  }

  residual(fe, mesh, pde, flux, bc, t, u, rhs);

  for (int e = 0; e < Ne; ++e)
  {
    for (int i = 0; i < fe.DoFs(); ++i)
    {
      CheckEqual(rhs(e, i), 0.0, tol, "Constant linear advection residual is not zero");
    }
  }
}

// -----------------------------------------------------------------------------
// Description: For Burgers' equation, test for zero residual with constant solution
// -----------------------------------------------------------------------------
void Test_residual_burgersConstant()
{
  const int p = 3;
  const int Ne = 4;
  const double tol = 1e-12;

  const double t = 0.0;

  BoundaryConditions1D bc;
  bc.left.type  = BoundaryConditionType::Periodic;
  bc.right.type = BoundaryConditionType::Periodic;

  FESpace1D fe(p);
  Mesh1D mesh(0.0, 1.0, Ne);
  Burgers1D pde;
  RusanovFlux flux;

  ModalVector u(Ne, fe.DoFs());
  ModalVector rhs(Ne, fe.DoFs());

  const double C = 2.0;

  for (int e = 0; e < Ne; ++e)
  {
    u(e, 0) = C;
    for (int i = 1; i < fe.DoFs(); ++i)
    {
      u(e, i) = 0.0;
    }
  }

  residual(fe, mesh, pde, flux, bc, t, u, rhs);

  for (int e = 0; e < Ne; ++e)
  {
    for (int i = 0; i < fe.DoFs(); ++i)
    {
      CheckEqual(rhs(e, i), 0.0, tol, "Constant Burgers residual is not zero");
    }
  }
}

// -----------------------------------------------------------------------------
// Description: Test DG reduced to upwind advection with P0 basis functions, 
//              which has an analytical solution. 
// -----------------------------------------------------------------------------
void Test_residual_linearAdvectionP0MatchesFV()
{
  const int p = 0;
  const int Ne = 3;
  const double tol = 1e-12;

  const double t = 0.0;

  BoundaryConditions1D bc;
  bc.left.type  = BoundaryConditionType::Periodic;
  bc.right.type = BoundaryConditionType::Periodic;

  FESpace1D fe(p);
  Mesh1D mesh(0.0, 3.0, Ne); // deltax = 1
  LinearAdvection1D pde(1.0);
  RusanovFlux flux;

  ModalVector u(Ne, fe.DoFs());
  ModalVector rhs(Ne, fe.DoFs());

  u(0,0) = 1.0;
  u(1,0) = 2.0;
  u(2,0) = 4.0;

  residual(fe, mesh, pde, flux, bc, t, u, rhs);

  CheckEqual(rhs(0,0),  3.0, tol, "p=0 FV equivalence failed on element 0");
  CheckEqual(rhs(1,0), -1.0, tol, "p=0 FV equivalence failed on element 1");
  CheckEqual(rhs(2,0), -2.0, tol, "p=0 FV equivalence failed on element 2");
}

// -----------------------------------------------------------------------------
// Description: For linear advection, test for linear scaling of the residual 
//              with respect to advection speed
// -----------------------------------------------------------------------------
void Test_residual_linearAdvectionScalingWithVelocity()
{
  const int p = 1;
  const int Ne = 3;

  const double t = 0.0;

  BoundaryConditions1D bc;
  bc.left.type  = BoundaryConditionType::Periodic;
  bc.right.type = BoundaryConditionType::Periodic;

  FESpace1D fe(p);
  Mesh1D mesh(0.0, 1.0, Ne);
  RusanovFlux flux;

  ModalVector u(Ne, fe.DoFs()); 
  ModalVector rhs1(Ne, fe.DoFs()); 
  ModalVector rhs2(Ne, fe.DoFs());

  // set some non-trivial u
  for (int e = 0; e < Ne; ++e) 
    u(e, 0) = e + 1.0;

  residual(fe, mesh, LinearAdvection1D(1.0), flux, bc, t, u, rhs1);
  residual(fe, mesh, LinearAdvection1D(2.0), flux, bc, t, u, rhs2);

  for (int e = 0; e < Ne; ++e)
  {
    for (int i = 0; i < fe.DoFs(); ++i)
    {
      CheckEqual(rhs2(e, i), 2.0 * rhs1(e, i), 1e-12, 
                 "residual must scale linearly with advection speed");
    }
  }
}

// -----------------------------------------------------------------------------
// Description: For positive advection, the left Dirichlet BC should use the
//              prescribed boundary expression evaluated at the current time
// -----------------------------------------------------------------------------
void Test_residual_leftDirichletUsesBoundaryValueAndTime()
{
  const int p = 0;
  const int Ne = 1;
  const double tol = 1e-12;
  const double t = 0.25;

  BoundaryConditions1D bc;
  bc.left.type = BoundaryConditionType::Dirichlet;
  bc.left.expression = std::make_unique<ExpressionFunction>("t");
  bc.right.type = BoundaryConditionType::Outflow;

  FESpace1D fe(p);
  Mesh1D mesh(0.0, 1.0, Ne);
  LinearAdvection1D pde(1.0);
  RusanovFlux flux;

  ModalVector u(Ne, fe.DoFs());
  ModalVector rhs(Ne, fe.DoFs());

  u(0, 0) = 1.5;

  residual(fe, mesh, pde, flux, bc, t, u, rhs);

  CheckEqual(rhs(0, 0), -1.25, tol,
             "Left Dirichlet BC should affect residual through its time-dependent value");
}

// -----------------------------------------------------------------------------
// Description: For negative advection, the right Dirichlet BC should use the
//              prescribed boundary expression evaluated at boundary position
//              and current time
// -----------------------------------------------------------------------------
void Test_residual_rightDirichletUsesBoundaryValueAndTime()
{
  const int p = 0;
  const int Ne = 1;
  const double tol = 1e-12;
  const double t = 0.5;

  BoundaryConditions1D bc;
  bc.left.type = BoundaryConditionType::Outflow;
  bc.right.type = BoundaryConditionType::Dirichlet;
  bc.right.expression = std::make_unique<ExpressionFunction>("x + t");

  FESpace1D fe(p);
  Mesh1D mesh(0.0, 1.0, Ne);
  LinearAdvection1D pde(-2.0);
  RusanovFlux flux;

  ModalVector u(Ne, fe.DoFs());
  ModalVector rhs(Ne, fe.DoFs());

  u(0, 0) = 0.75;

  residual(fe, mesh, pde, flux, bc, t, u, rhs);

  CheckEqual(rhs(0, 0), 1.5, tol,
             "Right Dirichlet BC should affect residual through x- and t-dependent value");
}
