// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for TimeIntegrator class
//
// -----------------------------------------------------------------------------

#ifndef TIME_INTEGRATOR_UT_H
#define TIME_INTEGRATOR_UT_H

#include "unittest.h"
#include "Spatial/modal_vector.h"
#include "Temporal/TimeIntegrator/time_integrator.h"
#include "PDE/pde.h"
#include "PDE/linear_advection1d.h"
#include "Spatial/NumericalFlux/numerical_flux.h"
#include "Spatial/NumericalFlux/rusanov.h"
#include "Spatial/l2_projection.h"

// -----------------------------------------------------------------------------
// Description: TimeIntegrator UTs
// -----------------------------------------------------------------------------
void Test_TimeIntegrator_ForwardEulerPreservesConstantSolution();
void Test_TimeIntegrator_RungeKutta2PreservesConstantSolution();
void Test_TimeIntegrator_RungeKutta3SSPPreservesConstantSolution();
void Test_TimeIntegrator_RungeKutta4PreservesConstantSolution();
void Test_TimeIntegrator_oneStepAccuracy();

// -----------------------------------------------------------------------------
// Description: TimeIntegrator UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_TimeIntegrator(TestRegistry& registry)
{
  registry.add("Test_TimeIntegrator_ForwardEulerPreservesConstantSolution",
                Test_TimeIntegrator_ForwardEulerPreservesConstantSolution);
  registry.add("Test_TimeIntegrator_RungeKutta2PreservesConstantSolution",
                Test_TimeIntegrator_RungeKutta2PreservesConstantSolution);
  registry.add("Test_TimeIntegrator_RungeKutta3SSPPreservesConstantSolution",
                Test_TimeIntegrator_RungeKutta3SSPPreservesConstantSolution);
  registry.add("Test_TimeIntegrator_RungeKutta4PreservesConstantSolution",
                Test_TimeIntegrator_RungeKutta4PreservesConstantSolution);
  registry.add("Test_TimeIntegrator_oneStepAccuracy", 
                Test_TimeIntegrator_oneStepAccuracy);
}

// ---------------------------------------------------------------------------
// Description: Forward declarations of helper functions for TimeIntegrator UTs
// -----------------------------------------------------------------------------
double modalErrorLinf(const ModalVector& a, const ModalVector& b);

// ---------------------------------------------------------------------------
// Description: Compute the one-step error of a time integrator for the linear
//              advection equation with a smooth initial condition, and check
//              that the error is consistent with the expected order of accuracy
// ---------------------------------------------------------------------------
template <typename Integrator>
double OneStepError(double dt)
{
  const int p = 5;
  const int Ne = 100;

  FESpace1D fe(p);
  Mesh1D mesh(0.0, 1.0, Ne);
  LinearAdvection1D pde(1.0);
  RusanovFlux flux;

  ModalVector u(mesh.Ne(), fe.DoFs());
  ModalVector u_exact(mesh.Ne(), fe.DoFs());

  auto u0 = [](double x)
  {
    return std::sin(2.0 * M_PI * x);
  };

  auto u_exact_fun = [dt](double x)
  {
    return std::sin(2.0 * M_PI * (x - dt));
  };

  for (int e = 0; e < mesh.Ne(); ++e)
  {
    L2ProjectionOnElement(fe, mesh.element(e), u0, u.elementPtr(e));
    L2ProjectionOnElement(fe, mesh.element(e), u_exact_fun, u_exact.elementPtr(e));
  }

  Integrator integrator;
  integrator.initialize(mesh, fe);
  integrator.doTimeStep(fe, mesh, pde, flux, dt, u);

  return modalErrorLinf(u, u_exact);
}

// -----------------------------------------------------------------------------
// Description: Test that the time integrator preserves a constant solution
// -----------------------------------------------------------------------------
template <typename Integrator>
void Test_TimeIntegrator_preservesConstantSolution()
{
  const int p = 2;
  const int Ne = 8;
  const double dt = 1.0e-3;
  const double tol = 1e-13;

  FESpace1D fe(p);
  Mesh1D mesh(0.0, 1.0, Ne);
  LinearAdvection1D pde(1.0);
  RusanovFlux flux;

  ModalVector u(mesh.Ne(), fe.DoFs());

  const double C = 2.5;
  for (int e = 0; e < mesh.Ne(); ++e)
  {
    u(e,0) = C;
    for (int i = 1; i < fe.DoFs(); ++i)
      u(e,i) = 0.0;
  }

  Integrator integrator;
  integrator.initialize(mesh, fe);
  integrator.doTimeStep(fe, mesh, pde, flux, dt, u);

  for (int e = 0; e < mesh.Ne(); ++e)
  {
    CheckEqual(u(e,0), C, tol, "Constant mode changed");
    for (int i = 1; i < fe.DoFs(); ++i)
      CheckEqual(u(e,i), 0.0, tol, "Higher mode changed for constant solution");
  }
}

#endif

