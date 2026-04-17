// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for TimeIntegrator class
//
// -----------------------------------------------------------------------------

#include <cmath>

#include "time_integrator_ut.h"
#include "Temporal/TimeIntegrator/forward_euler.h"
#include "Temporal/TimeIntegrator/runge_kutta_2.h"
#include "Temporal/TimeIntegrator/runge_kutta_3_ssp.h"
#include "Temporal/TimeIntegrator/runge_kutta_4.h"

// -----------------------------------------------------------------------------
// Description: Helper functions for TimeIntegrator UTs
// -----------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // Description: Compute the L-infinity error between two ModalVectors, which 
  //              is the maximum absolute difference of their modal coefficients 
  //              across all elements and modes
  // ---------------------------------------------------------------------------
  double modalErrorLinf(const ModalVector& a, const ModalVector& b)
  {
    assert(a.Ne() == b.Ne());
    assert(a.DoFs() == b.DoFs());

    double err = 0.0;
    for (int e = 0; e < a.Ne(); ++e)
      for (int i = 0; i < a.localDoFs(); ++i)
        err = std::max(err, std::abs(a(e,i) - b(e,i)));

    return err;
  }

// -----------------------------------------------------------------------------
// Description: Test that the time integrators preserve a constant solution
// -----------------------------------------------------------------------------
void Test_TimeIntegrator_ForwardEulerPreservesConstantSolution()
{
  Test_TimeIntegrator_preservesConstantSolution<ForwardEuler>();
}

void Test_TimeIntegrator_RungeKutta2PreservesConstantSolution()
{
  Test_TimeIntegrator_preservesConstantSolution<RungeKutta2>();
}

void Test_TimeIntegrator_RungeKutta3SSPPreservesConstantSolution()
{
  Test_TimeIntegrator_preservesConstantSolution<RungeKutta3SSP>();
}

void Test_TimeIntegrator_RungeKutta4PreservesConstantSolution()
{
  Test_TimeIntegrator_preservesConstantSolution<RungeKutta4>();
}

// -----------------------------------------------------------------------------
// Description: Sanity check that the one-step error of the time integrators is
//              small for a smooth problem. The full order of accuracy is tested 
//              in the regression tests.
// -----------------------------------------------------------------------------
void Test_TimeIntegrator_oneStepAccuracy()
{
  double err = OneStepError<ForwardEuler>(1.0e-4);
  Check(err < 1.0e-3, "Forward Euler one-step error too large");

  err = OneStepError<RungeKutta2>(1.0e-4);
  Check(err < 1.0e-3, "Runge-Kutta 2 one-step error too large");

  err = OneStepError<RungeKutta3SSP>(1.0e-4);
  Check(err < 1.0e-3, "Runge-Kutta 3 SSP one-step error too large");

  err = OneStepError<RungeKutta4>(1.0e-4);
  Check(err < 1.0e-3, "Runge-Kutta 4 one-step error too large");
}
