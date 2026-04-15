// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Compute effective CFL number for a given time step
//
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cassert>

#include "Temporal/cfl_number.h"

// -----------------------------------------------------------------------------
// Description: Compute effective CFL number for a given time step size dt, 
//              based on the maximum wave speed in the solution and the spatial 
//              discretization. The effective CFL number is computed as 
//
//              CFL = dt * (2p+1) * max_lambda / h
//
//              where p is the order of the DG method, max_lambda is the maximum 
//              wave speed, and h is the element size. This allows us to assess 
//              how close we are to the stability limit for a given time step.
// -----------------------------------------------------------------------------
double computeEffectiveCFL
(
  const FESpace1D& fe,
  const Mesh1D& mesh,
  const PDE& pde,
  const ModalVector& solution,
  double dt
)
{
  const double max_lambda = evalMaxConvectiveEigenvalue(fe, mesh, pde, solution);
  const double h = mesh.element(0).right() - mesh.element(0).left();

  return dt * ((2.0 * fe.order() + 1.0) * max_lambda) / h;
}

// -----------------------------------------------------------------------------
// Description: Evaluate the maximum convective eigenvalue (wave speed) in the
//              solution, which is needed to compute the effective CFL number. The 
//              maximum is taken over all elements and quadrature points, as well
//              as the left and right boundaries of each element, to account for
//              the fact that the wave speed may depend on the solution value.
// -----------------------------------------------------------------------------
double evalMaxConvectiveEigenvalue
(
  const FESpace1D& fe,
  const Mesh1D& mesh,
  const PDE& pde,
  const ModalVector& solution
)
{
  double max_lambda = 0.0;

  for (int e = 0; e < mesh.Ne(); ++e)
  {
    const double* u_e = solution.elementPtr(e);

    // Max over volume quadrature points
    for (int q = 0; q < fe.nip(); ++q)
    {
      double uq = 0.0;
      for (int j = 0; j < fe.DoFs(); ++j)
      {
        uq += u_e[j] * fe.phi(q, j);
      }
      max_lambda = std::max(max_lambda, pde.maxConvectiveEigenvalue(uq));
    }

    // Max over left and right boundaries
    double uL = 0.0;
    double uR = 0.0;
    for (int j = 0; j < fe.DoFs(); ++j)
    {
      uL += u_e[j] * fe.phiLeft(j);
      uR += u_e[j] * fe.phiRight(j);
    }

    max_lambda = std::max(max_lambda, pde.maxConvectiveEigenvalue(uL));
    max_lambda = std::max(max_lambda, pde.maxConvectiveEigenvalue(uR));
  }

  assert(max_lambda > 0.0);

  return max_lambda;
}