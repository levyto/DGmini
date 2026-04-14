// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: L2 projection of a prescribed function onto the space 
//              of piecewise polynomials
//
// -----------------------------------------------------------------------------

#include <cassert>

#include "Spatial/residual.h"

// -----------------------------------------------------------------------------
// Description: L2 projection of a prescribed function u0 onto the space of 
//              polynomials of degree p on a given element. The projection is 
//              computed by solving the linear system M u = b, where M is the 
//              mass matrix on the reference element and b is the vector of 
//              inner products of u0 with the basis functions.
//
//              Note that the mass matrix is diagonal for Legendre basis 
//              functions, which simplifies the solution of the linear system.
// -----------------------------------------------------------------------------
void residual
(
  const FESpace1D& fe,
  const Mesh1D& mesh,
  const PDE& pde,
  const NumericalFlux& flux,
  const ModalVector& solution,
  ModalVector& rhs
)
{
  assert(solution.Ne() == mesh.Ne());
  assert(rhs.Ne() == mesh.Ne());
  assert(solution.localDoFs() == fe.DoFs());
  assert(rhs.localDoFs() == fe.DoFs());

  const int Ne = mesh.Ne();
  const int ndof = fe.DoFs();
  const int nip = fe.nip();

  rhs.zero();

  for (int e = 0; e < Ne; ++e)
  {
    const Element1D& element = mesh.element(e);
    const double J = element.jacobian();

    const double* u_e = solution.elementPtr(e);
    double* rhs_e = rhs.elementPtr(e);

    // -------------------------------------------------------------------------
    // Volume term                                  int_element f(u) * dphi_i dx
    // -------------------------------------------------------------------------
    for (int q = 0; q < nip; ++q)
    {
      double uq = 0.0;

      // Evaluate solution at quadrature point q
      for (int j = 0; j < ndof; ++j)
      {
        uq += u_e[j] * fe.phi(q, j);
      }

      // Evaluate flux at quadrature point q
      const double fq = pde.convectiveFlux(uq);      
      const double wq = fe.quadrature().weight(q);

      for (int i = 0; i < ndof; ++i)
      {
        rhs_e[i] += wq * fq * fe.dphi(q, i);
      }
    }

    // -------------------------------------------------------------------------
    // Boundary states
    //                            uL- | uL+          uR- | uR+
    //            |-------------------|------------------|------------------|
    //                     e-1        L        e         R       e+1
    // -------------------------------------------------------------------------
    double uL_plus  = 0.0;
    double uR_minus = 0.0;

    // Evaluate left and right traces of the solution at the element boundaries
    for (int j = 0; j < ndof; ++j)
    {
      uL_plus  += u_e[j] * fe.phiLeft(j);
      uR_minus += u_e[j] * fe.phiRight(j);
    }

    double uL_minus = 0.0;
    double uR_plus = 0.0;

    // periodic BC for now
    const int eL = (e == 0) ? (Ne - 1) : (e - 1);
    const int eR = (e == Ne - 1) ? 0 : (e + 1);

    const double* u_em1 = solution.elementPtr(eL);
    const double* u_ep1 = solution.elementPtr(eR);

    for (int j = 0; j < ndof; ++j)
    {
      uL_minus += u_em1[j] * fe.phiRight(j); // right trace of left neighbor
      uR_plus  += u_ep1[j] * fe.phiLeft(j);  // left trace of right neighbor
    }

    // -------------------------------------------------------------------------
    // Numerical fluxes
    // -------------------------------------------------------------------------
    const double fhatL = flux.evaluate(pde, uL_minus, uL_plus);
    const double fhatR = flux.evaluate(pde, uR_minus, uR_plus);

    // -------------------------------------------------------------------------
    // Surface term                      - phi_i(xR) * fhatR + phi_i(xL) * fhatL
    // -------------------------------------------------------------------------
    for (int i = 0; i < ndof; ++i)
    {
      rhs_e[i] += - fhatR * fe.phiRight(i)
                  + fhatL * fe.phiLeft(i);
    }

    // -------------------------------------------------------------------------
    // Apply inverse mass matrix                     rhs <- (1/J) M_ref^{-1} rhs
    // (assuming mass matrix is diagonal)
    // -------------------------------------------------------------------------
    for (int i = 0; i < ndof; ++i)
    {
      rhs_e[i] = (1.0 / J) * fe.inverseMassMatrix()(i, i) * rhs_e[i];
    }
  }
}
