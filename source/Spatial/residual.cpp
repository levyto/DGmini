// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: du/dt = R residual assembly
//
// -----------------------------------------------------------------------------

#include <cassert>
#include <stdexcept>

#include "Spatial/residual.h"

// -----------------------------------------------------------------------------
// Description: Compute the residual R(u) for the given solution vector u, which 
//              represents the spatial discretization of the PDE
// -----------------------------------------------------------------------------
void residual
(
  const FESpace1D& fe,
  const Mesh1D& mesh,
  const PDE& pde,
  const NumericalFlux& flux,
  const BoundaryConditions1D& bc,
  const double time,
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

    /* Left state ----------------------------------------------------------- */
    double uL_minus = 0.0;

    if (e > 0) 
    {
      // INTERIOR
      const double* u_em1 = solution.elementPtr(e - 1);
      for (int j = 0; j < ndof; ++j)
      {
        uL_minus += u_em1[j] * fe.phiRight(j); // right trace of left neighbor
      }
    }
    else 
    {
      // LEFT BOUNDARY
      switch (bc.left.type)
      {
        case BoundaryConditionType::Periodic:
        {
          const double* u_em1 = solution.elementPtr(Ne - 1);
          for (int j = 0; j < ndof; ++j)
          {
            uL_minus += u_em1[j] * fe.phiRight(j); // right trace of left neighbor
          }
          break;
        }
        case BoundaryConditionType::Dirichlet:
        {
          uL_minus = (*(bc.left.expression))(element.left(), time);
          break;
        }
        case BoundaryConditionType::Outflow:
        {
          uL_minus = uL_plus;
          break;
        }
        default:
          throw std::runtime_error(
            "\n\nERROR: Unsupported left boundary condition type\n");
          break;
      }
    }

    /* Right state ---------------------------------------------------------- */
    double uR_plus = 0.0;

    if (e < Ne - 1) 
    {
      // INTERIOR
      const double* u_ep1 = solution.elementPtr(e + 1);
      for (int j = 0; j < ndof; ++j)
      {
        uR_plus += u_ep1[j] * fe.phiLeft(j); // left trace of right neighbor
      }
    }
    else 
    {
      // RIGHT BOUNDARY
      switch (bc.right.type)
      {
        case BoundaryConditionType::Periodic:
        {
          const double* u_ep1 = solution.elementPtr(0);

          for (int j = 0; j < ndof; ++j)
          {
            uR_plus += u_ep1[j] * fe.phiLeft(j); // left trace of right neighbor
          }
          break;
        }
        case BoundaryConditionType::Dirichlet:
        {
          uR_plus = (*(bc.right.expression))(element.right(), time);
          break;
        }
        case BoundaryConditionType::Outflow:
        {
          uR_plus = uR_minus;
          break;
        }
        default:
          throw std::runtime_error(
            "\n\nERROR: Unsupported right boundary condition type\n");
          break;
      }
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
