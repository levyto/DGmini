// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: L2 projection of a prescribed function onto the space 
//              of piecewise polynomials
//
// -----------------------------------------------------------------------------

#include <cassert>

#include "Spatial/l2_projection.h"
#include "FEM/basis1d.h"

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
void L2ProjectionOnElement(const FESpace1D& fe, 
                           const Element1D& el,
                           std::function<double(double)> u0,
                           Vec& u_e)
{
  const int ndof = fe.DoFs();

  assert(u_e.size() == ndof);

  Vec phi(ndof);
  Vec dphi(ndof);

  Vec b(ndof);
  b.fill(0.0);

  const double J = el.jacobian();

  for (int q = 0; q < fe.quadrature().nip(); ++q)
  {
    const double xi = fe.quadrature().point(q);
    const double w  = fe.quadrature().weight(q);
    const double x  = el.mapToPhysical(xi);

    const double u0_val = u0(x);

    evaluateLegendreBasis(fe.order(), xi, phi, dphi);

    for (int i = 0; i < ndof; ++i)
    {
      b[i] += w * u0_val * phi[i] * J;
    }
  }

  // solve M u = b (assuming M is diagonal!)
  for (int i = 0; i < ndof; ++i)
  {
    u_e[i] = b[i] / ( fe.massMatrix()(i, i) * J );
  }
}