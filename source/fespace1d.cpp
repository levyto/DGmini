// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Finite element space class
//
// -----------------------------------------------------------------------------

#include <cassert>

#include "fespace1d.h"
#include "basis1d.h"
#include "mass_matrix.h"
#include "stiffness_matrix.h"

// -----------------------------------------------------------------------------
// Description: FESpace1D constructor
// -----------------------------------------------------------------------------
FESpace1D::FESpace1D(int p)
  : p_(p),
    quadrature_(2 * p),
    phi_q_(quadrature_.nip(), p + 1),
    dphi_q_(quadrature_.nip(), p + 1),
    phi_left_(p + 1),
    phi_right_(p + 1),
    M_(buildMassMatrix1D(quadrature_, p)),
    M_inv_(buildMassMatrix1DInverse(p)),
    K_(buildStiffnessMatrix1D(p))
{
  assert(p >= 0);
  
  Vec  phi(p + 1);
  Vec dphi(p + 1);

  for (int q = 0; q < quadrature_.nip(); ++q)
  {
    const double xi = quadrature_.point(q);
    evaluateLegendreBasis(p_, xi, phi, dphi);

    for (int i = 0; i <= p_; ++i)
    {
       phi_q_(q, i) = phi[i];
      dphi_q_(q, i) = dphi[i];
    }
  }

  evaluateLegendreBasis(p_, -1.0, phi_left_,  dphi);
  evaluateLegendreBasis(p_,  1.0, phi_right_, dphi);
}