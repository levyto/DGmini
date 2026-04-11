// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Mass matrix on one element
//
// -----------------------------------------------------------------------------

#include <cassert>

#include "mass_matrix.h"
#include "basis1d.h"

// -----------------------------------------------------------------------------
// Description: Build local mass matrix on reference 1D element assuming the 
//              order of the basis functions is p.
//
//              Note that the mass matrix on the physical element can be 
//              obtained from the mass matrix on the reference element as: 
//                  M_e_physical = M_e_reference * J,
//              where J is the Jacobian of the element.
// -----------------------------------------------------------------------------
void buildMassMatrix1D(const Quadrature1D& quadrature, int p, Mat& M_e)
{
  assert(p >= 0);

  const int n_dofs = p + 1;
  
  assert(M_e.rows() == n_dofs);
  assert(M_e.cols() == n_dofs);

  Vec phi(n_dofs);
  Vec dphi(n_dofs);

  for (int q = 0; q < quadrature.nip(); ++q)
  {
    const double xi = quadrature.point(q);
    const double w  = quadrature.weight(q);

    evaluateLegendreBasis(p, xi, phi, dphi);

    for (int i = 0; i < n_dofs; ++i)
    {
      for (int j = 0; j < n_dofs; ++j)
      {
        M_e(i, j) += w * phi[i] * phi[j];
      }
    }
  }
}

Mat buildMassMatrix1D(const Quadrature1D& quadrature, int p)
{
  Mat M_e(p + 1, p + 1);
  buildMassMatrix1D(quadrature, p, M_e);
  return M_e;
}

// -----------------------------------------------------------------------------
// Description: Build local mass matrix on reference 1D element assuming the 
//              Legendre basis functions are used and the order of the basis 
//              functions is p for every element, i.e, the mass matrix is
//              diagonal and the diagonal entries can be computed as: 
//                  M_e(n,n) = 2/(2*n+1).
//                
//              Note that the mass matrix on the physical element can be 
//              obtained from the mass matrix on the reference element as: 
//                  M_e_physical = M_e_reference * J,
//              where J is the Jacobian of the element.
// -----------------------------------------------------------------------------
void buildMassMatrix1D(int p, Mat& M_e)
{
  assert(p >= 0);

  const int n_dofs = p + 1;

  assert(M_e.rows() == n_dofs);
  assert(M_e.cols() == n_dofs);

  for (int i = 0; i < n_dofs; i++)
  {
    M_e(i, i) = 2.0 / (2.0 * i + 1.0);
  }
}

Mat buildMassMatrix1D(int p)
{  
  Mat M_e(p + 1, p + 1);
  buildMassMatrix1D(p, M_e);
  return M_e;
}

// -----------------------------------------------------------------------------
// Description: Build local mass matrix inverse on reference 1D element assuming 
//              the Legendre basis functions are used and the order of the basis 
//              functions is p for every element, i.e, the mass matrix is 
//              diagonal and the diagonal entries can be computed as: 
//                  M_e(n,n) = 2/(2*n+1).
//                
//              Note that the mass matrix on the physical element can be 
//              obtained from the mass matrix on the reference element as: 
//                  M_e_physical = M_e_reference * J,
//              where J is the Jacobian of the element.
// -----------------------------------------------------------------------------
void buildMassMatrix1DInverse(int p, Mat& M_e_inv)
{
  assert(p >= 0);

  const int n_dofs = p + 1;

  assert(M_e_inv.rows() == n_dofs);
  assert(M_e_inv.cols() == n_dofs);

  for (int i = 0; i < n_dofs; i++)
  {
    M_e_inv(i, i) = (2.0 * i + 1.0) / 2.0;
  }
}

Mat buildMassMatrix1DInverse(int p)
{
  Mat M_e_inv(p + 1, p + 1);
  buildMassMatrix1DInverse(p, M_e_inv);
  
  return M_e_inv;
}