// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Stiffness matrix on one element
//
// -----------------------------------------------------------------------------

#include <cassert>

#include "FEM/stiffness_matrix.h"
#include "FEM/basis1d.h"

// -----------------------------------------------------------------------------
// Description: Build local stiffness matrix on one 1D element assuming the order 
//              of the basis functions is p for every element, i.e the Element1D 
//              does not store info about p, quadrature nor basis functions.
//
//              The stiffness matrix is computed as:
//                  K_e(i,j) = \int_Omega_e dphi_i(xi) * phi_j(xi)
// -----------------------------------------------------------------------------
void buildStiffnessMatrix1D(const Quadrature1D& quadrature, int p, Mat& K_e)
{
  assert(p >= 0);

  const int n_dofs = p + 1;
  
  assert(K_e.rows() == n_dofs);
  assert(K_e.cols() == n_dofs);

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
        K_e(i, j) += w * dphi[i] * phi[j];
      }
    }
  }
}

Mat buildStiffnessMatrix1D(const Quadrature1D& quadrature, int p)
{
  Mat K_e(p + 1, p + 1);
  buildStiffnessMatrix1D(quadrature, p, K_e);
  return K_e;
}

// -----------------------------------------------------------------------------
// Description: Build local stiffness matrix on one 1D element assuming the order 
//              of the basis functions is p for every element, i.e the Element1D 
//              does not store info about p, quadrature nor basis functions.
//
//              The stiffness matrix is computed as:
//                  K_e(i,j) = \int_{-1}^{1} dphi_i(xi) * phi_j(xi)
//
//              For Legendre basis functions, the stiffness matrix is lower 
//              triangular and nonzero entries are only on the positions i > j 
//              and (i-j) is odd (alternating parity of Legendre polynomials), 
//              where the value is 2.
//              
//              E.g, for p = 5, the stiffness matrix has this structure:
//                i\j   0   1   2   3   4   5
//                --------------------------------
//                0     0   0   0   0   0   0
//                1     2   0   0   0   0   0
//                2     0   2   0   0   0   0
//                3     2   0   2   0   0   0
//                4     0   2   0   2   0   0
//                5     2   0   2   0   2   0
// -----------------------------------------------------------------------------
void buildStiffnessMatrix1D(int p, Mat& K_e)
{
  assert(p >= 0);

  const int n_dofs = p + 1;

  assert(K_e.rows() == n_dofs);
  assert(K_e.cols() == n_dofs);

  for (int i = 0; i <= p; ++i)
  {
    for (int j = 0; j <= p; ++j)
    {
      if (i > j && ((i - j) % 2 == 1))
        K_e(i,j) = 2.0;
      else
        K_e(i,j) = 0.0;
    }
  }
}

Mat buildStiffnessMatrix1D(int p)
{
  Mat K_e(p + 1, p + 1);
  buildStiffnessMatrix1D(p, K_e);
  return K_e;
}
