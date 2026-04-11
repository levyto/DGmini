// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for 1D basis functions and their derivatives
//
// -----------------------------------------------------------------------------

#include <algorithm>
#include <string>
#include <vector>

#include "basis1d_ut.h"
#include "quadrature1d.h"

// -----------------------------------------------------------------------------
// Description: Check  P(1) == 1,         P(-1) == (-1)^n,
//                    dP(1) == n(n+1)/2, dP(-1) == (-1)^(n+1)*n(n+1)/2
// -----------------------------------------------------------------------------
void Test_basis1D_boundaryValues()
{
  const int pmax = 8;
  Vec P(pmax + 1), dP(pmax + 1);

  evaluateLegendreBasis(pmax, 1.0, P, dP);

  for (int p = 0; p < pmax+1; p++)
  {
    CheckEqual( P[p], 1.0,               1.e-14, 
               "Boundary value not correct for P("  + std::to_string(p) + ") at xi=1.0");
    CheckEqual(dP[p], 0.5 * p * (p + 1), 1.e-14, 
               "Boundary value not correct for P'(" + std::to_string(p) + ") at xi=1.0");
  }
  

  evaluateLegendreBasis(pmax, -1.0, P, dP);

  for (int p = 0; p < pmax+1; p++)
  {
    CheckEqual( P[p], (p % 2 == 0) ? 1.0 : -1.0,                     1.e-14, 
               "Boundary value not correct for P("  + std::to_string(p) + ") at xi=-1.0");
    CheckEqual(dP[p], (((p+1) % 2 == 0) ? 1.0 : -1.0) * p*(p+1)/2.0, 1.e-14, 
               "Boundary value not correct for P'(" + std::to_string(p) + ") at xi=-1.0");
  }
}

// -----------------------------------------------------------------------------
// Description: Check   (P_n,P_m) = 0 for n != m 
//              and     (P_n,P_n) = 2/(2n+1) 
//              where the inner product is defined as
//                          (f,g) = int_{-1}^{1} f(x)*g(x) dx
// -----------------------------------------------------------------------------
void Test_basis1D_orthogonality()
{
  const int pmax = 6;
  
  for (int p1 = 0; p1 < pmax; p1++)
  {
    for (int p2 = 0; p2 < pmax; p2++)
    {
      double integral = 0.0;

      Quadrature1D quadrature( p1 + p2 ); 
     
      for (int q = 0; q < quadrature.nip(); ++q)
      {
        double xi = quadrature.point(q);
        double w  = quadrature.weight(q);

        Vec P1(p1 + 1), dP1(p1 + 1);
        Vec P2(p2 + 1), dP2(p2 + 1);

        evaluateLegendreBasis(p1, xi, P1, dP1);
        evaluateLegendreBasis(p2, xi, P2, dP2);

        integral += w * P1[p1] * P2[p2];
      }

      if (p1 == p2)
      {
        CheckEqual(integral, 2.0 / (2.0 * p1 + 1.0), 1.e-14, 
                   "Orthogonality check failed for P(" + std::to_string(p1) + ")" 
                                           + " and P(" + std::to_string(p2) + ")");
      }
      else
      {
        CheckEqual(integral, 0.0,                    1.e-14, 
                   "Orthogonality check failed for P(" + std::to_string(p1) + ")" 
                                           + " and P(" + std::to_string(p2) + ")");
      }      
    }    
  } 
}

// -----------------------------------------------------------------------------
// Description: Check exactness of P0, P1, P2, P3 and their derivatives 
//              at a few points in [-1,1]
//
// Note: P0 = 1                dP0 = 0
//       P1 = x                dP1 = 1
//       P2 = 0.5*(3x^2 - 1)   dP2 = 3x  
//       P3 = 0.5*(5x^3 - 3x)  dP3 = 0.5*(15x^2 - 3) 
// -----------------------------------------------------------------------------
void Test_basis1D_lowOrderExactness()
{
  std::vector<double> xs = {-1.0, -0.5, 0.0, 0.5, 1.0};

  for (double x : xs)
  {
    Vec P(4), dP(4);
    evaluateLegendreBasis(3, x, P, dP);

    CheckEqual(P[0], 1.0,                   1.e-14, "Exactness check failed for P(0)");
    CheckEqual(P[1], x,                     1.e-14, "Exactness check failed for P(1)");
    CheckEqual(P[2], 0.5 * (3*x*x - 1.0),   1.e-14, "Exactness check failed for P(2)");
    CheckEqual(P[3], 0.5 * (5*x*x*x - 3*x), 1.e-14, "Exactness check failed for P(3)");

    CheckEqual(dP[0], 0.0,                  1.e-14, "Exactness check failed for P'(0)");
    CheckEqual(dP[1], 1.0,                  1.e-14, "Exactness check failed for P'(1)");
    CheckEqual(dP[2], 3.0 * x,              1.e-14, "Exactness check failed for P'(2)");
    CheckEqual(dP[3], 0.5 * (15*x*x - 3.0), 1.e-14, "Exactness check failed for P'(3)");
  }
}