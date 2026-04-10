// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Gauss-Legendre quadrature rule on reference interval [-1,1]. 
//
//              This integration rule integrates polynomials of degree Q = 2n-1 
//              exactly, where n is the number of integration points. 
//
//              The points and weights are hard-coded till 20 integration 
//              points, which corresponds to polynomials of 39th degree 
//              integrated exactly. 
//
//              To be more generic, one would need to compute the roots of 
//              the Legendre polynomials and their derivatives.
//
// -----------------------------------------------------------------------------

#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "Vec.h"

// -----------------------------------------------------------------------------
// Description: Gauss-Legendre quadrature on [-1,1]
// -----------------------------------------------------------------------------
class Quadrature1D
{
  public:
    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------
    explicit Quadrature1D(int order);

    // -------------------------------------------------------------------------
    // Access
    // -------------------------------------------------------------------------
    int    nip()         const { return points_.size(); }
    double point(int i)  const { return points_[i];     }
    double weight(int i) const { return weights_[i];    }

    // -------------------------------------------------------------------------
    // Modification
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Operations
    // -------------------------------------------------------------------------

  private:
    // -------------------------------------------------------------------------
    // Data
    // -------------------------------------------------------------------------
    Vec points_;
    Vec weights_;
};

#endif