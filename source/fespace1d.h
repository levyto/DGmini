// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Finite element space class
//
// -----------------------------------------------------------------------------

#ifndef FESPACE1D_H
#define FESPACE1D_H

#include "Vec.h"
#include "Mat.h"
#include "quadrature1d.h"

// -----------------------------------------------------------------------------
// Description: Finite element space on one 1D element. This class stores all 
//              the info about the basis functions, quadrature and mass 
//              and stiffness matrices on reference element. 
// -----------------------------------------------------------------------------
class FESpace1D
{
  public:
    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------
    explicit FESpace1D(int p);

    // -------------------------------------------------------------------------
    // Access
    // -------------------------------------------------------------------------
    int    order()            const { return p_;                 }
    int    DoFs()             const { return p_ + 1;             }
    int    nip()              const { return quadrature_.nip();  }
    double phi(int q, int i)  const { return phi_q_(q, i);       }
    double dphi(int q, int i) const { return dphi_q_(q, i);      }
    double phiLeft(int i)     const { return phi_left_[i];       }
    double phiRight(int i)    const { return phi_right_[i];      }

    const Quadrature1D& quadrature() const { return quadrature_; }

    const Mat& massMatrix()          const { return M_;          }
    const Mat& inverseMassMatrix()   const { return M_inv_;      }
    const Mat& stiffnessMatrix()     const { return K_;          }

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
    int p_ = 0;

    Quadrature1D quadrature_;

    Vec phi_left_;
    Vec phi_right_;

    Mat phi_q_;
    Mat dphi_q_;
    Mat M_;
    Mat M_inv_;
    Mat K_;
};

#endif
