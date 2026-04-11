// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Abstract PDE class, intended to be used as a base class for 
//              specific PDEs (e.g. advection, diffusion, Burgers, etc.) 
//              and to be extended
//
// -----------------------------------------------------------------------------

#ifndef PDE_H
#define PDE_H

class PDE
{
  public:
    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------
    virtual ~PDE() = default;

    // -------------------------------------------------------------------------
    // Access
    // -------------------------------------------------------------------------
    virtual double convectiveFlux(double state) const { return 0.0; }
    virtual double diffusiveFlux(double state) const { return 0.0; }
    virtual double sourceTerm(double state, double x, double t) const { return 0.0; }

    // -------------------------------------------------------------------------
    // Modification
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Operations
    // -------------------------------------------------------------------------
    virtual double maxConvectiveEigenvalue(double state) const = 0;

  private:
    // -------------------------------------------------------------------------
    // Data
    // -------------------------------------------------------------------------
};

#endif