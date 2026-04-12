// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver
// -----------------------------------------------------------------------------
//
// Description: Modal solution data structure
// -----------------------------------------------------------------------------

#ifndef MODAL_VECTOR_1D_H
#define MODAL_VECTOR_1D_H

#include <cassert>

#include "Algebra/Vec.h"

class ModalVector
{
  public:
    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------
    ModalVector() = default;

    ModalVector(int n_mesh_elements, int n_dofs)
      : n_mesh_elements_(n_mesh_elements), 
        n_dofs_(n_dofs), 
        data_(n_mesh_elements * n_dofs)
    {
      assert(n_mesh_elements > 0);
      assert(n_dofs > 0);
    }

    // -------------------------------------------------------------------------
    // Access
    // -------------------------------------------------------------------------
    int Ne()          const { return n_mesh_elements_;           }
    int DoFs()        const { return n_mesh_elements_ * n_dofs_; }
    int localDoFs()   const { return n_dofs_;                    }
    const Vec& data() const { return data_;                      }

    double operator()(int e, int i) const
    {
      assert(e >= 0 && e < n_mesh_elements_);
      assert(i >= 0 && i < n_dofs_);
      return data_[e * n_dofs_ + i];
    }

    int elementOffset(int e) const
    {
      assert(e >= 0 && e < n_mesh_elements_);
      return e * n_dofs_;
    }

    const double* elementPtr(int e) const
    {
      assert(e >= 0 && e < n_mesh_elements_);
      return &data_[e * n_dofs_];
    }

    // -------------------------------------------------------------------------
    // Modification
    // -------------------------------------------------------------------------
    Vec& data() { return data_; }

    double& operator()(int e, int i)
    {
      assert(e >= 0 && e < n_mesh_elements_);
      assert(i >= 0 && i < n_dofs_);
      return data_[e * n_dofs_ + i];
    }

    double* elementPtr(int e)
    {
      assert(e >= 0 && e < n_mesh_elements_);
      return &data_[e * n_dofs_];
    }

    void fill(double value) { data_.fill(value);}
    void zero()             { data_.fill(0.0);  }

    void axpy(double alpha, const ModalVector& x)
    {
      assert(n_mesh_elements_ == x.n_mesh_elements_);
      assert(n_dofs_ == x.n_dofs_);

      for (int k = 0; k < data_.size(); ++k)
      {
        data_[k] += alpha * x.data_[k];
      }
    }

    // -------------------------------------------------------------------------
    // Operations
    // -------------------------------------------------------------------------

  private:
    // -------------------------------------------------------------------------
    // Data
    // -------------------------------------------------------------------------
    int n_mesh_elements_ = 0;
    int n_dofs_ = 0;
    Vec data_;
};

#endif
