// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Dense matrix wrapper for small local linear algebra objects.
//              Indexation is row-major, i.e., 
//              i * n_cols_ + j corresponds to (i,j) entry of the matrix.
//
// -----------------------------------------------------------------------------

#ifndef MAT_H
#define MAT_H

#include <algorithm>
#include <cassert>
#include <vector>

class Mat
{
  public:
    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------
    Mat() = default;
    explicit Mat(int m, int n) : n_rows_(m), n_cols_(n), data_(m * n, 0.0) {}

    // -------------------------------------------------------------------------
    // Access
    // -------------------------------------------------------------------------
    int rows() const { return n_rows_; }
    int cols() const { return n_cols_; }

    double operator()(int i, int j) const
    {
      assert(i >= 0 && i < rows());
      assert(j >= 0 && j < cols());
      return data_[i * n_cols_ + j];
    }

    // -------------------------------------------------------------------------
    // Modification
    // -------------------------------------------------------------------------
    double& operator()(int i, int j)
    {
      assert(i >= 0 && i < rows());
      assert(j >= 0 && j < cols());
      return data_[i * n_cols_ + j];
    }

    void fill(double value)
    {
      std::fill(data_.begin(), data_.end(), value);
    }

    // -------------------------------------------------------------------------
    // Operations
    // -------------------------------------------------------------------------
    
  private:

    // -------------------------------------------------------------------------
    // Data
    // -------------------------------------------------------------------------
    int n_rows_ = 0;
    int n_cols_ = 0;
    std::vector<double> data_;
};

#endif