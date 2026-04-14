// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Definitions of outputs
//
// -----------------------------------------------------------------------------

#ifndef OUTPUT_H
#define OUTPUT_H

#include <filesystem>
#include <string>
#include <vector>

#include "Algebra/Vec.h"
#include "Mesh/mesh1d.h"
#include "Spatial/modal_vector.h"

void writeModalSolution1D
(
  const std::string& filename,
  const Mesh1D& mesh,
  const ModalVector& solution,
  const double time,
  bool verbose = true
);

void writeSolution1D
(
  const std::string& filename,
  const Vec& x, 
  const Vec& solution,
  const double time,
  bool verbose = true
);

class TimeSeriesWriter
{
  public:
    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------
    TimeSeriesWriter
    (
      const std::string& directory,
      const std::string& prefix,
      double output_dt
    )
    : directory_(directory),
      prefix_(prefix),
      output_dt_(output_dt),
      next_output_time_(0.0),
      output_id_(0)
    {
      std::filesystem::create_directories(directory_);
    }

    // -------------------------------------------------------------------------
    // Access
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Modification
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Operations
    // -------------------------------------------------------------------------
    bool writeNow(double time) const { return time >= next_output_time_; }

    void write
    (
      const Mesh1D& mesh, 
      const ModalVector& solution, 
      double time,
    bool verbose = true
    );
    
    void writeFinal
    (
      const Mesh1D& mesh, 
      const ModalVector& solution, 
      double time, 
      bool verbose = true
    );

  private:
    // -------------------------------------------------------------------------
    // Data
    // -------------------------------------------------------------------------
    std::string directory_;
    std::string prefix_;

    double output_dt_ = 0.0;
    double next_output_time_ = 0.0;

    int output_id_ = 0;
};

#endif

