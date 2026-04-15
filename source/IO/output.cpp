// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Definitions of outputs
//
// -----------------------------------------------------------------------------

#include <cassert>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "IO/output.h"

// -----------------------------------------------------------------------------
// Description: Write element-wise modal coefficients to a text file
//
// Format:
//   # DGmini output
//   # time = <value>
//   # Ne <value>
//   # p <value>
//   # columns: element x_left x_right c0 c1 ... c_{p}
//   0 <x_left> <x_right> <c0> <c1> ... <c_p>
// -----------------------------------------------------------------------------
void writeModalSolution1D
(
  const std::string& filename,
  const Mesh1D& mesh,
  const ModalVector& solution,
  const double time,
  bool verbose
)
{
  assert(static_cast<int>(solution.Ne()) == mesh.Ne());

  const int n_coeffs = solution.localDoFs();

  std::filesystem::path path(filename);
  std::filesystem::path directory = path.parent_path();

  if (!directory.empty())
  {
    std::filesystem::create_directories(directory);
  }

  std::ofstream out(filename);
  assert(out.is_open());

  out << std::setprecision(16);

  out << "# DGmini output\n";
  out << "# time " << time << "\n";
  out << "# Ne "   << mesh.Ne()  << "\n";
  out << "# p "    << n_coeffs-1 << "\n";
  out << "# columns: element\tx_left\tx_right";
  for (int j = 0; j < n_coeffs; ++j)
    out << "\tc" << j;
  out << "\n";

  for (int e = 0; e < mesh.Ne(); ++e)
  {
    const Element1D& element = mesh.element(e);

    out << e                << "\t"
        << element.left()   << "\t"
        << element.right();

    for (int j = 0; j < n_coeffs; ++j)
      out << "\t" << solution(e,j);

    out << "\n";
  }

  if (verbose)
  {
    std::cout << "\nModal solution written to " << filename;
  }
}

// -----------------------------------------------------------------------------
// Description: Output solution u at points x to text file
// -----------------------------------------------------------------------------
void writeSolution1D
(
  const std::string& filename, 
  const Vec& x, 
  const Vec& solution, 
  const double time,
  bool verbose
)
{
  assert(x.size() == solution.size());

  std::filesystem::path path(filename);
  std::filesystem::path directory = path.parent_path();
  
  if (!directory.empty())
  {
    std::filesystem::create_directories(directory);
  }

  std::ofstream out(filename);
  assert(out.is_open());

  out << std::setprecision(16);

  out << "# DGmini output\n";
  out << "# time " << time << "\n";
  out << "# columns: \tx\tsolution_value\n";

  for (int i = 0; i < x.size(); ++i)
  {
    out << x[i] << " " << solution[i] << "\n";
  }

  if (verbose)
  {
    std::cout << "\nPointwise solution written to " << filename << "\n\n";
  }
}

// -----------------------------------------------------------------------------
// Description: Write time series of solutions to files with names 
//              <prefix>_0000.dat, <prefix>_0001.dat, etc.
// -----------------------------------------------------------------------------
void TimeSeriesWriter::write
(
  const Mesh1D& mesh, 
  const ModalVector& solution, 
  double time,
  bool verbose
)
{
  if (!writeNow(time))
  {
    return;
  }

  std::ostringstream name;
  name << prefix_
       << "_"
       << std::setw(4) << std::setfill('0') << output_id_
       << ".dat";

  writeModalSolution1D
  (
    directory_ + "/" + name.str(), 
    mesh, 
    solution, 
    time,
    verbose
  );

  ++output_id_;
  next_output_time_ += output_dt_;
}

// -----------------------------------------------------------------------------
// Description: Write final solution to file with name <prefix>_final.dat
// -----------------------------------------------------------------------------
void TimeSeriesWriter::writeFinal
(
  const Mesh1D& mesh,
  const ModalVector& solution,
  double time,
  bool verbose
)
{
  writeModalSolution1D
  (
    directory_ + "/" + prefix_ + "_final.dat", 
    mesh, 
    solution, 
    time,
    verbose
  );
}