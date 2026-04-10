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

#include "output.h"

// -----------------------------------------------------------------------------
// Description: Output solution u at points x to text file
// -----------------------------------------------------------------------------
void writeSolution1D(const std::string& filename, const Vec& x, const Vec& u)
{
  assert(x.size() == u.size());

  std::filesystem::create_directories("output");
  std::ofstream out("output/" + filename);
  assert(out.is_open());

  out << std::setprecision(16);

  out << "# DGmini output\n";
  out << "# columns: \tx\tsolution_value\n";

  for (int i = 0; i < x.size(); ++i)
  {
    out << x[i] << " " << u[i] << "\n";
  }

  std::cout << "\nPointwise solution written to output/" << filename << "\n\n";
}

// -----------------------------------------------------------------------------
// Description: Write element-wise modal coefficients to a text file
//
// Format:
//   # DGmini output
//   # Ne <value>
//   # p <value>
//   # columns: element x_left x_right c0 c1 ... c_{p}
//   0 <x_left> <x_right> <c0> <c1> ... <c_p>
// -----------------------------------------------------------------------------
void writeModalSolution1D(const std::string& filename,
                          const Mesh1D& mesh,
                          const std::vector<Vec>& coefficients)
{
  assert(static_cast<int>(coefficients.size()) == mesh.Ne());

  const int n_coeffs = coefficients[0].size();

  for (int e = 0; e < mesh.Ne(); ++e)
  {
    assert(coefficients[e].size() == n_coeffs);
  }

  std::filesystem::create_directories("output");
  std::ofstream out("output/" + filename);
  assert(out.is_open());

  out << std::setprecision(16);

  out << "# DGmini output\n";
  out << "# Ne \t" << mesh.Ne()  << "\n";
  out << "# p  \t"  << n_coeffs-1 << "\n";
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
      out << "\t" << coefficients[e][j];

    out << "\n";
  }

  std::cout << "\nModal solution written to output/" << filename << "\n\n";
}