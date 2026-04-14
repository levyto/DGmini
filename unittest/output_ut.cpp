// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for output functions  
//
// -----------------------------------------------------------------------------

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "output_ut.h"
#include "Mesh/mesh1d.h"
#include "Spatial/modal_vector.h"

// -----------------------------------------------------------------------------
// Description: Helper functions for output UTs
// -----------------------------------------------------------------------------
namespace
{
  // ---------------------------------------------------------------------------
  // Description: Read all lines from a text file into a vector of strings
  // ---------------------------------------------------------------------------
  std::vector<std::string> readAllLines(const std::string& filename)
  {
    std::ifstream in(filename);
    if (!in.is_open())
    {
      throw std::runtime_error("Failed to open file in unit test: " + filename);
    }

    std::vector<std::string> lines;
    std::string line;

    while (std::getline(in, line))
    {
      lines.push_back(line);
    }

    return lines;
  }

  // ---------------------------------------------------------------------------
  // Description: Check if a file exists
  // ---------------------------------------------------------------------------
  bool fileExists(const std::string& filename)
  {
    return std::filesystem::exists(filename);
  }

  // ---------------------------------------------------------------------------
  // Description: Remove a file if it exists
  // ---------------------------------------------------------------------------
  void removeFileIfExists(const std::string& filename)
  {
    if (std::filesystem::exists(filename))
    {
      std::filesystem::remove(filename);
    }
  }

  // ---------------------------------------------------------------------------
  // Description: Remove a directory if it exists
  // ---------------------------------------------------------------------------
  void removeDirectoryIfExists(const std::string& dirname)
  {
    if (std::filesystem::exists(dirname))
    {
      std::filesystem::remove_all(dirname);
    }
  }
}

// -----------------------------------------------------------------------------
// Description: Check that the header of the output file from 
//              writeModalSolution1D is correct
// -----------------------------------------------------------------------------
void Test_output_writeModalSolution1DHeader()
{
  const std::string filename = "output/unittest_modal_header.dat";
  removeFileIfExists(filename);

  const int Ne = 2;
  const int p = 2;
  Mesh1D mesh(0.0, 1.0, Ne);
  ModalVector solution(mesh.Ne(), p + 1);

  solution.zero();

  const double time = 0.125;
  writeModalSolution1D(filename, mesh, solution, time, false);

  const std::vector<std::string> lines = readAllLines(filename);

  Check(lines.size() >= 5,                        "Modal output header has too few lines");
  Check(lines[0] == "# DGmini output",            "Wrong first header line");
  Check(lines[1] == "# time 0.125",               "Wrong time header line");
  Check(lines[2] == "# Ne 2",                     "Wrong Ne header line");
  Check(lines[3] == "# p 2",                      "Wrong p header line");
  Check(lines[4].find("# columns: element") == 0, "Wrong columns header line");
}

// -----------------------------------------------------------------------------
// Description: Check that the data of the output file from 
//              writeModalSolution1D is correct
// -----------------------------------------------------------------------------
void Test_output_writeModalSolution1DData()
{
  const std::string filename = "output/unittest_modal_data.dat";
  removeFileIfExists(filename);

  const int Ne = 2;
  const int p = 1;
  Mesh1D mesh(0.0, 1.0, Ne);
  ModalVector solution(mesh.Ne(), p + 1);

  solution(0,0) = 1.0;
  solution(0,1) = 2.0;
  solution(1,0) = 3.0;
  solution(1,1) = 4.0;

  writeModalSolution1D(filename, mesh, solution, 0.0, false);

  const std::vector<std::string> lines = readAllLines(filename);

  Check(lines.size(), "Unexpected number of lines in modal output");

  std::istringstream row0(lines[5]);
  std::istringstream row1(lines[6]);

  int e;
  double xl, xr, c0, c1;

  row0 >> e >> xl >> xr >> c0 >> c1;
  Check(e == 0, "Wrong element id in first data row");
  CheckEqual(xl, 0.0, 1e-14, "Wrong x_left in first data row");
  CheckEqual(xr, 0.5, 1e-14, "Wrong x_right in first data row");
  CheckEqual(c0, 1.0, 1e-14, "Wrong c0 in first data row");
  CheckEqual(c1, 2.0, 1e-14, "Wrong c1 in first data row");

  row1 >> e >> xl >> xr >> c0 >> c1;
  Check(e == 1, "Wrong element id in second data row");
  CheckEqual(xl, 0.5, 1e-14, "Wrong x_left in second data row");
  CheckEqual(xr, 1.0, 1e-14, "Wrong x_right in second data row");
  CheckEqual(c0, 3.0, 1e-14, "Wrong c0 in second data row");
  CheckEqual(c1, 4.0, 1e-14, "Wrong c1 in second data row");
}

// -----------------------------------------------------------------------------
// Description: Check that the output file from writeSolution1D has correct 
//              header and data
// -----------------------------------------------------------------------------
void Test_output_writeSolution1D()
{
  const std::string filename = "output/unittest_pointwise.dat";
  removeFileIfExists(filename);

  Vec x(3);
  Vec u(3);

  x[0] = 0.0;  x[1] = 0.5;  x[2] = 1.0;
  u[0] = 1.0;  u[1] = 2.0;  u[2] = 3.0;

  writeSolution1D(filename, x, u, 0.25, false);

  const std::vector<std::string> lines = readAllLines(filename);

  Check(lines.size(),                     "Unexpected number of lines in pointwise output");
  Check(lines[0] == "# DGmini output",    "Wrong first header line in pointwise output");
  Check(lines[1] == "# time 0.25",        "Wrong time header line in pointwise output");
  Check(lines[2].find("# columns:") == 0, "Wrong columns header line in pointwise output");

  std::istringstream row0(lines[3]);
  std::istringstream row2(lines[5]);

  double xv, uv;

  row0 >> xv >> uv;
  CheckEqual(xv, 0.0, 1e-14, "Wrong x in first pointwise row");
  CheckEqual(uv, 1.0, 1e-14, "Wrong u in first pointwise row");

  row2 >> xv >> uv;
  CheckEqual(xv, 1.0, 1e-14, "Wrong x in last pointwise row");
  CheckEqual(uv, 3.0, 1e-14, "Wrong u in last pointwise row");
}

// -----------------------------------------------------------------------------
// Description: Check that TimeSeriesWriter writes output files at correct times 
//              with correct names
// -----------------------------------------------------------------------------
void Test_output_TimeSeriesWriterWriteSequence()
{
  const std::string directory = "output/unittest_timeseries";
  removeDirectoryIfExists(directory);

  const int Ne = 2;
  const int p = 1;
  Mesh1D mesh(0.0, 1.0, Ne);
  ModalVector solution(mesh.Ne(), p+1);

  solution.zero();

  TimeSeriesWriter writer(directory, "solution", 0.1);

  writer.write(mesh, solution, 0.0, false);   // should write 0000
  writer.write(mesh, solution, 0.05, false);  // should not write
  writer.write(mesh, solution, 0.10, false);  // should write 0001
  writer.write(mesh, solution, 0.19, false);  // should not write
  writer.write(mesh, solution, 0.20, false);  // should write 0002
  writer.writeFinal(mesh, solution, 0.23, false);

  Check(fileExists(directory + "/solution_0000.dat"), "Missing solution_0000.dat");
  Check(fileExists(directory + "/solution_0001.dat"), "Missing solution_0001.dat");
  Check(fileExists(directory + "/solution_0002.dat"), "Missing solution_0002.dat");
  Check(!fileExists(directory + "/solution_0003.dat"), "Unexpected solution_0003.dat was written");
  Check(fileExists(directory + "/solution_final.dat"), "Missing solution_final.dat");

  const std::vector<std::string> final_lines = readAllLines(directory + "/solution_final.dat");

  Check(final_lines[1] == "# time 0.23", "Wrong final time written to solution_final.dat");
}
