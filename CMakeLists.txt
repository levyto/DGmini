cmake_minimum_required(VERSION 3.16)
project(dgmini LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(yaml-cpp REQUIRED)
find_package(muparser REQUIRED)

add_executable(dgmini 
  source/main.cpp

  source/Algebra/linalg.cpp

  source/FEM/basis1d.cpp
  source/FEM/fespace1d.cpp
  source/FEM/mass_matrix.cpp
  source/FEM/quadrature1d.cpp

  source/IO/cli_parser.cpp
  source/IO/config_reader.cpp
  source/IO/output.cpp

  source/Mesh/element1d.cpp
  source/Mesh/mesh1d.cpp

  source/Spatial/residual.cpp

  source/Temporal/cfl_number.cpp

  source/Solver/solver.cpp
  source/Solver/runtime_factory.cpp
)
target_include_directories(dgmini PRIVATE source)
target_link_libraries(dgmini PRIVATE yaml-cpp muparser)

add_executable(dgmini_unittest
  tests/unit/unittest_main.cpp
  tests/unit/basis1d_ut.cpp
  tests/unit/cli_parser_ut.cpp
  tests/unit/config_reader_ut.cpp
  tests/unit/cfl_number_ut.cpp
  tests/unit/element1d_ut.cpp
  tests/unit/expression_function_ut.cpp
  tests/unit/fespace1d_ut.cpp
  tests/unit/l2_projection_ut.cpp
  tests/unit/linalg_ut.cpp
  tests/unit/mass_matrix_ut.cpp
  tests/unit/Mat_ut.cpp
  tests/unit/mesh1d_ut.cpp
  tests/unit/modal_vector_ut.cpp
  tests/unit/numerical_flux_ut.cpp
  tests/unit/output_ut.cpp
  tests/unit/pde_ut.cpp
  tests/unit/quadrature1d_ut.cpp
  tests/unit/residual_ut.cpp
  tests/unit/runtime_factory_ut.cpp
  tests/unit/time_integrator_ut.cpp
  tests/unit/time_step_controller_ut.cpp
  tests/unit/Vec_ut.cpp

  source/Algebra/linalg.cpp

  source/FEM/basis1d.cpp
  source/FEM/fespace1d.cpp
  source/FEM/mass_matrix.cpp
  source/FEM/quadrature1d.cpp

  source/IO/cli_parser.cpp
  source/IO/config_reader.cpp
  source/IO/output.cpp

  source/Mesh/element1d.cpp
  source/Mesh/mesh1d.cpp

  source/Spatial/residual.cpp

  source/Temporal/cfl_number.cpp

  source/Solver/runtime_factory.cpp
)
target_include_directories(dgmini_unittest PRIVATE source tests/unit)
target_link_libraries(dgmini_unittest PRIVATE yaml-cpp muparser)
