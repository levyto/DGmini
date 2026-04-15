cmake_minimum_required(VERSION 3.16)
project(dgmini LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(dgmini 
  source/main.cpp

  source/Algebra/linalg.cpp

  source/FEM/basis1d.cpp
  source/FEM/fespace1d.cpp
  source/FEM/mass_matrix.cpp
  source/FEM/quadrature1d.cpp
  source/FEM/stiffness_matrix.cpp

  source/IO/input.cpp
  source/IO/output.cpp

  source/Mesh/element1d.cpp
  source/Mesh/mesh1d.cpp

  source/Spatial/l2_projection.cpp
  source/Spatial/residual.cpp

  source/Temporal/cfl_number.cpp
)
target_include_directories(dgmini PRIVATE source)

add_executable(dgmini_unittest
  unittest/unittest_main.cpp
  unittest/basis1d_ut.cpp
  unittest/cfl_number_ut.cpp
  unittest/element1d_ut.cpp
  unittest/fespace1d_ut.cpp
  unittest/l2_projection_ut.cpp
  unittest/linalg_ut.cpp
  unittest/mass_matrix_ut.cpp
  unittest/Mat_ut.cpp
  unittest/mesh1d_ut.cpp
  unittest/modal_vector_ut.cpp
  unittest/numerical_flux_ut.cpp
  unittest/output_ut.cpp
  unittest/pde_ut.cpp
  unittest/quadrature1d_ut.cpp
  unittest/residual_ut.cpp
  unittest/stiffness_matrix_ut.cpp
  unittest/time_integrator_ut.cpp
  unittest/time_step_controller_ut.cpp
  unittest/Vec_ut.cpp
  
  source/Algebra/linalg.cpp

  source/FEM/basis1d.cpp
  source/FEM/fespace1d.cpp
  source/FEM/mass_matrix.cpp
  source/FEM/quadrature1d.cpp
  source/FEM/stiffness_matrix.cpp

  source/IO/output.cpp

  source/Mesh/element1d.cpp
  source/Mesh/mesh1d.cpp

  source/Spatial/l2_projection.cpp
  source/Spatial/residual.cpp

  source/Temporal/cfl_number.cpp
)
target_include_directories(dgmini_unittest PRIVATE source unittest)