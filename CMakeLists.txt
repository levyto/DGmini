cmake_minimum_required(VERSION 3.16)
project(dgmini LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(dgmini 
  source/main.cpp
  source/basis1d.cpp
  source/element1d.cpp
  source/fespace1d.cpp
  source/input.cpp
  source/linalg.cpp
  source/mass_matrix.cpp
  source/mesh1d.cpp
  source/output.cpp
  source/quadrature1d.cpp
  source/stiffness_matrix.cpp
)
target_include_directories(dgmini PRIVATE source)

add_executable(dgmini_unittest
  unittest/unittest_main.cpp
  unittest/basis1d_ut.cpp
  unittest/element1d_ut.cpp
  unittest/fespace1d_ut.cpp
  unittest/linalg_ut.cpp
  unittest/mass_matrix_ut.cpp
  unittest/Mat_ut.cpp
  unittest/mesh1d_ut.cpp
  unittest/quadrature1d_ut.cpp
  unittest/stiffness_matrix_ut.cpp
  unittest/Vec_ut.cpp
  source/basis1d.cpp
  source/element1d.cpp
  source/fespace1d.cpp
  source/linalg.cpp
  source/mass_matrix.cpp
  source/mesh1d.cpp
  source/quadrature1d.cpp
  source/stiffness_matrix.cpp
)
target_include_directories(dgmini_unittest PRIVATE source unittest)