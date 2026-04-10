cmake_minimum_required(VERSION 3.16)
project(dgmini LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(dgmini 
  source/main.cpp
  source/linalg.cpp
)
target_include_directories(dgmini PRIVATE source)

add_executable(dgmini_unittest
    unittest/main_ut.cpp
    unittest/Vec_ut.cpp
    unittest/linalg_ut.cpp
    source/linalg.cpp
)
target_include_directories(dgmini_unittest PRIVATE source unittest)