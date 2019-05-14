cmake_minimum_required(VERSION 3.6)

project(LinearAlgebra CXX)

find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_C_COMPILER mpicc)

add_subdirectory(LinearAlgebra)

add_executable(test main.cpp)
target_link_libraries(test LinearAlgebra)
