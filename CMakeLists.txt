cmake_minimum_required(VERSION 3.6)

project(LinearAlgebra CXX)

#####################################
# GTest
#####################################
enable_testing()
add_subdirectory(contrib/gtest)

#####################################
# MPI
#####################################
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_C_COMPILER mpicc)

#####################################
# Project
#####################################
add_subdirectory(LinearAlgebra)

add_executable(test_exec main.cpp)
target_link_libraries(test_exec LinearAlgebra)
