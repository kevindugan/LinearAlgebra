cmake_minimum_required(VERSION 3.6)

project(LinearAlgebra CXX)

#####################################
# GTest
#####################################
enable_testing()
add_subdirectory(contrib/gtest)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/Modules")

#####################################
# MPI
#####################################
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_C_COMPILER mpicc)

#####################################
# Project
#####################################
add_subdirectory(Nucleus)
add_subdirectory(Parallel)
add_subdirectory(LinearAlgebra)
add_subdirectory(Algorithms)

add_executable(test_exec main.cpp)
target_link_libraries(test_exec Vector)
