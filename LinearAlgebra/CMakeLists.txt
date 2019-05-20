
enable_testing()

add_library(LinearAlgebra LinearAlgebra.cpp)
target_include_directories(LinearAlgebra INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

add_subdirectory(Vector)
add_subdirectory(Matrix)
add_subdirectory(LinearSolver)