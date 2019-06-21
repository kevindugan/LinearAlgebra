
enable_testing()

add_library(LinearAlgebra INTERFACE)
target_include_directories(LinearAlgebra INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(LinearAlgebra INTERFACE Parallel)

add_subdirectory(Vector)
add_subdirectory(Matrix)
add_subdirectory(LinearSolver)