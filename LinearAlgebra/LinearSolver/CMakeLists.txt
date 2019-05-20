
enable_testing()

add_library(LinearSolver GaussianElimination.cpp)

target_include_directories(LinearSolver INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

target_link_libraries(LinearSolver Matrix)

add_subdirectory(tests)