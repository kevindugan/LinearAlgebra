include(Nucleus)
enable_testing()

set(libraries LinearSolver)

add_mpi_test(gaussian_elimination 5 "${libraries}")

add_executable(time_GE_cycle_partition timing_gaussian_elimination_cycle.cpp)
target_link_libraries(time_GE_cycle_partition LinearSolver)