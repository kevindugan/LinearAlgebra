include(Nucleus)
enable_testing()

set(libraries LinearSolver)

add_mpi_test(gaussian_elimination 5 "${libraries}")