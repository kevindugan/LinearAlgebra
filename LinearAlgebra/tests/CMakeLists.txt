include(Nucleus)
enable_testing()

set(libraries LinearAlgebra)

add_mpi_test(parallel_vector 3 "${libraries}")
