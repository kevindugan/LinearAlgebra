include(Nucleus)
enable_testing()

set(libraries Vector)

add_mpi_test(parallel_vector 3 "${libraries}")
