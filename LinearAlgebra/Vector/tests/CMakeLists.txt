include(Nucleus)
enable_testing()

set(libraries Vector)

add_mpi_test(parallel_vector_blockPartition 3 "${libraries}")
