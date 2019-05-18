include(Nucleus)
enable_testing()

set(libraries Matrix)

add_mpi_test(parallel_matrix_blockRowPartition 5 "${libraries}")
