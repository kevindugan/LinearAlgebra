include(Nucleus)
enable_testing()

set(libraries Parallel)

add_mpi_test(parallel_blockPartition 5 "${libraries}")
add_mpi_test(parallel_cyclePartition 5 "${libraries}")