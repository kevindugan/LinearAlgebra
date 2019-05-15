
enable_testing()

function(add_mpi_test test_name num_procs)
    message(STATUS "Adding Test: ${test_name}")
    add_executable(${test_name} test_${test_name}.cpp)
    target_link_libraries(${test_name} LinearAlgebra gmock_main)
    set(test_parameters ${MPIEXEC_NUMPROC_FLAG} ${num_procs} "./${test_name}")
    add_test(NAME ${test_name} COMMAND ${MPIEXEC} ${test_parameters})
endfunction(add_mpi_test)

add_mpi_test(parallel_vector 3)

# add_executable(test_parallel_vector test_parallel_vector.cpp)
# target_link_libraries(test_parallel_vector LinearAlgebra gtest_main)
# add_test(ParallelVector test_parallel_vector)
