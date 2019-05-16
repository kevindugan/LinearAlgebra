
function(add_mpi_test test_name num_procs link_libraries)
    message(STATUS "Adding Test: ${test_name}")
    add_executable(${test_name} test_${test_name}.cpp)
    target_link_libraries(${test_name} ${link_libraries} Nucleus)
    set(test_parameters ${MPIEXEC_NUMPROC_FLAG} ${num_procs} "./${test_name}")
    add_test(NAME ${test_name} COMMAND ${MPIEXEC} ${test_parameters})
endfunction(add_mpi_test)