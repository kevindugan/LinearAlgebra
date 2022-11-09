
enable_testing()

add_library(Matrix Matrix_BlockRowPartition.cpp
                   Matrix_CycleRowPartition.cpp
                   SparseMatrix_BlockRowPartition.cpp)
target_include_directories(Matrix INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

target_link_libraries(Matrix Vector)

add_subdirectory(tests)