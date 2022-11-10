enable_testing()

add_library(Parallel Parallel.cpp
                     BlockPartition1D.cpp
                     CyclePartition1D.cpp)
target_include_directories(Parallel INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

add_subdirectory(tests)