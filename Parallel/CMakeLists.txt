enable_testing()

add_library(Parallel Parallel.cpp)
target_include_directories(Parallel INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")