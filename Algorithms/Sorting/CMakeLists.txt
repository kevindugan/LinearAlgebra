enable_testing()

add_library(Sorting INTERFACE)
target_include_directories(Sorting INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(Sorting INTERFACE Parallel)

add_subdirectory(tests)