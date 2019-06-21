enable_testing()

add_library(Sorting INTERFACE)
target_include_directories(Sorting INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

add_subdirectory(tests)