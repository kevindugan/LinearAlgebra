
enable_testing()

add_library(Vector Vector.cpp)

target_include_directories(Vector INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

target_link_libraries(Vector LinearAlgebra)

add_subdirectory(tests)