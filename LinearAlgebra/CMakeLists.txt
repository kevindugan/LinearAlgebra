
enable_testing()

add_library(LinearAlgebra LinearAlgebra.cpp
                          Vector.cpp
                          Matrix.cpp)
target_include_directories(LinearAlgebra INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

add_subdirectory(tests)
