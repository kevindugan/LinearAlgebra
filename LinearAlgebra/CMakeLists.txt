
enable_testing()

add_library(LinearAlgebra LinearAlgebra.cpp
                          Vector.cpp)

add_subdirectory(tests)