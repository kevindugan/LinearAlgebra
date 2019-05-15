
add_library(Nucleus Nucleus.cpp)

target_link_libraries(Nucleus gmock_main)
target_include_directories(Nucleus INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")
