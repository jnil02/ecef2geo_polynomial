
# Catch2 unit test build target.
# For testing reference transformations which the polynomial transformation
# generations are based on.

cmake_minimum_required(VERSION 3.12)
project(transformation_tests CXX)

set(CMAKE_CXX_STANDARD 17)  # At least 17 required.

add_executable(transformation_tests
        transformation_tests.cpp
        )

target_link_libraries(transformation_tests
        mpfr
        )