
# Catch2 unit test build target.
# For testing reference transformations which the polynomial transformation
# generations are based on.

# Not sure what the minimum is. A check is required to avoid a warning. Any
# reasonably recent cmake should do.
cmake_minimum_required(VERSION 3.0)
project(transformation_tests)

add_executable(transformation_tests
        transformation_tests.cpp
        )

target_link_libraries(transformation_tests
        mpfr
        )