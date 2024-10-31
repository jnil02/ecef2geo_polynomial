
# These build targets are for the coefficien computations (code generation).
# The transformation approximation are used by simply including the header
# file transformation_implementations/ecef2geo.hpp.

# Not sure what the minimum is. A check is required to avoid a warning. Any
# reasonably recent cmake should do.
cmake_minimum_required(VERSION 3.0)
project(ecef2geo_polynomial)

set(CMAKE_CXX_STANDARD 23)  # At least 17 required.
set(CMAKE_CXX_FLAGS -O3)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

add_executable(main_trigonometric
        coefficient_computations/main_trigonometric.cpp
        )
target_link_libraries(main_trigonometric
        mpfr
        sollya
        )

add_executable(main_multiplicative
        coefficient_computations/main_multiplicative.cpp
        )
target_link_libraries(main_multiplicative
        mpfr
        mplapack_mpfr
        mpblas_mpfr
        sollya
        )

add_executable(main_additive
        coefficient_computations/main_additive.cpp
        )
target_link_libraries(main_additive
        mpfr
        mplapack_mpfr
        mpblas_mpfr
        tbb  # For std::execution::par
        ginac  # Symbolic expressions.
        cln  # Numbers in ginac.
        )
