
# These build targets are for the coefficient computations (code generation).
# The transformation approximations themselves are header-only. Include a
# coefficient set plus the transformation headers (e.g. via
# ecef2geo_default.hpp, or the nilsson.hpp pattern in a consuming project).

# 3.12 is required by list(TRANSFORM ...) (used by copy_insource_coefficients).
cmake_minimum_required(VERSION 3.12)
project(ecef2geo_polynomial CXX)

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
        symengine
        gmp
        )

add_executable(main_fourier
        coefficient_computations/main_fourier.cpp
        )
target_link_libraries(main_fourier
        mpfr
        ptecoef  # From the https://github.com/jnil02/point2ellipse_series repo.
        symengine
        gmp
        )

add_executable(main_taylor
        coefficient_computations/main_taylor.cpp
        )
target_link_libraries(main_taylor
        mpfr
        )

add_executable(main_additive_inv
        coefficient_computations/main_additive_inv.cpp
        )
target_link_libraries(main_additive_inv
        mpfr
        ptecoef  # From the https://github.com/jnil02/point2ellipse_series repo.
        symengine
        gmp
        sollya
        )

# ===========================================================================
# Operating-range / degree-limit configuration for coefficient generation.
#
# The approximations are optimized for a SPECIFIC altitude range -- there is no
# universal set. Override these per range when generating, e.g.:
#   cmake -DECEF2GEO_ALT_LO_LIMIT=-10000 -DECEF2GEO_ALT_HI_LIMIT=500000 ...
# A blank value keeps the header default (the "common" range) in settings.hpp.
# ===========================================================================
set(ECEF2GEO_ALT_REFERENCE "" CACHE STRING "Altitude reference h_0 [m from Earth center] (blank = default 0.0)")
set(ECEF2GEO_ALT_LO_LIMIT  "" CACHE STRING "Lower altitude limit [m] (blank = default -5000)")
set(ECEF2GEO_ALT_HI_LIMIT  "" CACHE STRING "Upper altitude limit [m] (blank = default 100000)")
set(ECEF2GEO_N_MAX "" CACHE STRING "Max N index (blank = default 8)")
set(ECEF2GEO_M_MAX "" CACHE STRING "Max M index (blank = default 6)")
set(ECEF2GEO_L_MAX "" CACHE STRING "Max L index (blank = default 8)")

set(_ecef2geo_codegen_targets
        main_trigonometric main_multiplicative main_additive
        main_fourier main_taylor main_additive_inv)
foreach(_t ${_ecef2geo_codegen_targets})
    foreach(_v ALT_REFERENCE ALT_LO_LIMIT ALT_HI_LIMIT N_MAX M_MAX L_MAX)
        if(NOT "${ECEF2GEO_${_v}}" STREQUAL "")
            target_compile_definitions(${_t} PRIVATE "ECEF2GEO_${_v}=${ECEF2GEO_${_v}}")
        endif()
    endforeach()
endforeach()

# ===========================================================================
# Installable, header-only library of the generic transformations.
#
# This installs ONLY the generic transformation headers. A user installs this
# normally (e.g. system-wide) and supplies its own range-specific coefficients
# in-tree (see the generate_coefficients / copy_insource_coefficients targets
# below and the nilsson.hpp pattern in the ecef2geo_benchmark project).
# ===========================================================================
include(GNUInstallDirs)
add_library(ecef2geo_polynomial INTERFACE)
add_library(ecef2geo_polynomial::ecef2geo_polynomial ALIAS ecef2geo_polynomial)
target_include_directories(ecef2geo_polynomial INTERFACE
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# ecef2geo_reference.hpp is only for internal usage is is not included.
set(_ecef2geo_transformation_headers
        transformation_implementations/coordinate_structs.hpp
        transformation_implementations/ecef2geo_trigonometric.hpp
        transformation_implementations/ecef2geo_polynomial.hpp
        transformation_implementations/ecef2geo_fourier.hpp)
install(FILES ${_ecef2geo_transformation_headers}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ecef2geo)
install(TARGETS ecef2geo_polynomial EXPORT ecef2geo_polynomialTargets)
install(EXPORT ecef2geo_polynomialTargets
        NAMESPACE ecef2geo_polynomial::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ecef2geo_polynomial
        FILE ecef2geo_polynomialConfig.cmake)

# ===========================================================================
# Range-specific coefficient targets.
#
# Both are run FROM THIS PROJECT, pointing at a destination directory. The
# destination defaults to this repo's own in-source generated/ folder but is
# normally set to a consuming project, e.g.:
#   cmake -DECEF2GEO_COEF_DEST=/path/to/ecef2geo/benchmark/ecef_to_geodetic/nilsson ...
#   cmake --build . --target generate_coefficients       # optimize for the configured range
#   cmake --build . --target copy_insource_coefficients  # reuse the in-source set, no optimizer
# Choosing between them is deliberate: there is no universal range.
#
# The in-source generated/ folder ships with a COMMITTED coefficient set for the
# default ("common") range. This is intentional. It lets the library be built,
# tested and used WITHOUT getting the (heavy) coefficient optimizer toolchain to
# work. By design the committed set is just the output of the default-range
# generation, so running generate_coefficients with the default settings (range
# AND destination) simply reproduces those same committed files in place.
#
# Consequently the in-source set is "the committed default-range set UNLESS it
# has been regenerated in place for a DIFFERENT range" (i.e. with ECEF2GEO_ALT_*
# overridden but the destination left at the default). copy_insource_coefficients
# copies whatever is there now.
# ===========================================================================
set(ECEF2GEO_COEF_DEST
        "${CMAKE_CURRENT_SOURCE_DIR}/transformation_implementations/generated"
        CACHE PATH "Destination directory for range-specific coefficient headers (default: the in-source generated/ folder)")

# Optimize a fresh coefficient set for the configured range into ECEF2GEO_COEF_DEST.
# With the default destination, this overwrites the in-source generated/ folder.
add_custom_target(generate_coefficients
        COMMAND ${CMAKE_COMMAND} -E make_directory "${ECEF2GEO_COEF_DEST}"
        COMMAND ${CMAKE_COMMAND} -E echo "Optimizing a fresh coefficient set into ${ECEF2GEO_COEF_DEST} for the configured range (ECEF2GEO_ALT_* cache variables). With the default destination this overwrites the in-source generated/ folder."
        COMMAND $<TARGET_FILE:main_trigonometric> "${ECEF2GEO_COEF_DEST}"
        COMMAND $<TARGET_FILE:main_multiplicative> "${ECEF2GEO_COEF_DEST}"
        COMMAND $<TARGET_FILE:main_additive> "${ECEF2GEO_COEF_DEST}"
        COMMAND $<TARGET_FILE:main_taylor> "${ECEF2GEO_COEF_DEST}"
        COMMAND $<TARGET_FILE:main_fourier> "${ECEF2GEO_COEF_DEST}"
        COMMAND $<TARGET_FILE:main_additive_inv> "${ECEF2GEO_COEF_DEST}"
        DEPENDS ${_ecef2geo_codegen_targets}
        VERBATIM
        USES_TERMINAL)

# The canonical set of coefficient headers a current generation produces (i.e.
# what generate_coefficients writes). Listed explicitly so the copy below does
# not drag along stale/old artifacts that may sit in the generated/ folder.
set(_ecef2geo_coef_headers
        xi.hpp chi.hpp                      # main_trigonometric
        sigma.hpp tau.hpp                   # main_multiplicative
        omega.hpp mu.hpp                    # main_additive
        ptepoly_sigma.hpp ptepoly_tau.hpp   # main_taylor
        ptepoly_phi.hpp ptepoly_h.hpp ptepoly_cos.hpp ptepoly_sin.hpp  # main_fourier
        phi.hpp h.hpp)                      # main_additive_inv
list(TRANSFORM _ecef2geo_coef_headers PREPEND
        "${CMAKE_CURRENT_SOURCE_DIR}/transformation_implementations/generated/"
        OUTPUT_VARIABLE _ecef2geo_coef_header_paths)

# Reuse the coefficient set currently in the in-source generated/ folder (no
# optimizer run). That is the committed common-range set UNLESS generate_coefficients
# has regenerated it in place -- this copies whatever is there now.
add_custom_target(copy_insource_coefficients
        COMMAND ${CMAKE_COMMAND} -E make_directory "${ECEF2GEO_COEF_DEST}"
        COMMAND ${CMAKE_COMMAND} -E echo "Copying the in-source coefficient set from transformation_implementations/generated/ into ${ECEF2GEO_COEF_DEST} (no optimizer run). This is the committed common-range set unless generate_coefficients has overwritten it in place. It is NOT universal -- generate a dedicated set for any other range."
        COMMAND ${CMAKE_COMMAND} -E copy ${_ecef2geo_coef_header_paths}
                "${ECEF2GEO_COEF_DEST}/"
        VERBATIM)
