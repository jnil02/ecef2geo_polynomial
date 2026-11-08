// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the BSD 2-Clause License.

#ifndef ECEF2GEO_POLYNOMIAL_SETTINGS_HPP
#define ECEF2GEO_POLYNOMIAL_SETTINGS_HPP

/*
 * General settings for polynomial approximation code generation.
 * Settings which are specific to one of the code generation mains are found
 * in the main-files themselves.
 */

#include "evaluation_schema.hpp"  // For polynomial evaluation schema.
#include <string>  // std::string

// Polynomial approximation altitude ranges and altitude reference point.
// Each may be overridden at configure/compile time via -D (e.g. the CMake cache
// variables ECEF2GEO_ALT_LO_LIMIT etc.) so a consumer can generate a set for
// its own operating range. The defaults below are the "common" range.
#ifndef ECEF2GEO_ALT_REFERENCE
#define ECEF2GEO_ALT_REFERENCE 0.0  // I.e. h_0. Meters from the center of the earth.
#endif
#ifndef ECEF2GEO_ALT_LO_LIMIT
#define ECEF2GEO_ALT_LO_LIMIT -5000  // Meter above sea level.
#endif
#ifndef ECEF2GEO_ALT_HI_LIMIT
#define ECEF2GEO_ALT_HI_LIMIT 100000  // Meter above sea level.
#endif
constexpr double ALT_REFERENCE = ECEF2GEO_ALT_REFERENCE;
constexpr double ALT_LO_LIMIT = ECEF2GEO_ALT_LO_LIMIT;
constexpr double ALT_HI_LIMIT = ECEF2GEO_ALT_HI_LIMIT;

// Where to place the generated files. This is the DEFAULT; it can be overridden
// at runtime by passing an output directory as the first command-line argument
// (see set_codegen_folder_from_args), which is how a generation is pointed at a
// consuming project. There is no universal range, so choosing where results go
// (and for which range) is always a deliberate, explicit choice.
#ifndef ECEF2GEO_CODEGEN_FOLDER
#define ECEF2GEO_CODEGEN_FOLDER "./../transformation_implementations/generated/"
#endif
inline std::string CODEGEN_FOLDER = ECEF2GEO_CODEGEN_FOLDER;

// Override CODEGEN_FOLDER from argv[1] when provided (a trailing '/' is added if
// missing). Call this first thing in each generator's main().
inline void set_codegen_folder_from_args(int argc, char** argv) {
	if (argc > 1) {
		CODEGEN_FOLDER = argv[1];
		if (!CODEGEN_FOLDER.empty() && CODEGEN_FOLDER.back() != '/')
			CODEGEN_FOLDER += '/';
	}
}

// Polynomials are generated for up to these index limits.
#ifndef ECEF2GEO_N_MAX
#define ECEF2GEO_N_MAX 8
#endif
#ifndef ECEF2GEO_M_MAX
#define ECEF2GEO_M_MAX 6
#endif
#ifndef ECEF2GEO_L_MAX
#define ECEF2GEO_L_MAX 8
#endif
constexpr int N_MAX = ECEF2GEO_N_MAX;
constexpr int M_MAX = ECEF2GEO_M_MAX;
constexpr int L_MAX = ECEF2GEO_L_MAX;

// Polynomial evaluation schema used throughout.
constexpr int (*POLY_EVAL_SCHEMA)(int, int) = evaluation_schema::estrin_schema;

#endif // ECEF2GEO_POLYNOMIAL_SETTINGS_HPP
