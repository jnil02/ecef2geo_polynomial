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
constexpr double ALT_REFERENCE = 0.0;  // I.e. h_0. Meters from the center of the earth.
constexpr double ALT_LO_LIMIT = -5000;  // Meter above sea level.
constexpr double ALT_HI_LIMIT = 100000;  // Meter above sea level.

// Where to place the generated files.
const std::string CODEGEN_FOLDER = "./../transformation_implementations/generated/";

// Polynomials are generated for up to these index limits.
constexpr int N_MAX = 8;
constexpr int M_MAX = 6;
constexpr int L_MAX = 8;

// Polynomial evaluation schema used throughout.
constexpr int (*POLY_EVAL_SCHEMA)(int, int) = evaluation_schema::estrin_schema;

#endif // ECEF2GEO_POLYNOMIAL_SETTINGS_HPP
