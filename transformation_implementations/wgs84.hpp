// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the BSD 2-Clause License.

/*
 * WGS84 reference ellipsoid constants.
 */

#ifndef ECEF2GEO_WGS84_HPP
#define ECEF2GEO_WGS84_HPP

#include <cmath>  // std::sqrt (for the derived eccentricity).

namespace wgs84 {

// Primary WGS84 constants.
constexpr double a = 6378137.0;  // Semi-major axis / equatorial radius.
constexpr double f = 1.0 / 298.257223563;  // Flattening. f = (a - b) / a
// Derived WGS84 constants.
constexpr double b = a - f * a;  // Semi-minor axis / polar radius.
constexpr double e2 = 1.0 - (b * b) / (a * a);  // First eccentricity squared.
// Note, constexpr sqrt is a gcc extension. If not available you may use:
//constexpr double e = 0.08181919084262157;
// However, in this case, you will have to make sure it matches the primary
// constants a and f.
constexpr double e = std::sqrt(e2);  // First eccentricity.

}  // namespace wgs84

#endif  // ECEF2GEO_WGS84_HPP
