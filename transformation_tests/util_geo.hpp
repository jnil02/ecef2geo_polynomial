// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the BSD 2-Clause License.

#ifndef TRANSFORMATION_TESTS_UTIL_GEO_HPP
#define TRANSFORMATION_TESTS_UTIL_GEO_HPP

/*
 * Various ECEF to geodetic calculation utility function including wgs84
 * constants and reference ECEF to geodetic coordinate transformations.
 */

#include "./../transformation_implementations/coordinate_structs.hpp"
#include "./../transformation_implementations/ecef2geo_reference.hpp"

#define _USE_MATH_DEFINES
#include <cmath>  // sqrt, sin, cos
#include <cassert>

namespace util_geo {

/** Convert {latitude, longitude, altitude} to {x, y, z}.
 *
 * @param geo {latitude, longitude, altitude}
 * @return {x, y, z}
 */
static inline xyz lla2ecef(lla geo) {
	double sin_lat = std::sin(geo.lat);
	double cos_lat = std::cos(geo.lat);
	double sin_lon = std::sin(geo.lon);
	double cos_lon = std::cos(geo.lon);
	double N = wgs84::a / std::sqrt(1.0 - wgs84::e2 * sin_lat * sin_lat);
	double x = (N + geo.alt) * cos_lat * cos_lon;
	double y = (N + geo.alt) * cos_lat * sin_lon;
#ifndef NDEBUG
	// Check that altitude is within its defined range.
	assert((1.0 - wgs84::e2) * N + geo.alt >= 0);
#endif
	double z = ((1.0 - wgs84::e2) * N + geo.alt) * sin_lat;
	return {x, y, z};
}

/** Convert {n-vector, altitude} to {x, y, z}.
 *
 * @param geo {n-vector, altitude}
 * @return {x, y, z}
 */
static inline xyz nva2ecef(nva geo) {
	double N = wgs84::a / std::sqrt(1.0 - wgs84::e2 * geo.n.z * geo.n.z);
	double x = (N + geo.alt) * geo.n.x;
	double y = (N + geo.alt) * geo.n.y;
#ifndef NDEBUG
	// Check that altitude is within its defined range.
	assert((1.0 - wgs84::e2) * N + geo.alt >= 0);
#endif
	double z = ((1.0 - wgs84::e2) * N + geo.alt) * geo.n.z;
	return {x, y, z};
}

/** Convert {latitude, longitude, altitude} to {n-vector, altitude}.
 *
 * @param geo {latitude, longitude, altitude}
 * @return {n-vector, altitude}
 */
static inline nva lla2nva(lla geo) {
	double sin_lat = std::sin(geo.lat);
	double cos_lat = std::cos(geo.lat);
	double sin_lon = std::sin(geo.lon);
	double cos_lon = std::cos(geo.lon);
#ifndef NDEBUG
	// Check that altitude is within its defined range.
	double  N = wgs84::a / std::sqrt(1.0 - wgs84::e2*sin_lat*sin_lat);
	assert((1.0 - wgs84::e2) * N + geo.alt >= 0);
#endif
	return {{cos_lat * cos_lon, cos_lat * sin_lon, sin_lat}, geo.alt};
}

/** Convert normalized altitude to altitude within the valid range.
 *
 * This is useful when sampling coordinates since the valid altitude range is
 * dependent on the latitude. Hence, to cover the full altitude range, one
 * first have to sample the latitude and a normalized altitude after which the
 * normalized altitude can be converted to altitude.
 *
 * @param normAlt Normalized altitude in range [0,1].
 * @param maxAlt Maximum altitude.
 * @param lat The latitude.
 * @return An altitude within the valid range.
 */
static inline double normAlt2Alt(double normAlt, double maxAlt, double lat) {
	double sin_lat = std::sin(lat);
	double N = wgs84::a / std::sqrt(1.0 - wgs84::e2 * sin_lat * sin_lat);
	double minAlt = -(1.0 - wgs84::e2) * N;
	double rangeAlt = maxAlt - minAlt;
	return minAlt + normAlt * rangeAlt;
}

}  // namespace util_geo


#endif // TRANSFORMATION_TESTS_UTIL_GEO_HPP
