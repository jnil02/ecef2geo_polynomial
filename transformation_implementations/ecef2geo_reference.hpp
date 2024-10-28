// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the BSD 2-Clause License.

#ifndef ECEF2GEO_REFERENCE_HPP
#define ECEF2GEO_REFERENCE_HPP

/*
 * Reference ECEF to geodetic coordinate transformations.
 */

#include "coordinate_structs.hpp"

#define _USE_MATH_DEFINES
#include <cmath>  // sqrt, sin, cos, atan, atan2, cbrt, abs

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

namespace vermeille {

using wgs84::a;
using wgs84::b;
using wgs84::e;
using wgs84::e2;

// Some constants for Vermeille's 2011 transformation.
constexpr double a2 = a * a;
constexpr double b2 = b * b;
constexpr double e4 = e2 * e2;
constexpr double b2am2 = b2 / a2;
constexpr double a3 = a2 * a;
constexpr double a4 = a2 * a2;
constexpr double am2 = 1. / a2;
constexpr double bam2 = b / a2;
constexpr double b2am4 = b2 / a4;
constexpr double inv_6 = 1. / 6.;
constexpr double e2_2 = e2 / 2;
constexpr double e2bam3 = e2 * b / a3;
constexpr double M_PI_6 = M_PI / 6.;
constexpr double M_2_3 = (2. / 3.);
constexpr double bam1 = b / a; // = 1 - f
constexpr double bem1 = b / e;
constexpr double e2bm1 = e2 / b;
constexpr double c = M_SQRT1_2 - 1.;

/** ECEF to geodetic coordinate transformation by Vermeille's 2011 method.
 *
 * Closed form ECEF to geodetic coordinate transformation as per:
 *
 * Vermeille, H, An analytical method to transform geocentric into geodetic
 * coorindates. J. Geod. (2011) 85:105-117. DOI 10.1007/s00190-010-0419-x
 *
 * Valid for all finite input ECEF coordinates.
 * ECEF coordinates on the singular disk are taken to be positive.
 *
 * @param ecef {x, y, z}
 * @return {latitude, longitude, altitude}
 */
static inline lla ecef2lla(xyz ecef) {
	double& x = ecef.x;
	double& y = ecef.y;
	double& z = ecef.z;

	double t2 = x * x + y * y;
	double t = std::sqrt(t2);
	double p = t2 * am2;
	double q = z * z * b2am4;
	double r = (p + q - e4) * inv_6;
	double r38 = r * r * r * 8.;
	double ev = r38 + e4 * p * q;

	// This is used rather than just a plain atan2 to ensure numerical stability
	// around the date line.
	double lon = c * std::abs(y) < t + x ?
				 2. * std::atan(y / (t + x)) : t + y < c * std::abs(x)
											   ? -M_PI_2 + 2. * std::atan(x / (t - y))
											   : M_PI_2 - 2. * std::atan(x / (t + y));

	if (ev > 0) {
		// Outside evolute.
		double s = e2bam3 * std::abs(z) * t;
		double sqrt_ev = std::sqrt(ev);
		double c1 = sqrt_ev + s;
		double c2 = sqrt_ev - s;
		double u = r + 0.5 * std::cbrt(c1 * c1) + 0.5 * std::cbrt(c2 * c2);

		double v = std::sqrt(u * u + e4 * q);
		double w = e2_2 * (u + v - q) / v;
		double k = (u + v) / (std::sqrt(w * w + u + v) + w);
		double D = k / (k + e2) * t;
		double d = std::sqrt(D * D + z * z);

		double h = (k - b2am2) * d / k;
		double lat = 2. * std::atan(z / (d + D));
		return {lat, lon, h};
	} else if (q != 0) {
		// On or inside evolute and not on singular disc.
		double s = e2bam3 * std::abs(z) * t;  // std::sqrt(e4*p*q)
		double up = M_2_3 * std::atan(s / (std::sqrt(-ev) + std::sqrt(-r38)));
		double u = -4. * r * std::sin(up) * std::cos(M_PI_6 + up);

		double v = std::sqrt(u * u + e4 * q);
		double w = e2_2 * (u + v - q) / v;
		double k = (u + v) / (std::sqrt(w * w + u + v) + w);
		double D = k * t / (k + e2);
		double d = std::sqrt(D * D + z * z);

		double h = (k - b2am2) * d / k;
		double lat = 2. * std::atan(z / (d + D));
		return {lat, lon, h};
	} else {
		// On the singular disc (including center of earth).
		// Values are taken to have positive latitude.
		double h = -bem1 * std::sqrt(e2 - p);
		double lat = 2. * std::atan(std::sqrt(e4 - p)
									/ (-e2bm1 * h + bam1 * std::sqrt(p)));
		return {lat, lon, h};
	}
}

/** ECEF to geodetic coordinate transformation by Vermeille
 *
 * Closed form ECEF to geodetic coordinate transformation as per:
 *
 * Vermeille, H, An analytical method to transform geocentric into geodetic
 * coorindates. J. Geod. (2011) 85:105-117. DOI 10.1007/s00190-010-0419-x
 *
 * Valid for all finite input ECEF coordinates.
 * ECEF coordinates on the singular disk are taken to be positive.
 *
 * @param ecef {x, y, z}
 * @return {n-vector, altitude}
 */
static inline nva ecef2nva(xyz ecef) {
	double& x = ecef.x;
	double& y = ecef.y;
	double& z = ecef.z;

	double t2 = x * x + y * y;
	double t = std::sqrt(t2);
	double p = t2 * am2;
	double q = z * z * b2am4;
	double r = (p + q - e4) * inv_6;
	double r38 = r * r * r * 8.;
	double ev = r38 + e4 * p * q;

	if (ev > 0) {
		// Outside evolute.
		double s = e2bam3 * std::abs(z) * t;
		double sqrt_evo = std::sqrt(ev);
		double c1 = sqrt_evo + s;
		double c2 = sqrt_evo - s;
		double u = r + 0.5 * std::cbrt(c1 * c1) + 0.5 * std::cbrt(c2 * c2);

		double v = std::sqrt(u * u + e4 * q);
		double w = e2_2 * (u + v - q) / v;
		double k = (u + v) / (std::sqrt(w * w + u + v) + w);
		double l = k / (k + e2);
		double D = l * t;
		double d = std::sqrt(D * D + z * z);

		double h = (k - b2am2) * d / k;
		double n = 1. / d;
		return {{n * l * x, n * l * y, n * z}, h};
	} else if (q != 0) {
		// On or inside evolute and not on singular disc.
		double s = e2bam3 * std::abs(z) * t;
		double up = M_2_3 * std::atan2(s, std::sqrt(-ev) + std::sqrt(-r38));
		double u = -4. * r * std::sin(up) * std::cos(M_PI_6 + up);

		double v = std::sqrt(u * u + e4 * q);
		double w = e2_2 * (u + v - q) / v;
		double k = (u + v) / (std::sqrt(w * w + u + v) + w);
		double l = k / (k + e2);
		double D = l * t;
		double d = std::sqrt(D * D + z * z);

		double h = (k - b2am2) * d / k;
		double n = 1. / d;
		return {{x * l * n, y * l * n, z * n}, h};
	} else {
		// On the singular disc (including center of earth).
		// Values are taken to have positive latitude.
		double h = -bem1 * std::sqrt(e2 - p);
		double n = 1. / (-e2bm1 * h);
		return {{x * bam2 * n, y * bam2 * n, std::sqrt(e4 - p) * n}, h};
	}
}
} // namespace vermeille

#endif // ECEF2GEO_REFERENCE_HPP
