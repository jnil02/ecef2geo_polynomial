// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the BSD 2-Clause License.

/*
 * Polynomial ECEF coordinate (x,y,z) to geodetic coordinate
 * (latitude, longitude, altitude) or (n-vector, altitude) transformation
 * approximations.
 */

#ifndef ECEF2GEO_COMMON_HPP_
#define ECEF2GEO_COMMON_HPP_

// Input and output (coordinate) structs.
#include "coordinate_structs.hpp"

// Generated polynomial approximation for trigonometric functions.
#include "generated/xi.hpp"
#include "generated/chi.hpp"

#include <cstdint>  // int64_t, uint64_t
#define _USE_MATH_DEFINES
#include <cmath>  // std::copysign, std::sqrt, std::abs.

// Make sure pi constants are defined. These should normally be defined in
// cmath but in case they are not.
#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2		1.57079632679489661923
#endif
#ifndef M_PI_4
#define M_PI_4		0.78539816339744830962
#endif

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


namespace ecef2geo {

/*
 * Trigonometric function approximations used by the transformation
 * approximations. Placed in the private namespace to avoid bloating the outer
 * namespace.
 */
namespace priv {

/** Asin approximation.
 *
 * Has a discontinuous derivative at 0.
 *
 * @tparam N Polynomial degree of the approximation.
 * @param z argument within the range [0,1]
 * @return Arcsine approximation in the range [pi/2,-pi/2].
 */
template<int N>
inline double asin_approx(double z) {
	double x = std::abs(z);
	return std::copysign(M_PI_2 - std::sqrt(1.0 - x) * chi<N>(x), z);
}

/** Two argument atan minimax approximation.
 *
 * If the argument to atan can be split into a divident and a positive divisor,
 * an atan approximation for [-1,1] can be mapped to the range (-inf,inf)
 * without any conditional code or extra divide at the cost of two extra
 * additions. This is similar to atan2 but where the two arguments are used to
 * improve implementation performance.
 *
 * Given a minimax approximation of atan over [-1,1]
 *   atan(X) \approx \sum_{i=0}^{N}c_{2i+1}X^{2i+1} \forall x\in [-1,1]
 *           = f(x)
 * An equally good approximation of [0,inf) is
 *   atan(X) \approx f((X-1)/X+1)) \forall x\in [0,inf]
 * Now with X=|y|/x
 *   (X-1)/X+1) = (|y|/x-1)/(|y|/x+1)
 *              = ((|y|-x)/x)/((|y|+x)/x)
 *              = (|y|-x)/(|y|+x)
 * Finally, the output can be mapped to the full output range with the sign of
 * y, i.e.
 * atan(y/x) = copysign(f((|y| - x) / (|y| + x)), y) \forall x\in [0,inf]
 *
 * @tparam N Polynomial order of minimax polynomial.
 * @param y dividend of atan argument.
 * @param x __positive__ divisor of atan argument.
 * @return Approximation of atan(y/x)
 */
template<int N>
inline double atan_approx(double y, double x) {
	double abs_y = std::abs(y);
	double X = (abs_y - x) / (abs_y + x);
	return std::copysign(M_PI_4 + xi<N>(X), y);
}

/** Atan2 approximation from two argument atan approximation on [-inf,inf)
 *
 * Given an two argument atan approximation atan(y/x) = f(y,x) on [0,inf), the
 * approximation can be mapped the full angular range, i.e. atan2, by
 *
 *                      f(y,x) : y>0, x>0
 *   atan2(y,x) =  pi - f(y,x) : y>0, x<0
 *                -pi + f(y,x) : y<0, x<0
 *                    - f(y,x) : y<0, x>0
 *
 * All bits and pieces apart from f(y,x) can easily be implemented with bit
 * manipulations:
 *  - The sign of the f(y,x) is given by the xor of the signbits of x and y.
 *  - Further, the sign of the offset is the "signbit of x" & "signbit of y".
 *  - Finally, the magnitude of the offset is "signbit of x" times pi.
 *
 *  Note that y=0 and x=0 gives NaN.
 *
 * @tparam N Polynomial order of minimax polynomial.
 * @param y dividend of atan2 argument.
 * @param x divisor of atan2 argument.
 * @return Approximation of atan2(y,x)
 */
template<int N>
inline double atan2_approx(double y, double x) {
	// Just some constants for bit manipulation of double.
	constexpr int _63 = (sizeof(double) * 8 - 1);
	constexpr int _62 = _63 - 1;

	// Offset for the four quadrants.
	int64_t offset =
			((((*((int64_t *) &x)) & (*((int64_t *) &y))) >> _62) & -2) |
			((*((uint64_t *) &x)) >> _63);

	// XOR y and x to get a double with the sign of atan(z).
	uint64_t xXORy = (*((uint64_t *) &x)) ^ (*((uint64_t *) &y));

	double _xa = std::abs(x);
	double _ya = std::abs(y);
	double _x = (_ya - _xa) / (_ya + _xa);
	return M_PI * (double) offset +
		   std::copysign(M_PI_4 + xi<N>(_x), *((double *) &xXORy));;
}

}  // namespace priv
}  // namespace

#endif // ECEF2GEO_COMMON_HPP_
