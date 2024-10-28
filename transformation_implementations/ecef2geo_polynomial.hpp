// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the BSD 2-Clause License.

/*
 * Polynomial ECEF coordinate (x,y,z) to geodetic coordinate
 * (latitude, longitude, altitude) or (n-vector, altitude) transformation
 * approximations.
 */

#ifndef ECEF2GEO_POLYNOMIAL_HPP_
#define ECEF2GEO_POLYNOMIAL_HPP_

// Input and output (coordinate) structs.
#include "coordinate_structs.hpp"

// Generated polynomial approximations components.
#include "generated/mu.hpp"
#include "generated/omega.hpp"
#include "generated/sigma.hpp"
#include "generated/tau.hpp"

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

namespace ecef2geo {

// Assert that all polynomials are generated with the same ranges.
static_assert(OMEGA_H_MAX == MU_H_MAX, "Nonmatching H_max.");
static_assert(OMEGA_H_MAX == SIGMA_H_MAX, "Nonmatching H_max.");
static_assert(OMEGA_H_MAX == TAU_H_MAX, "Nonmatching H_max.");
static_assert(OMEGA_H_MIN == MU_H_MIN, "Nonmatching H_min.");
static_assert(OMEGA_H_MIN == SIGMA_H_MIN, "Nonmatching H_min.");
static_assert(OMEGA_H_MIN == TAU_H_MIN, "Nonmatching H_min.");
static_assert(SIGMA_DELTA_MIN == TAU_DELTA_MIN, "Nonmatching delta_max.");
static_assert(SIGMA_DELTA_MAX == TAU_DELTA_MAX, "Nonmatching delta_max.");
static_assert(OMEGA_H_0 == MU_H_0, "Nonmatching H_0.");

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


/** ECEF coordinate to geodetic coordinate (n-vector, altitude) transformation.
 *
 * The special case for which L=0 and N=0 is handled separately. This is just
 * for improved performance for this case and can be removed if it is not
 * relevant.
 *
 * @tparam L Polynomial order of multiplicative correction.
 * @tparam N Polynomial order of additive latitude correction for latitude.
 * @tparam M Polynomial order of additive altitude correction for latitude.
 * @tparam Nh Polynomial order of additive latitude correction for altitude.
 * @tparam Mh Polynomial order of additive altitude correction for altitude.
 * @param ecef The ECEF coordinate (x,y,z).
 * @return The geodetic coordinate (n-vector, altitude).
 */
template<int L, int N, int M, int Nh, int Mh>
static std::enable_if_t<L != 0 && N != 0, nva> ecef2nva(xyz ecef) {
	double &x = ecef.x;
	double &y = ecef.y;
	double &z = ecef.z;

	double p = std::sqrt(x * x + y * y + z * z);
	// This is to ensure that -0 is optimized away.
#ifndef H_0_IS_ZERO
	double h_c = p - OMEGA_H_0;
#else
	double &h_c = p;
#endif
	double t = z / p;
	double t2 = t * t;
	double alt = priv::mu<Nh, Mh>(h_c, t2);
	double w = priv::omega<N, M>(h_c, t2);
	double de = t2 * (1.0 - t2) * w * w;
	double si = priv::sigma<L>(de);
	double ta = priv::tau<L>(w, de);
	double inv_p_rho = (si - t2 * ta) / p;
	double t_eta = t * (si + (1.0 - t2) * ta);
	return {{x * inv_p_rho, y * inv_p_rho, t_eta}, alt};
}

/** Specialization to handle poor compiler performance.
 *
 * For some reason, some compilers appear to fail at generating efficient
 * code for the simplest cases. Potentially, it has to do with +0
 * analytically being a no-op but but for floating points (IEEE-754) it is
 * not.
 *
 * Special cases:
 *  - si == 1 for L in [0,1] or w = 0.
 *  - ta == 0 for L == 0 or w == 0.
 *
 *  w == 0 happens for N = 0.
 *  Is is only relevant to optimize for si == 1 when ta == 0.
 *
 * @tparam L Polynomial order of multiplicative correction.
 * @tparam N Polynomial order of additive latitude correction for latitude.
 * @tparam M Polynomial order of additive altitude correction for latitude.
 * @tparam Nh Polynomial order of additive latitude correction for altitude.
 * @tparam Mh Polynomial order of additive altitude correction for altitude.
 * @param ecef The ECEF coordinate (x,y,z).
 * @return The geodetic coordinate (n-vector, altitude).
 */
template<int L, int N, int M, int Nh, int Mh>
static std::enable_if_t<L == 0 || N == 0, nva> ecef2nva(xyz ecef) {
	double &x = ecef.x;
	double &y = ecef.y;
	double &z = ecef.z;

	double p = std::sqrt(x * x + y * y + z * z);
	// This is to ensure that -0 is optimized away.
#ifndef H_0_IS_ZERO
	double h_c = p - OMEGA_H_0;
#else
	double &h_c = p;
#endif
	double t = z / p;
	double t2 = t * t;
	double alt = priv::mu<Nh, Mh>(h_c, t2);
	return {{x / p, y / p, t}, alt};
}

/** ECEF to geodetic coordinate transformation by multiplicative arcsin correction.
 *
 * The special case for which L=0 and N=0 is handled separately. This is just
 * for improved performance for this case and can be removed if it is not
 * relevant.
 *
 * @tparam L Polynomial order of multiplicative correction.
 * @tparam N Polynomial order of additive latitude correction for latitude.
 * @tparam M Polynomial order of additive altitude correction for latitude.
 * @tparam Nh Polynomial order of additive latitude correction for altitude.
 * @tparam Mh Polynomial order of additive altitude correction for altitude.
 * @tparam I Polynomial order of arcsin approximation.
 * @tparam J Polynomial order of arctan2 approximation.
 * @param ecef The ECEF coordinate (x,y,z).
 * @return The geodetic coordinate (latitude, longitude, altitude).
 */
template<int L, int N, int M, int Nh, int Mh, int I, int J>
static std::enable_if_t<L != 0 && N != 0, lla> ecef2lla_asinm(xyz ecef) {
	double &x = ecef.x;
	double &y = ecef.y;
	double &z = ecef.z;

	double lon = priv::atan2_approx<J>(y, x);
	double p = std::sqrt(x * x + y * y + z * z);
	// This is to ensure that -0 is optimized away.
#ifndef H_0_IS_ZERO
	double h_c = p - OMEGA_H_0;
#else
	double &h_c = p;
#endif
	// Important to divide here rather than to divide by one and multiply.
	// For N>=4, the latter gives significant errors around the poles.
	double t = z / p;
	double t2 = t * t;
	double alt = priv::mu<Nh, Mh>(h_c, t2);
	double w = priv::omega<N, M>(h_c, t2);

	double de = t2 * (1.0 - t2) * w * w;
	double si = priv::sigma<L>(de);
	double ta = priv::tau<L>(w, de);
	double t_eta = t * (si + (1.0 - t2) * ta);

	double lat = priv::asin_approx<I>(t_eta);
	return {lat, lon, alt};
}

/** Specialization to handle poor compiler performance.
 *
 * For some reason, some compilers appear to fail at generating efficient
 * code for the simplest cases. Potentially, it has to do with +0
 * analytically being a no-op but but for floating points (IEEE-754) it is
 * not.
 *
 * Special cases:
 *  - si == 1 for L in [0,1] or w = 0.
 *  - ta == 0 for L == 0 or w == 0.
 *
 *  w == 0 happens for N = 0.
 *  Is is only relevant to optimize for si == 1 when ta == 0.
 *
 * @tparam L Polynomial order of multiplicative correction.
 * @tparam N Polynomial order of additive latitude correction for latitude.
 * @tparam M Polynomial order of additive altitude correction for latitude.
 * @tparam Nh Polynomial order of additive latitude correction for altitude.
 * @tparam Mh Polynomial order of additive altitude correction for altitude.
 * @tparam I Polynomial order of arcsin approximation.
 * @tparam J Polynomial order of arctan2 approximation.
 * @param ecef The ECEF coordinate (x,y,z).
 * @return The geodetic coordinate (latitude, longitude, altitude).
 */
template<int L, int N, int M, int Nh, int Mh, int I, int J>
static std::enable_if_t<L == 0 || N == 0, lla> ecef2lla_asinm(xyz ecef) {
	double &x = ecef.x;
	double &y = ecef.y;
	double &z = ecef.z;

	double lon = priv::atan2_approx<J>(y, x);
	double p = std::sqrt(x * x + y * y + z * z);
	// This is to ensure that -0 is optimized away.
#ifndef H_0_IS_ZERO
	double h_c = p - OMEGA_H_0;
#else
	double &h_c = p;
#endif
	// Important to divide here rather than to divide by one and multiply.
	// For N>=4, the latter gives significant errors around the poles.
	double t = z / p;
	double t2 = t * t;
	double alt = priv::mu<Nh, Mh>(h_c, t2);

	double lat = priv::asin_approx<I>(t);
	return {lat, lon, alt};
}

/** ECEF to geodetic coordinate transformation by multiplicative arctan correction.
 *
 * @tparam L Polynomial order of multiplicative correction.
 * @tparam N Polynomial order of additive latitude correction for latitude.
 * @tparam M Polynomial order of additive altitude correction for latitude.
 * @tparam Nh Polynomial order of additive latitude correction for altitude.
 * @tparam Mh Polynomial order of additive altitude correction for altitude.
 * @tparam I Polynomial order of arctan approximation.
 * @tparam J Polynomial order of arctan2 approximation.
 * @param ecef The ECEF coordinate (x,y,z).
 * @return The geodetic coordinate (latitude, longitude, altitude).
 */
template<int L, int N, int M, int Nh, int Mh, int I, int J>
static lla ecef2lla_atanm(xyz ecef) {
	double &x = ecef.x;
	double &y = ecef.y;
	double &z = ecef.z;

	double lon = priv::atan2_approx<J>(y, x);
	double p = std::sqrt(x * x + y * y + z * z);
	// This is to ensure that -0 is optimized away.
#ifndef H_0_IS_ZERO
	double h_c = p - OMEGA_H_0;
#else
	double &h_c = p;
#endif
	double t = z / p;
	double t2 = t * t;
	double alt = priv::mu<Nh, Mh>(h_c, t2);
	double w = priv::omega<N, M>(h_c, t2);
	double de = t2 * (1.0 - t2) * w * w;
	double si = priv::sigma<L>(de);
	double ta = priv::tau<L>(w, de);
	double rho = si - t2 * ta;
	double eta = si + (1.0 - t2) * ta;
	double lat = priv::atan_approx<I>(z * eta, rho * std::sqrt(x * x + y * y));
	return {lat, lon, alt};
}

/** ECEF to geodetic coordinate transformation by additive arcsin correction.
 *
 * @tparam N Polynomial order of additive latitude correction for latitude.
 * @tparam M Polynomial order of additive altitude correction for latitude.
 * @tparam Nh Polynomial order of additive latitude correction for altitude.
 * @tparam Mh Polynomial order of additive altitude correction for altitude.
 * @tparam I Polynomial order of arcsin approximation.
 * @tparam J Polynomial order of arctan2 approximation.
 * @param ecef The ECEF coordinate (x,y,z).
 * @return The geodetic coordinate (latitude, longitude, altitude).
 */
template<int N, int M, int Nh, int Mh, int I, int J>
static lla ecef2lla_asina(xyz ecef) {
	double &x = ecef.x;
	double &y = ecef.y;
	double &z = ecef.z;

	double lon = priv::atan2_approx<J>(y, x);
	double p = std::sqrt(x * x + y * y + z * z);
	// This is to ensure that -0 is optimized away.
#ifndef H_0_IS_ZERO
	double h_c = p - OMEGA_H_0;
#else
	double &h_c = p;
#endif
	double t = z / p;
	double t2 = t * t;
	double alt = priv::mu<Nh, Mh>(h_c, t2);
	double w = priv::omega<N, M>(h_c, t2);
	double lat = priv::asin_approx<I>(t) + t * std::sqrt(1.0 - t2) * w;
	return {lat, lon, alt};
}

/** ECEF to geodetic coordinate transformation by additive arctan correction.
 *
 * @tparam N Polynomial order of additive latitude correction for latitude.
 * @tparam M Polynomial order of additive altitude correction for latitude.
 * @tparam Nh Polynomial order of additive latitude correction for altitude.
 * @tparam Mh Polynomial order of additive altitude correction for altitude.
 * @tparam I Polynomial order of arcsin approximation.
 * @tparam J Polynomial order of arctan2 approximation.
 * @param ecef The ECEF coordinate (x,y,z).
 * @return The geodetic coordinate (latitude, longitude, altitude).
 */
template<int N, int M, int Nh, int Mh, int I, int J>
static lla ecef2lla_atana(xyz ecef) {
	double &x = ecef.x;
	double &y = ecef.y;
	double &z = ecef.z;

	double lon = priv::atan2_approx<J>(y, x);
	double p = std::sqrt(x * x + y * y + z * z);
	// This is to ensure that -0 is optimized away.
#ifndef H_0_IS_ZERO
	double h_c = p - OMEGA_H_0;
#else
	double &h_c = p;
#endif
	double t = z / p;
	double t2 = t * t;
	double alt = priv::mu<Nh, Mh>(h_c, t2);
	double w = priv::omega<N, M>(h_c, t2);
	double lat = priv::atan_approx<I>(z, std::sqrt(x * x + y * y)) +
				 t * std::sqrt(1.0 - t2) * w;
	return {lat, lon, alt};
}

}  // namespace

#endif // ECEF2GEO_POLYNOMIAL_HPP_
