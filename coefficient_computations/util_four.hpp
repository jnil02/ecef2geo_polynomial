// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the BSD 2-Clause License.

#ifndef ECEF2GEO_POLYNOMIAL_UTIL_FOUR_HPP
#define ECEF2GEO_POLYNOMIAL_UTIL_FOUR_HPP

#include "consts.hpp"

#include <mpreal.h>
#include <coefficients.hpp>


// Just a helper to convert numerator and denominator to an mpreal.
inline static mpfr::mpreal rc_expr(const point_to_ellipse_series::rc &d)
{
	return mpfr::mpreal(d.num) / mpfr::mpreal(d.den);
}

inline mpfr::mpreal d_phi(int n, int k) {
	mpfr::mpreal d(0);
	for (int l = std::max(k, n + 1); l <= k + n; ++l)
		d += rc_expr(point_to_ellipse_series::d_phi(n, k, l)) * mpfr::pow(mp_e2(), l);
	return d;
}

inline mpfr::mpreal d_sin(int n, int k) {
	mpfr::mpreal d(0);
	for (int l = std::max(k, n); l <= n + k; ++l)
		d += rc_expr(point_to_ellipse_series::d_sin(n, k, l)) * mpfr::pow(mp_e2(), l);
	return d;
}

inline mpfr::mpreal d_cos(int n, int k) {
	mpfr::mpreal d(0);
	for (int l = std::max(k, n); l < n + k; ++l)
		d += rc_expr(point_to_ellipse_series::d_cos(n, k, l)) * mpfr::pow(mp_e2(), l);
	return d;
}

inline mpfr::mpreal d_h(int n, int k) {
	mpfr::mpreal d(0);
	for (int l = std::max(k + 1, n); l <= n + k; ++l)
		d += rc_expr(point_to_ellipse_series::d_h(n, k, l)) * mpfr::pow(mp_e2(), l);
	return d;
}

#endif //ECEF2GEO_POLYNOMIAL_UTIL_FOUR_HPP
