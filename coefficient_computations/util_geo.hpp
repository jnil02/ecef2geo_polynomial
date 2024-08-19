// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the GPLv3 License.

#ifndef ECEF2GEO_POLYNOMIAL_UTIL_GEO_HPP
#define ECEF2GEO_POLYNOMIAL_UTIL_GEO_HPP

/*
 * Utility functions for multiprecision (mpfr) geodetic calculations.
 */

#include <mpfr.h>
#include "mplapack/mpreal.h"

#include <functional>  // std::functional.

using namespace mpfr;

// Latitude and altitude.
struct la_mpfr {
	mpfr_t lat;
	mpfr_t alt;
};

// "xyz" to contain an Earth-Centered Earth Fixed (ECEF) coordinate.
struct xyz_mpfr {
	mpfr_t x;
	mpfr_t y;
	mpfr_t z;
};

struct wgs84_mpfr_constants {
	// Primary WGS84 constants.
	mpfr_t a_mpfr;  // Semi-major axis.
	mpfr_t f_inv_mpfr;  // Inverse flattening.

	// Derived WGS84 constants.
	mpfr_t a2_mpfr;  // Semi-major axis squared.
	mpfr_t a2_inv_mpfr; // Inverse semi-major axis squared.
	mpfr_t f_mpfr;  // Flattening.
	mpfr_t b_mpfr;  // Semi-minor axis / polar radius.
	mpfr_t b2_mpfr;  // Semi-minor axis / polar radius squared.
	mpfr_t e2_mpfr;  // First eccentricity squared.
	mpfr_t e4_mpfr;  // First eccentricity to the power of 4.
	mpfr_t ec2_mpfr;  // Axis ratio squared.
	mpfr_t ec2_a2_mpfr;

	// Initialization of WGS84 constant.
	wgs84_mpfr_constants() {
		// These literal values are from the WGS84 standard and are exact.
		// To use GRS 80, instead use "298.257222100882711243" for the inverse
		// flattening. (Or better, compete it will full precision from the
		// underlying constants of GRS 80.)
		mpfr_init_set_str(a_mpfr, "6378137.0", 10, MPFR_RNDN);
		mpfr_init_set_str(f_inv_mpfr, "298.257223563", 10, MPFR_RNDN);

		mpfr_init(a2_mpfr);
		mpfr_mul(a2_mpfr, a_mpfr, a_mpfr, MPFR_RNDN);
		mpfr_init(a2_inv_mpfr);
		mpfr_ui_div(a2_inv_mpfr, 1, a2_mpfr, MPFR_RNDN);

		mpfr_init(f_mpfr);
		mpfr_ui_div(f_mpfr, 1, f_inv_mpfr, MPFR_RNDN);

		mpfr_init(b_mpfr);
		mpfr_mul(b_mpfr, f_mpfr, a_mpfr, MPFR_RNDN);
		mpfr_sub(b_mpfr, a_mpfr, b_mpfr, MPFR_RNDN);
		mpfr_init(b2_mpfr);
		mpfr_mul(b2_mpfr, b_mpfr, b_mpfr, MPFR_RNDN);

		// 1.0 - (b * b) / (a * a)
		mpfr_init(e2_mpfr);
		mpfr_div(e2_mpfr, b2_mpfr, a2_mpfr, MPFR_RNDN);
		mpfr_ui_sub(e2_mpfr, 1, e2_mpfr, MPFR_RNDN);
		mpfr_init(e4_mpfr);
		mpfr_mul(e4_mpfr, e2_mpfr, e2_mpfr, MPFR_RNDN);

		mpfr_init(ec2_mpfr);
		mpfr_si_sub(ec2_mpfr, 1, e2_mpfr, MPFR_RNDN);
		mpfr_init(ec2_a2_mpfr);
		mpfr_div(ec2_a2_mpfr, ec2_mpfr, a2_mpfr, MPFR_RNDN);
	}
};

// Struct and not namespace to enable having the static mpfr-stuff.
struct util_geo {

	// Statically allocated mpfr constants.
	static wgs84_mpfr_constants _static_stuff;

	/** Computes latitude and altitude from an ECEF coordinate.
	 *
	 * The computations are done according to Vermielle's 2004 exact method.
	 *
	 * @param la latitude-altitude output struct.
	 * @param ecef Ecef coordinate.
	 */
	static void ecef_to_lla(la_mpfr &la, const xyz_mpfr &ecef) {
		mp_prec_t prec = mpfr_get_prec(la.lat);
		mpfr_t x2y2, sqrtx2y2, z2, tmp;
		mpfr_init2(tmp, prec);
		mpfr_init2(x2y2, prec);
		mpfr_init2(sqrtx2y2, prec);
		mpfr_init2(z2, prec);

		mpfr_sqr(tmp, ecef.x, MPFR_RNDN);
		mpfr_sqr(x2y2, ecef.y, MPFR_RNDN);
		mpfr_add(x2y2, x2y2, tmp, MPFR_RNDN);
		mpfr_sqrt(sqrtx2y2, x2y2, MPFR_RNDN);
		mpfr_sqr(z2, ecef.z, MPFR_RNDN);
//        double x2y2 = ecef.x * ecef.x + ecef.y * ecef.y;
//        double sqrtx2y2 = std::sqrt(x2y2);
//        double z2 = ecef.z * ecef.z;

		mpfr_t q, p, r, s, t, u, v, w, k, k_k_e2, d, tmp2;
		mpfr_init2(q, prec);
		mpfr_init2(p, prec);
		mpfr_init2(r, prec);
		mpfr_init2(s, prec);
		mpfr_init2(t, prec);
		mpfr_init2(u, prec);
		mpfr_init2(v, prec);
		mpfr_init2(w, prec);
		mpfr_init2(k, prec);
		mpfr_init2(k_k_e2, prec);
		mpfr_init2(d, prec);
		mpfr_init2(tmp2, prec);
		mpfr_mul(q, z2, _static_stuff.ec2_a2_mpfr, MPFR_RNDN);
//        double q = ec2_a2 * z2;
		mpfr_mul(p, x2y2, _static_stuff.a2_inv_mpfr, MPFR_RNDN);
//        double p = x2y2 * a2_inv;
		mpfr_add(r, p, q, MPFR_RNDN);
		mpfr_sub(r, r, _static_stuff.e4_mpfr, MPFR_RNDN);
		mpfr_div_ui(r, r, 6, MPFR_RNDN);
//        double r = (p + q - e4) / 6.0;
		mpfr_mul(s, _static_stuff.e4_mpfr, p, MPFR_RNDN);
		mpfr_mul(s, s, q, MPFR_RNDN);
		mpfr_sqr(tmp, r, MPFR_RNDN);
		mpfr_mul(tmp, tmp, r, MPFR_RNDN);
		mpfr_mul_ui(tmp, tmp, 4, MPFR_RNDN);
		mpfr_div(s, s, tmp, MPFR_RNDN);
//        double s = e4 * p * q / (4.0 * r * r * r);
		mpfr_add_ui(t, s, 2, MPFR_RNDN);
		mpfr_mul(t, s, t, MPFR_RNDN);
		mpfr_sqrt(t, t, MPFR_RNDN);
		mpfr_add(t, t, s, MPFR_RNDN);
		mpfr_add_ui(t, t, 1, MPFR_RNDN);
		mpfr_cbrt(t, t, MPFR_RNDN);
//        double t = std::cbrt(1.0 + s + std::sqrt(s * (2.0 + s)));
		mpfr_ui_div(u, 1, t, MPFR_RNDN);
		mpfr_add(u, u, t, MPFR_RNDN);
		mpfr_add_ui(u, u, 1, MPFR_RNDN);
		mpfr_mul(u, u, r, MPFR_RNDN);
//        double u = r * (1.0 + t + 1.0 / t);
		mpfr_sqr(v, u, MPFR_RNDN);
		mpfr_mul(tmp, q, _static_stuff.e4_mpfr, MPFR_RNDN);
		mpfr_add(v, v, tmp, MPFR_RNDN);
		mpfr_sqrt(v, v, MPFR_RNDN);
//        double v = std::sqrt(u * u + e4 * q);
		mpfr_mul_ui(tmp, v, 2, MPFR_RNDN);
		mpfr_add(w, u, v, MPFR_RNDN);
		mpfr_sub(w, w, q, MPFR_RNDN);
		mpfr_div(w, w, tmp, MPFR_RNDN);
		mpfr_mul(w, w, _static_stuff.e2_mpfr, MPFR_RNDN);
//        double w = e2 * (u + v - q) / (2.0 * v);
		mpfr_sqr(k, w, MPFR_RNDN);
		mpfr_add(k, k, v, MPFR_RNDN);
		mpfr_add(k, k, u, MPFR_RNDN);
		mpfr_sqrt(k, k, MPFR_RNDN);
		mpfr_sub(k, k, w, MPFR_RNDN);
//        double k = std::sqrt(u + v + w * w) - w;
		mpfr_add(k_k_e2, k, _static_stuff.e2_mpfr, MPFR_RNDN);
		mpfr_div(k_k_e2, k, k_k_e2, MPFR_RNDN);
//        double k_k_e2 = k / (k + e2);
		mpfr_mul(d, k_k_e2, sqrtx2y2, MPFR_RNDN);
//        double d = k_k_e2 * sqrtx2y2;

		mpfr_sqr(tmp2, d, MPFR_RNDN);
		mpfr_add(tmp2, tmp2, z2, MPFR_RNDN);
		mpfr_sqrt(tmp2, tmp2, MPFR_RNDN);
//        double tmp = std::sqrt(d * d + z2);
		mpfr_add(la.lat, d, tmp2, MPFR_RNDN);
		mpfr_div(la.lat, ecef.z, la.lat, MPFR_RNDN);
		mpfr_atan(la.lat, la.lat, MPFR_RNDN);
		mpfr_mul_ui(la.lat, la.lat, 2, MPFR_RNDN);
//        double lat = 2.0 * std::atan(ecef.z / (d + tmp));
//        double lon = ecef.y >= 0.0 ?
//                     M_PI_2 - 2.0 * std::atan(ecef.x / (sqrtx2y2 + ecef.y)) :
//                     -M_PI_2 + 2.0 * std::atan(ecef.x / (sqrtx2y2 - ecef.y));
		mpfr_add(la.alt, k, _static_stuff.e2_mpfr, MPFR_RNDN);
		mpfr_sub_ui(la.alt, la.alt, 1, MPFR_RNDN);
		mpfr_div(la.alt, la.alt, k, MPFR_RNDN);
		mpfr_mul(la.alt, la.alt, tmp2, MPFR_RNDN);
//        double alt = (k + e2 - 1.0) / k * tmp;
//        return {lat, lon, alt};
		mpfr_clear(x2y2);
		mpfr_clear(sqrtx2y2);
		mpfr_clear(z2);
		mpfr_clear(tmp);
		mpfr_clear(q);
		mpfr_clear(p);
		mpfr_clear(r);
		mpfr_clear(s);
		mpfr_clear(t);
		mpfr_clear(u);
		mpfr_clear(v);
		mpfr_clear(w);
		mpfr_clear(k);
		mpfr_clear(k_k_e2);
		mpfr_clear(d);
		mpfr_clear(tmp2);
	}

	/** Computes latitude from an ECEF coordinate.
	 *
	 * The computations are done according to Vermielle's 2004 exact method.
	 *
	 * @param la latitude(-altitude) output struct.
	 * @param ecef ECEF coordinate.
	 */
	static void ecef_to_lat(la_mpfr &la, const xyz_mpfr &ecef) {
		mp_prec_t prec = mpfr_get_prec(la.lat);
		mpfr_t x2y2, sqrtx2y2, z2, tmp;
		mpfr_init2(tmp, prec);
		mpfr_init2(x2y2, prec);
		mpfr_init2(sqrtx2y2, prec);
		mpfr_init2(z2, prec);

		mpfr_sqr(tmp, ecef.x, MPFR_RNDN);
		mpfr_sqr(x2y2, ecef.y, MPFR_RNDN);
		mpfr_add(x2y2, x2y2, tmp, MPFR_RNDN);
		mpfr_sqrt(sqrtx2y2, x2y2, MPFR_RNDN);
		mpfr_sqr(z2, ecef.z, MPFR_RNDN);
//        double x2y2 = ecef.x * ecef.x + ecef.y * ecef.y;
//        double sqrtx2y2 = std::sqrt(x2y2);
//        double z2 = ecef.z * ecef.z;

		mpfr_t q, p, r, s, t, u, v, w, k, k_k_e2, d, tmp2;
		mpfr_init2(q, prec);
		mpfr_init2(p, prec);
		mpfr_init2(r, prec);
		mpfr_init2(s, prec);
		mpfr_init2(t, prec);
		mpfr_init2(u, prec);
		mpfr_init2(v, prec);
		mpfr_init2(w, prec);
		mpfr_init2(k, prec);
		mpfr_init2(k_k_e2, prec);
		mpfr_init2(d, prec);
		mpfr_init2(tmp2, prec);
		mpfr_mul(q, z2, _static_stuff.ec2_a2_mpfr, MPFR_RNDN);
//        double q = ec2_a2 * z2;
		mpfr_mul(p, x2y2, _static_stuff.a2_inv_mpfr, MPFR_RNDN);
//        double p = x2y2 * a2_inv;
		mpfr_add(r, p, q, MPFR_RNDN);
		mpfr_sub(r, r, _static_stuff.e4_mpfr, MPFR_RNDN);
		mpfr_div_ui(r, r, 6, MPFR_RNDN);
//        double r = (p + q - e4) / 6.0;
		mpfr_mul(s, _static_stuff.e4_mpfr, p, MPFR_RNDN);
		mpfr_mul(s, s, q, MPFR_RNDN);
		mpfr_sqr(tmp, r, MPFR_RNDN);
		mpfr_mul(tmp, tmp, r, MPFR_RNDN);
		mpfr_mul_ui(tmp, tmp, 4, MPFR_RNDN);
		mpfr_div(s, s, tmp, MPFR_RNDN);
//        double s = e4 * p * q / (4.0 * r * r * r);
		mpfr_add_ui(t, s, 2, MPFR_RNDN);
		mpfr_mul(t, s, t, MPFR_RNDN);
		mpfr_sqrt(t, t, MPFR_RNDN);
		mpfr_add(t, t, s, MPFR_RNDN);
		mpfr_add_ui(t, t, 1, MPFR_RNDN);
		mpfr_cbrt(t, t, MPFR_RNDN);
//        double t = std::cbrt(1.0 + s + std::sqrt(s * (2.0 + s)));
		mpfr_ui_div(u, 1, t, MPFR_RNDN);
		mpfr_add(u, u, t, MPFR_RNDN);
		mpfr_add_ui(u, u, 1, MPFR_RNDN);
		mpfr_mul(u, u, r, MPFR_RNDN);
//        double u = r * (1.0 + t + 1.0 / t);
		mpfr_sqr(v, u, MPFR_RNDN);
		mpfr_mul(tmp, q, _static_stuff.e4_mpfr, MPFR_RNDN);
		mpfr_add(v, v, tmp, MPFR_RNDN);
		mpfr_sqrt(v, v, MPFR_RNDN);
//        double v = std::sqrt(u * u + e4 * q);
		mpfr_mul_ui(tmp, v, 2, MPFR_RNDN);
		mpfr_add(w, u, v, MPFR_RNDN);
		mpfr_sub(w, w, q, MPFR_RNDN);
		mpfr_div(w, w, tmp, MPFR_RNDN);
		mpfr_mul(w, w, _static_stuff.e2_mpfr, MPFR_RNDN);
//        double w = e2 * (u + v - q) / (2.0 * v);
		mpfr_sqr(k, w, MPFR_RNDN);
		mpfr_add(k, k, v, MPFR_RNDN);
		mpfr_add(k, k, u, MPFR_RNDN);
		mpfr_sqrt(k, k, MPFR_RNDN);
		mpfr_sub(k, k, w, MPFR_RNDN);
//        double k = std::sqrt(u + v + w * w) - w;
		mpfr_add(k_k_e2, k, _static_stuff.e2_mpfr, MPFR_RNDN);
		mpfr_div(k_k_e2, k, k_k_e2, MPFR_RNDN);
//        double k_k_e2 = k / (k + e2);
		mpfr_mul(d, k_k_e2, sqrtx2y2, MPFR_RNDN);
//        double d = k_k_e2 * sqrtx2y2;

		mpfr_sqr(tmp2, d, MPFR_RNDN);
		mpfr_add(tmp2, tmp2, z2, MPFR_RNDN);
		mpfr_sqrt(tmp2, tmp2, MPFR_RNDN);
//        double tmp = std::sqrt(d * d + z2);
		mpfr_add(la.lat, d, tmp2, MPFR_RNDN);
		mpfr_div(la.lat, ecef.z, la.lat, MPFR_RNDN);
		mpfr_atan(la.lat, la.lat, MPFR_RNDN);
		mpfr_mul_ui(la.lat, la.lat, 2, MPFR_RNDN);
//        double lat = 2.0 * std::atan(ecef.z / (d + tmp));
//        double lon = ecef.y >= 0.0 ?
//                     M_PI_2 - 2.0 * std::atan(ecef.x / (sqrtx2y2 + ecef.y)) :
//                     -M_PI_2 + 2.0 * std::atan(ecef.x / (sqrtx2y2 - ecef.y));
//        mpfr_add(la.alt, k, _static_stuff.e2_mpfr, MPFR_RNDN);
//        mpfr_sub_ui(la.alt, la.alt, 1, MPFR_RNDN);
//        mpfr_div(la.alt, la.alt, k, MPFR_RNDN);
//        mpfr_mul(la.alt, la.alt, tmp2, MPFR_RNDN);
//        double alt = (k + e2 - 1.0) / k * tmp;
//        return {lat, lon, alt};
		mpfr_clear(x2y2);
		mpfr_clear(sqrtx2y2);
		mpfr_clear(z2);
		mpfr_clear(tmp);
		mpfr_clear(q);
		mpfr_clear(p);
		mpfr_clear(r);
		mpfr_clear(s);
		mpfr_clear(t);
		mpfr_clear(u);
		mpfr_clear(v);
		mpfr_clear(w);
		mpfr_clear(k);
		mpfr_clear(k_k_e2);
		mpfr_clear(d);
		mpfr_clear(tmp2);
	}

	/** Computes altitude from an ECEF coordinate.
	 *
	 * The computations are done according to Vermielle's 2004 exact method.
	 *
	 * @param la (latitude-)altitude output struct.
	 * @param ecef ECEF coordinate.
	 */
	static void ecef_to_alt(la_mpfr &la, const xyz_mpfr &ecef) {
		mp_prec_t prec = mpfr_get_prec(la.lat);
		mpfr_t x2y2, sqrtx2y2, z2, tmp;
		mpfr_init2(tmp, prec);
		mpfr_init2(x2y2, prec);
		mpfr_init2(sqrtx2y2, prec);
		mpfr_init2(z2, prec);

		mpfr_sqr(tmp, ecef.x, MPFR_RNDN);
		mpfr_sqr(x2y2, ecef.y, MPFR_RNDN);
		mpfr_add(x2y2, x2y2, tmp, MPFR_RNDN);
		mpfr_sqrt(sqrtx2y2, x2y2, MPFR_RNDN);
		mpfr_sqr(z2, ecef.z, MPFR_RNDN);
//        double x2y2 = ecef.x * ecef.x + ecef.y * ecef.y;
//        double sqrtx2y2 = std::sqrt(x2y2);
//        double z2 = ecef.z * ecef.z;

		mpfr_t q, p, r, s, t, u, v, w, k, k_k_e2, d, tmp2;
		mpfr_init2(q, prec);
		mpfr_init2(p, prec);
		mpfr_init2(r, prec);
		mpfr_init2(s, prec);
		mpfr_init2(t, prec);
		mpfr_init2(u, prec);
		mpfr_init2(v, prec);
		mpfr_init2(w, prec);
		mpfr_init2(k, prec);
		mpfr_init2(k_k_e2, prec);
		mpfr_init2(d, prec);
		mpfr_init2(tmp2, prec);
		mpfr_mul(q, z2, _static_stuff.ec2_a2_mpfr, MPFR_RNDN);
//        double q = ec2_a2 * z2;
		mpfr_mul(p, x2y2, _static_stuff.a2_inv_mpfr, MPFR_RNDN);
//        double p = x2y2 * a2_inv;
		mpfr_add(r, p, q, MPFR_RNDN);
		mpfr_sub(r, r, _static_stuff.e4_mpfr, MPFR_RNDN);
		mpfr_div_ui(r, r, 6, MPFR_RNDN);
//        double r = (p + q - e4) / 6.0;
		mpfr_mul(s, _static_stuff.e4_mpfr, p, MPFR_RNDN);
		mpfr_mul(s, s, q, MPFR_RNDN);
		mpfr_sqr(tmp, r, MPFR_RNDN);
		mpfr_mul(tmp, tmp, r, MPFR_RNDN);
		mpfr_mul_ui(tmp, tmp, 4, MPFR_RNDN);
		mpfr_div(s, s, tmp, MPFR_RNDN);
//        double s = e4 * p * q / (4.0 * r * r * r);
		mpfr_add_ui(t, s, 2, MPFR_RNDN);
		mpfr_mul(t, s, t, MPFR_RNDN);
		mpfr_sqrt(t, t, MPFR_RNDN);
		mpfr_add(t, t, s, MPFR_RNDN);
		mpfr_add_ui(t, t, 1, MPFR_RNDN);
		mpfr_cbrt(t, t, MPFR_RNDN);
//        double t = std::cbrt(1.0 + s + std::sqrt(s * (2.0 + s)));
		mpfr_ui_div(u, 1, t, MPFR_RNDN);
		mpfr_add(u, u, t, MPFR_RNDN);
		mpfr_add_ui(u, u, 1, MPFR_RNDN);
		mpfr_mul(u, u, r, MPFR_RNDN);
//        double u = r * (1.0 + t + 1.0 / t);
		mpfr_sqr(v, u, MPFR_RNDN);
		mpfr_mul(tmp, q, _static_stuff.e4_mpfr, MPFR_RNDN);
		mpfr_add(v, v, tmp, MPFR_RNDN);
		mpfr_sqrt(v, v, MPFR_RNDN);
//        double v = std::sqrt(u * u + e4 * q);
		mpfr_mul_ui(tmp, v, 2, MPFR_RNDN);
		mpfr_add(w, u, v, MPFR_RNDN);
		mpfr_sub(w, w, q, MPFR_RNDN);
		mpfr_div(w, w, tmp, MPFR_RNDN);
		mpfr_mul(w, w, _static_stuff.e2_mpfr, MPFR_RNDN);
//        double w = e2 * (u + v - q) / (2.0 * v);
		mpfr_sqr(k, w, MPFR_RNDN);
		mpfr_add(k, k, v, MPFR_RNDN);
		mpfr_add(k, k, u, MPFR_RNDN);
		mpfr_sqrt(k, k, MPFR_RNDN);
		mpfr_sub(k, k, w, MPFR_RNDN);
//        double k = std::sqrt(u + v + w * w) - w;
		mpfr_add(k_k_e2, k, _static_stuff.e2_mpfr, MPFR_RNDN);
		mpfr_div(k_k_e2, k, k_k_e2, MPFR_RNDN);
//        double k_k_e2 = k / (k + e2);
		mpfr_mul(d, k_k_e2, sqrtx2y2, MPFR_RNDN);
//        double d = k_k_e2 * sqrtx2y2;

		mpfr_sqr(tmp2, d, MPFR_RNDN);
		mpfr_add(tmp2, tmp2, z2, MPFR_RNDN);
		mpfr_sqrt(tmp2, tmp2, MPFR_RNDN);
//        double tmp = std::sqrt(d * d + z2);
//        mpfr_add(la.lat, d, tmp2, MPFR_RNDN);
//        mpfr_div(la.lat, ecef.z, la.lat, MPFR_RNDN);
//        mpfr_atan(la.lat, la.lat, MPFR_RNDN);
//        mpfr_mul_ui(la.lat, la.lat, 2, MPFR_RNDN);
//        double lat = 2.0 * std::atan(ecef.z / (d + tmp));
//        double lon = ecef.y >= 0.0 ?
//                     M_PI_2 - 2.0 * std::atan(ecef.x / (sqrtx2y2 + ecef.y)) :
//                     -M_PI_2 + 2.0 * std::atan(ecef.x / (sqrtx2y2 - ecef.y));
		mpfr_add(la.alt, k, _static_stuff.e2_mpfr, MPFR_RNDN);
		mpfr_sub_ui(la.alt, la.alt, 1, MPFR_RNDN);
		mpfr_div(la.alt, la.alt, k, MPFR_RNDN);
		mpfr_mul(la.alt, la.alt, tmp2, MPFR_RNDN);
//        double alt = (k + e2 - 1.0) / k * tmp;
//        return {lat, lon, alt};
		mpfr_clear(x2y2);
		mpfr_clear(sqrtx2y2);
		mpfr_clear(z2);
		mpfr_clear(tmp);
		mpfr_clear(q);
		mpfr_clear(p);
		mpfr_clear(r);
		mpfr_clear(s);
		mpfr_clear(t);
		mpfr_clear(u);
		mpfr_clear(v);
		mpfr_clear(w);
		mpfr_clear(k);
		mpfr_clear(k_k_e2);
		mpfr_clear(d);
		mpfr_clear(tmp2);
	}

	/** Computes latitude from spherical latitude and altitude.
	 *
	 * The computations are done according to Vermielle's 2004 method.
	 *
	 * @param lat_c spherical latitude.
	 * @param h_c spherical altitude.
	 * @param h_0 Reference for the spherical altitude.
	 * @param prec Precision of the arithmetics.
	 * @return the latitude.
	 */
	static mpreal f_c(const mpreal &lat_c, const mpreal &h_c, const mpreal &h_0,
					  mp_prec_t prec) {
		mpreal p;
		p.set_prec(prec);
		p = h_c + h_0;
		mpreal r;
		r.set_prec(prec);
		r = p * cos(lat_c);
		mpreal z;
		z.set_prec(prec);
		z = p * sin(lat_c);

		la_mpfr la;
		mpfr_init2(la.lat, prec);
		mpfr_init2(la.alt, prec);
		xyz_mpfr ecef;
		mpfr_init2(ecef.x, prec);
		mpfr_init2(ecef.y, prec);
		mpfr_init2(ecef.z, prec);

		mpfr_ptr r_ptr(r);
		mpfr_set(ecef.x, r_ptr, MPFR_RNDN);
		mpfr_set_d(ecef.y, 0.0, MPFR_RNDN);
		mpfr_ptr z_ptr(z);
		mpfr_set(ecef.z, z_ptr, MPFR_RNDN);

		util_geo::ecef_to_lat(la, ecef);
		mpreal res(la.lat);
		mpfr_clear(ecef.x);
		mpfr_clear(ecef.y);
		mpfr_clear(ecef.z);
		mpfr_clear(la.lat);
		mpfr_clear(la.alt);
		return res;
	}

	/** Computes altitude from spherical latitude and altitude.
	 *
	 * The computations are done according to Vermielle's 2004 method.
	 *
	 * @param lat_c spherical latitude.
	 * @param h_c spherical altitude.
	 * @param h_0 Reference for the spherical altitude.
	 * @param prec Precision of the arithmetics.
	 * @return the altitude.
	 */
	static mpreal g_c(const mpreal &lat_c, const mpreal &h_c, const mpreal &h_0,
					  mp_prec_t prec) {
		mpreal p;
		p.set_prec(prec);
		p = h_c + h_0;
		mpreal r;
		r.set_prec(prec);
		r = p * cos(lat_c);
		mpreal z;
		z.set_prec(prec);
		z = p * sin(lat_c);

		la_mpfr la;
		mpfr_init2(la.lat, prec);
		mpfr_init2(la.alt, prec);
		xyz_mpfr ecef;
		mpfr_init2(ecef.x, prec);
		mpfr_init2(ecef.y, prec);
		mpfr_init2(ecef.z, prec);

		mpfr_ptr r_ptr(r);
		mpfr_set(ecef.x, r_ptr, MPFR_RNDN);
		mpfr_set_d(ecef.y, 0.0, MPFR_RNDN);
		mpfr_ptr z_ptr(z);
		mpfr_set(ecef.z, z_ptr, MPFR_RNDN);

		util_geo::ecef_to_alt(la, ecef);
		mpreal res(la.alt);
		mpfr_clear(ecef.x);
		mpfr_clear(ecef.y);
		mpfr_clear(ecef.z);
		mpfr_clear(la.lat);
		mpfr_clear(la.alt);
		return res;
	}

	/** Integration by Simpson's 1/3 rule.
	 *
	 * @param f Function to integrate.
	 * @param a Lower integration limit.
	 * @param b Upper integration limit.
	 * @param N Number of integration sub-intervals.
	 * @param prec Precision of arithmetics.
	 * @return the integral value.
	 */
	static mpreal
	simpson_integration(const std::function<mpreal(const mpreal &)> &f,
						const mpreal &a, const mpreal &b, int N,
						mp_prec_t prec) {
		// N+1 evenly spaced latitude samples including both end-points.
		mpreal *lats = new mpreal[N + 1];
		mpreal step;
		step.set_prec(prec);
		step = (b - a) / N;
		for (int i = 0; i < N; ++i) {
			lats[i].set_prec(prec);
			lats[i] = step * i + a;
		}
		lats[N].set_prec(prec);
		lats[N] = b;  // To ensure that last point is exactly "b".

		// Weights for composite Simpson's 1/3 rule.
		mpreal h3;
		h3.set_prec(prec);
		h3 = (b - a) / (3 * N);  // h/3
		mpreal *ws = new mpreal[N + 1];
		ws[0].set_prec(prec);
		ws[0] = h3;
		for (int i = 1; i < N; ++i) {
			ws[i].set_prec(prec);
			if (i % 2 == 0)
				ws[i] = h3 * 2;
			else
				ws[i] = h3 * 4;
		}
		ws[N].set_prec(prec);
		ws[N] = h3;

		// Compute the integral.
		mpreal *vals = new mpreal[N + 1];
		for (int i = 0; i < N + 1; ++i) {
			vals[i].set_prec(prec);
			vals[i] = ws[i] * f(lats[i]);
		}
		mpreal integral;
		integral.set_prec(prec);
		integral = sum(vals, N + 1, MPFR_RNDN);

		delete[] lats;
		delete[] ws;
		delete[] vals;

		return integral;
	}
};

#endif // ECEF2GEO_POLYNOMIAL_UTIL_GEO_HPP
