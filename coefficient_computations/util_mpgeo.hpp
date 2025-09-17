// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the BSD 2-Clause License.

#ifndef ECEF2GEO_POLYNOMIAL_UTIL_MPGEO_HPP
#define ECEF2GEO_POLYNOMIAL_UTIL_MPGEO_HPP

/*
 * Utility functions for multiprecision (mpfr) geodetic calculations.
 *
 * NOTE, ANYONE INCLUDING THIS FILE SHOULD ADD:
 *
 *   wgs84_mpfr_constants util_geo::mpfr_consts = wgs84_mpfr_constants(prec);
 *
 * TO INITIALIZE THE MPFR CONSTANTS TO SOME APPROPRIATE VALUE. THE ARGUMENT prec
 * IS THE DESIRED PRECISION.
 */

#include <mpfr.h>
#include "mplapack/mpreal.h"

#include <functional>  // std::functional.

using namespace mpfr;

struct lla_mpfr {
	mpfr_t lat;
	mpfr_t lon;
	mpfr_t alt;
};

struct xyz_mpfr {
	mpfr_t x;
	mpfr_t y;
	mpfr_t z;
};

struct nva_mpfr {
	xyz_mpfr n;
	mpfr_t alt;
};

struct wgs84_mpfr_constants {
	// Primary WGS84 constants.
	mpfr_t a;  // Semi-major axis.
	mpfr_t inv_f;  // Inverse flattening.

	// Derived WGS84 constants.
	mpfr_t a2;  // Semi-major axis squared.
	mpfr_t a3;
	mpfr_t a4;
	mpfr_t f;  // Flattening.
	mpfr_t b;  // Semi-minor axis / polar radius.
	mpfr_t b2;  // Semi-minor axis / polar radius squared.
	mpfr_t e;  // First eccentricity
	mpfr_t e2;  // First eccentricity squared.
	mpfr_t e4;  // First eccentricity to the power of 4.
	mpfr_t e2_2;
	mpfr_t am2; // Inverse semi-major axis squared.
	mpfr_t b2am2;  // Axis ratio squared.
	mpfr_t b2am4;
	mpfr_t e2bam3;
	mpfr_t bam1;
	mpfr_t bem1;
	mpfr_t neg_bem1;
	mpfr_t neg_e2bm1;
	mpfr_t bam2;

	// Other constants.
	mpfr_t M_1_6;
	mpfr_t M_2_3;
	mpfr_t M_PI_6;

	// Initialization of WGS84 constant.
	wgs84_mpfr_constants(mp_prec_t prec) {
		// These literal values are from the WGS84 standard and are exact.
		// To use GRS 80, instead use "298.257222100882711243" for the inverse
		// flattening. (Or better, compute it will full precision from the
		// underlying constants of GRS 80.)
		mpfr_init_set_str(a, "6378137.0", 10, MPFR_RNDN);
		mpfr_init_set_str(inv_f, "298.257223563", 10, MPFR_RNDN);

		mpfr_init2(a2, prec);
		mpfr_mul(a2, a, a, MPFR_RNDN);
		mpfr_init2(a3, prec);
		mpfr_mul(a3, a2, a, MPFR_RNDN);
		mpfr_init2(a4, prec);
		mpfr_mul(a4, a2, a2, MPFR_RNDN);
		mpfr_init2(am2, prec);
		mpfr_ui_div(am2, 1, a2, MPFR_RNDN);

		mpfr_init2(f, prec);
		mpfr_ui_div(f, 1, inv_f, MPFR_RNDN);

		mpfr_init2(b, prec);
		mpfr_mul(b, f, a, MPFR_RNDN);
		mpfr_sub(b, a, b, MPFR_RNDN);
		mpfr_init2(b2, prec);
		mpfr_mul(b2, b, b, MPFR_RNDN);

		// 1.0 - (b * b) / (a * a)
		mpfr_init2(e2, prec);
		mpfr_div(e2, b2, a2, MPFR_RNDN);
		mpfr_ui_sub(e2, 1, e2, MPFR_RNDN);
		mpfr_init2(e, prec);
		mpfr_sqrt(e, e2, MPFR_RNDN);
		mpfr_init2(e4, prec);
		mpfr_mul(e4, e2, e2, MPFR_RNDN);
		mpfr_init2(e2_2, prec);
		mpfr_div_si(e2_2, e2, 2, MPFR_RNDN);

		mpfr_init2(b2am2, prec);
		mpfr_si_sub(b2am2, 1, e2, MPFR_RNDN);
		mpfr_init2(b2am4, prec);
		mpfr_div(b2am4, b2am2, a2, MPFR_RNDN);
		mpfr_init2(e2bam3, prec);
		mpfr_mul(e2bam3, e2, b, MPFR_RNDN);
		mpfr_div(e2bam3, e2bam3, a3, MPFR_RNDN);
		mpfr_init2(bam1, prec);
		mpfr_div(bam1, b, a, MPFR_RNDN);
		mpfr_init2(bem1, prec);
		mpfr_div(bem1, b, e, MPFR_RNDN);
		mpfr_init2(neg_bem1, prec);
		mpfr_neg(neg_bem1, bem1, MPFR_RNDN);
		mpfr_init2(neg_e2bm1, prec);
		mpfr_div(neg_e2bm1, e2, b, MPFR_RNDN);
		mpfr_neg(neg_e2bm1, neg_e2bm1, MPFR_RNDN);
		mpfr_init2(bam2, prec);
		mpfr_div(bam2, b, a2, MPFR_RNDN);

		mpfr_init2(M_1_6, prec);
		mpfr_set_ui(M_1_6, 1, MPFR_RNDN);
		mpfr_div_si(M_1_6, M_1_6, 6, MPFR_RNDN);
		mpfr_init2(M_2_3, prec);
		mpfr_set_ui(M_2_3, 2, MPFR_RNDN);
		mpfr_div_si(M_2_3, M_2_3, 3, MPFR_RNDN);
		mpfr_init2(M_PI_6, prec);
		mpfr_const_pi(M_PI_6, MPFR_RNDN);
		mpfr_div_ui(M_PI_6, M_PI_6, 6, MPFR_RNDN);
	}
};

// Struct and not namespace to enable having the static mpfr constants.
struct util_mpgeo {

	// Statically allocated mpfr constants.
	static wgs84_mpfr_constants mpfr_consts;

	/** Computes latitude and altitude from an ECEF coordinate.
	 *
	 * The computations are done according to Vermielle's 2011 exact method.
	 *
	 * @param lla latitude-longitude-altitude output struct.
	 * @param ecef Ecef coordinate.
	 */
	static inline void ecef2lla(lla_mpfr &lla, const xyz_mpfr &ecef) {
		mp_prec_t prec = mpfr_get_prec(lla.lat);
		mpfr_t y2, t2, t, z2, q, p, r, r38, ev;
		mpfr_init2(y2, prec);
		mpfr_init2(t2, prec);
		mpfr_init2(t, prec);
		mpfr_init2(z2, prec);
		mpfr_init2(q, prec);
		mpfr_init2(p, prec);
		mpfr_init2(r, prec);
		mpfr_init2(r38, prec);
		mpfr_init2(ev, prec);
		mpfr_sqr(y2, ecef.y, MPFR_RNDN);
		mpfr_sqr(t2, ecef.x, MPFR_RNDN);
		mpfr_add(t2, t2, y2, MPFR_RNDN);
		mpfr_sqrt(t, t2, MPFR_RNDN);
		mpfr_sqr(z2, ecef.z, MPFR_RNDN);
		mpfr_mul(q, z2, mpfr_consts.b2am4, MPFR_RNDN);
		mpfr_mul(p, t2, mpfr_consts.am2, MPFR_RNDN);
		mpfr_add(r, p, q, MPFR_RNDN);
		mpfr_sub(r, r, mpfr_consts.e4, MPFR_RNDN);
		mpfr_div_ui(r, r, 6, MPFR_RNDN);
		mpfr_sqr(r38, r, MPFR_RNDN);
		mpfr_mul(r38, r38, r, MPFR_RNDN);
		mpfr_mul_ui(r38, r38, 8, MPFR_RNDN);
		mpfr_mul(ev, mpfr_consts.e4, p, MPFR_RNDN);
		mpfr_mul(ev, ev, q, MPFR_RNDN);
		mpfr_add(ev, ev, r38, MPFR_RNDN);

		mpfr_atan2(lla.lon, ecef.y, ecef.x, MPFR_RNDN);

		if (mpfr_sgn(ev) > 0) {
			mpfr_t s, sqrt_ev, c1, c2, u, u2, v, upv, w, k, l, D, d;
			mpfr_init2(s, prec);
			mpfr_init2(sqrt_ev, prec);
			mpfr_init2(c1, prec);
			mpfr_init2(c2, prec);
			mpfr_init2(u, prec);
			mpfr_init2(u2, prec);
			mpfr_init2(v, prec);
			mpfr_init2(upv, prec);
			mpfr_init2(w, prec);
			mpfr_init2(k, prec);
			mpfr_init2(l, prec);
			mpfr_init2(D, prec);
			mpfr_init2(d, prec);

			mpfr_abs(s, ecef.z, MPFR_RNDN);
			mpfr_mul(s, mpfr_consts.e2bam3, s, MPFR_RNDN);
			mpfr_mul(s, s, t, MPFR_RNDN);
			mpfr_sqrt(sqrt_ev, ev, MPFR_RNDN);
			mpfr_add(c1, sqrt_ev, s, MPFR_RNDN);
			mpfr_sub(c2, sqrt_ev, s, MPFR_RNDN);
			mpfr_sqr(c1, c1, MPFR_RNDN);
			mpfr_cbrt(c1, c1, MPFR_RNDN);
			mpfr_div_ui(c1, c1, 2, MPFR_RNDN);
			mpfr_sqr(c2, c2, MPFR_RNDN);
			mpfr_cbrt(c2, c2, MPFR_RNDN);
			mpfr_div_ui(c2, c2, 2, MPFR_RNDN);
			mpfr_add(u, r, c1, MPFR_RNDN);
			mpfr_add(u, u, c2, MPFR_RNDN);

			mpfr_sqr(u2, u, MPFR_RNDN);
			mpfr_mul(v, mpfr_consts.e4, q, MPFR_RNDN);
			mpfr_add(v, v, u2, MPFR_RNDN);
			mpfr_sqrt(v, v, MPFR_RNDN);
			mpfr_add(upv, u, v, MPFR_RNDN);
			mpfr_add(w, u, v, MPFR_RNDN);
			mpfr_sub(w, w, q, MPFR_RNDN);
			mpfr_div(w, w, v, MPFR_RNDN);
			mpfr_mul(w, w, mpfr_consts.e2_2, MPFR_RNDN);
			mpfr_sqr(k, w, MPFR_RNDN);
			mpfr_add(k, k, v, MPFR_RNDN);
			mpfr_add(k, k, u, MPFR_RNDN);
			mpfr_sqrt(k, k, MPFR_RNDN);
			mpfr_add(k, k, w, MPFR_RNDN);
			mpfr_div(k, upv, k, MPFR_RNDN);
			mpfr_add(l, k, mpfr_consts.e2, MPFR_RNDN);
			mpfr_div(l, k, l, MPFR_RNDN);
			mpfr_mul(D, l, t, MPFR_RNDN);
			mpfr_sqr(d, D, MPFR_RNDN);
			mpfr_add(d, d, z2, MPFR_RNDN);
			mpfr_sqrt(d, d, MPFR_RNDN);

			mpfr_sub(lla.alt, k, mpfr_consts.b2am2, MPFR_RNDN);
			mpfr_mul(lla.alt, lla.alt, d, MPFR_RNDN);
			mpfr_div(lla.alt, lla.alt, k, MPFR_RNDN);

			mpfr_add(lla.lat, D, d, MPFR_RNDN);
			mpfr_div(lla.lat, ecef.z, lla.lat, MPFR_RNDN);
			mpfr_atan(lla.lat, lla.lat, MPFR_RNDN);
			mpfr_mul_ui(lla.lat, lla.lat, 2, MPFR_RNDN);

			mpfr_clear(s);
			mpfr_clear(sqrt_ev);
			mpfr_clear(c1);
			mpfr_clear(c2);
			mpfr_clear(u);
			mpfr_clear(u2);
			mpfr_clear(v);
			mpfr_clear(upv);
			mpfr_clear(w);
			mpfr_clear(k);
			mpfr_clear(l);
			mpfr_clear(D);
			mpfr_clear(d);
		} else if (!mpfr_zero_p(q)) {
			mpfr_t s, up, sqrt_r38, sin_up, u, u2, v, upv, w, k, l, D, d;
			mpfr_init2(s, prec);
			mpfr_init2(up, prec);
			mpfr_init2(sqrt_r38, prec);
			mpfr_init2(sin_up, prec);
			mpfr_init2(u, prec);
			mpfr_init2(u2, prec);
			mpfr_init2(v, prec);
			mpfr_init2(upv, prec);
			mpfr_init2(w, prec);
			mpfr_init2(k, prec);
			mpfr_init2(l, prec);
			mpfr_init2(D, prec);
			mpfr_init2(d, prec);

			mpfr_abs(s, ecef.z, MPFR_RNDN);
			mpfr_mul(s, mpfr_consts.e2bam3, s, MPFR_RNDN);
			mpfr_mul(s, s, t, MPFR_RNDN);
			mpfr_neg(sqrt_r38, r38, MPFR_RNDN);
			mpfr_sqrt(sqrt_r38, sqrt_r38, MPFR_RNDN);
			mpfr_neg(up, ev, MPFR_RNDN);
			mpfr_sqrt(up, up, MPFR_RNDN);
			mpfr_add(up, up, sqrt_r38, MPFR_RNDN);
			mpfr_div(up, s, up, MPFR_RNDN);
			mpfr_atan(up, up, MPFR_RNDN);
			mpfr_mul(up, mpfr_consts.M_2_3, up, MPFR_RNDN);
			mpfr_sin(sin_up, up, MPFR_RNDN);
			mpfr_add(u, mpfr_consts.M_PI_6, up, MPFR_RNDN);
			mpfr_cos(u, u, MPFR_RNDN);
			mpfr_mul(u, sin_up, u, MPFR_RNDN);
			mpfr_mul(u, u, r, MPFR_RNDN);
			mpfr_mul_si(u, u, -4, MPFR_RNDN);

			mpfr_sqr(u2, u, MPFR_RNDN);
			mpfr_mul(v, mpfr_consts.e4, q, MPFR_RNDN);
			mpfr_add(v, v, u2, MPFR_RNDN);
			mpfr_sqrt(v, v, MPFR_RNDN);
			mpfr_add(upv, u, v, MPFR_RNDN);
			mpfr_add(w, u, v, MPFR_RNDN);
			mpfr_sub(w, w, q, MPFR_RNDN);
			mpfr_div(w, w, v, MPFR_RNDN);
			mpfr_mul(w, w, mpfr_consts.e2_2, MPFR_RNDN);
			mpfr_sqr(k, w, MPFR_RNDN);
			mpfr_add(k, k, v, MPFR_RNDN);
			mpfr_add(k, k, u, MPFR_RNDN);
			mpfr_sqrt(k, k, MPFR_RNDN);
			mpfr_add(k, k, w, MPFR_RNDN);
			mpfr_div(k, upv, k, MPFR_RNDN);
			mpfr_add(l, k, mpfr_consts.e2, MPFR_RNDN);
			mpfr_div(l, k, l, MPFR_RNDN);
			mpfr_mul(D, l, t, MPFR_RNDN);
			mpfr_sqr(d, D, MPFR_RNDN);
			mpfr_add(d, d, z2, MPFR_RNDN);
			mpfr_sqrt(d, d, MPFR_RNDN);

			mpfr_sub(lla.alt, k, mpfr_consts.b2am2, MPFR_RNDN);
			mpfr_mul(lla.alt, lla.alt, d, MPFR_RNDN);
			mpfr_div(lla.alt, lla.alt, k, MPFR_RNDN);

			mpfr_add(lla.lat, D, d, MPFR_RNDN);
			mpfr_div(lla.lat, ecef.z, lla.lat, MPFR_RNDN);
			mpfr_atan(lla.lat, lla.lat, MPFR_RNDN);
			mpfr_mul_ui(lla.lat, lla.lat, 2, MPFR_RNDN);

			mpfr_clear(s);
			mpfr_clear(up);
			mpfr_clear(sqrt_r38);
			mpfr_clear(sin_up);
			mpfr_clear(u);
			mpfr_clear(u2);
			mpfr_clear(v);
			mpfr_clear(upv);
			mpfr_clear(w);
			mpfr_clear(k);
			mpfr_clear(l);
			mpfr_clear(D);
			mpfr_clear(d);
		} else {
			// On the singular disc (including center of earth).
			// Values are taken to have positive latitude.
			mpfr_t hp, sqrt_e4_p;
			mpfr_init2(hp, prec);
			mpfr_init2(sqrt_e4_p, prec);

			mpfr_sub(lla.alt, mpfr_consts.e2, p, MPFR_RNDN);
			mpfr_sqrt(lla.alt, lla.alt, MPFR_RNDN);
			mpfr_mul(lla.alt, lla.alt, mpfr_consts.neg_bem1, MPFR_RNDN);
			mpfr_mul(hp, lla.alt, mpfr_consts.neg_e2bm1, MPFR_RNDN);
			mpfr_sub(sqrt_e4_p, mpfr_consts.e4, p, MPFR_RNDN);
			mpfr_sqrt(sqrt_e4_p, sqrt_e4_p, MPFR_RNDN);

			mpfr_sqrt(lla.lat, p, MPFR_RNDN);
			mpfr_mul(lla.lat, lla.lat, mpfr_consts.bam1, MPFR_RNDN);
			mpfr_add(lla.lat, hp, lla.lat, MPFR_RNDN);
			mpfr_div(lla.lat, sqrt_e4_p, lla.lat, MPFR_RNDN);
			mpfr_atan(lla.lat, lla.lat, MPFR_RNDN);
			mpfr_mul_ui(lla.lat, lla.lat, 2, MPFR_RNDN);

			mpfr_clear(hp);
			mpfr_clear(sqrt_e4_p);
		}
		mpfr_clear(t2);
		mpfr_clear(t);
		mpfr_clear(z2);
		mpfr_clear(y2);
		mpfr_clear(q);
		mpfr_clear(p);
		mpfr_clear(r);
		mpfr_clear(r38);
		mpfr_clear(ev);
	}

	static inline void ecef2nva(nva_mpfr &nva, const xyz_mpfr &ecef) {
		mp_prec_t prec = mpfr_get_prec(nva.n.x);
		mpfr_t y2, t2, t, z2, q, p, r, r38, ev;
		mpfr_init2(y2, prec);
		mpfr_init2(t2, prec);
		mpfr_init2(t, prec);
		mpfr_init2(z2, prec);
		mpfr_init2(q, prec);
		mpfr_init2(p, prec);
		mpfr_init2(r, prec);
		mpfr_init2(r38, prec);
		mpfr_init2(ev, prec);

		mpfr_sqr(y2, ecef.y, MPFR_RNDN);
		mpfr_sqr(t2, ecef.x, MPFR_RNDN);
		mpfr_add(t2, t2, y2, MPFR_RNDN);
		mpfr_sqrt(t, t2, MPFR_RNDN);
		mpfr_sqr(z2, ecef.z, MPFR_RNDN);
		mpfr_mul(q, z2, mpfr_consts.b2am4, MPFR_RNDN);
		mpfr_mul(p, t2, mpfr_consts.am2, MPFR_RNDN);
		mpfr_add(r, p, q, MPFR_RNDN);
		mpfr_sub(r, r, mpfr_consts.e4, MPFR_RNDN);
		mpfr_div_ui(r, r, 6, MPFR_RNDN);
		mpfr_sqr(r38, r, MPFR_RNDN);
		mpfr_mul(r38, r38, r, MPFR_RNDN);
		mpfr_mul_ui(r38, r38, 8, MPFR_RNDN);
		mpfr_mul(ev, mpfr_consts.e4, p, MPFR_RNDN);
		mpfr_mul(ev, ev, q, MPFR_RNDN);
		mpfr_add(ev, ev, r38, MPFR_RNDN);

		if (mpfr_sgn(ev) > 0) {
			mpfr_t s, sqrt_ev, c1, c2, u, u2, v, upv, w, k, l, D, d, n, nl;
			mpfr_init2(s, prec);
			mpfr_init2(sqrt_ev, prec);
			mpfr_init2(c1, prec);
			mpfr_init2(c2, prec);
			mpfr_init2(u, prec);
			mpfr_init2(u2, prec);
			mpfr_init2(v, prec);
			mpfr_init2(upv, prec);
			mpfr_init2(w, prec);
			mpfr_init2(k, prec);
			mpfr_init2(l, prec);
			mpfr_init2(D, prec);
			mpfr_init2(d, prec);
			mpfr_init2(n, prec);
			mpfr_init2(nl, prec);

			mpfr_abs(s, ecef.z, MPFR_RNDN);
			mpfr_mul(s, mpfr_consts.e2bam3, s, MPFR_RNDN);
			mpfr_mul(s, s, t, MPFR_RNDN);
			mpfr_sqrt(sqrt_ev, ev, MPFR_RNDN);
			mpfr_add(c1, sqrt_ev, s, MPFR_RNDN);
			mpfr_sub(c2, sqrt_ev, s, MPFR_RNDN);
			mpfr_sqr(c1, c1, MPFR_RNDN);
			mpfr_cbrt(c1, c1, MPFR_RNDN);
			mpfr_div_ui(c1, c1, 2, MPFR_RNDN);
			mpfr_sqr(c2, c2, MPFR_RNDN);
			mpfr_cbrt(c2, c2, MPFR_RNDN);
			mpfr_div_ui(c2, c2, 2, MPFR_RNDN);
			mpfr_add(u, r, c1, MPFR_RNDN);
			mpfr_add(u, u, c2, MPFR_RNDN);

			mpfr_sqr(u2, u, MPFR_RNDN);
			mpfr_mul(v, mpfr_consts.e4, q, MPFR_RNDN);
			mpfr_add(v, v, u2, MPFR_RNDN);
			mpfr_sqrt(v, v, MPFR_RNDN);
			mpfr_add(upv, u, v, MPFR_RNDN);
			mpfr_add(w, u, v, MPFR_RNDN);
			mpfr_sub(w, w, q, MPFR_RNDN);
			mpfr_div(w, w, v, MPFR_RNDN);
			mpfr_mul(w, w, mpfr_consts.e2_2, MPFR_RNDN);
			mpfr_sqr(k, w, MPFR_RNDN);
			mpfr_add(k, k, v, MPFR_RNDN);
			mpfr_add(k, k, u, MPFR_RNDN);
			mpfr_sqrt(k, k, MPFR_RNDN);
			mpfr_add(k, k, w, MPFR_RNDN);
			mpfr_div(k, upv, k, MPFR_RNDN);
			mpfr_add(l, k, mpfr_consts.e2, MPFR_RNDN);
			mpfr_div(l, k, l, MPFR_RNDN);
			mpfr_mul(D, l, t, MPFR_RNDN);
			mpfr_sqr(d, D, MPFR_RNDN);
			mpfr_add(d, d, z2, MPFR_RNDN);
			mpfr_sqrt(d, d, MPFR_RNDN);

			mpfr_sub(nva.alt, k, mpfr_consts.b2am2, MPFR_RNDN);
			mpfr_mul(nva.alt, nva.alt, d, MPFR_RNDN);
			mpfr_div(nva.alt, nva.alt, k, MPFR_RNDN);
			mpfr_ui_div(n, 1, d, MPFR_RNDN);
			mpfr_mul(nl, n, l, MPFR_RNDN);
			mpfr_mul(nva.n.x, nl, ecef.x, MPFR_RNDN);
			mpfr_mul(nva.n.y, nl, ecef.y, MPFR_RNDN);
			mpfr_mul(nva.n.z, n, ecef.z, MPFR_RNDN);

			mpfr_clear(s);
			mpfr_clear(sqrt_ev);
			mpfr_clear(c1);
			mpfr_clear(c2);
			mpfr_clear(u);
			mpfr_clear(u2);
			mpfr_clear(v);
			mpfr_clear(upv);
			mpfr_clear(w);
			mpfr_clear(k);
			mpfr_clear(l);
			mpfr_clear(D);
			mpfr_clear(d);
			mpfr_clear(n);
			mpfr_clear(nl);
		} else if (!mpfr_zero_p(q)) {
			mpfr_t s, up, sqrt_r38, sin_up, u, u2, v, upv, w, k, l, D, d, n, nl;
			mpfr_init2(s, prec);
			mpfr_init2(up, prec);
			mpfr_init2(sqrt_r38, prec);
			mpfr_init2(sin_up, prec);
			mpfr_init2(u, prec);
			mpfr_init2(u2, prec);
			mpfr_init2(v, prec);
			mpfr_init2(upv, prec);
			mpfr_init2(w, prec);
			mpfr_init2(k, prec);
			mpfr_init2(l, prec);
			mpfr_init2(D, prec);
			mpfr_init2(d, prec);
			mpfr_init2(n, prec);
			mpfr_init2(nl, prec);

			mpfr_abs(s, ecef.z, MPFR_RNDN);
			mpfr_mul(s, mpfr_consts.e2bam3, s, MPFR_RNDN);
			mpfr_mul(s, s, t, MPFR_RNDN);
			mpfr_neg(sqrt_r38, r38, MPFR_RNDN);
			mpfr_sqrt(sqrt_r38, sqrt_r38, MPFR_RNDN);
			mpfr_neg(up, ev, MPFR_RNDN);
			mpfr_sqrt(up, up, MPFR_RNDN);
			mpfr_add(up, up, sqrt_r38, MPFR_RNDN);
			mpfr_div(up, s, up, MPFR_RNDN);
			mpfr_atan(up, up, MPFR_RNDN);
			mpfr_mul(up, mpfr_consts.M_2_3, up, MPFR_RNDN);
			mpfr_sin(sin_up, up, MPFR_RNDN);
			mpfr_add(u, mpfr_consts.M_PI_6, up, MPFR_RNDN);
			mpfr_cos(u, u, MPFR_RNDN);
			mpfr_mul(u, sin_up, u, MPFR_RNDN);
			mpfr_mul(u, u, r, MPFR_RNDN);
			mpfr_mul_si(u, u, -4, MPFR_RNDN);

			mpfr_sqr(u2, u, MPFR_RNDN);
			mpfr_mul(v, mpfr_consts.e4, q, MPFR_RNDN);
			mpfr_add(v, v, u2, MPFR_RNDN);
			mpfr_sqrt(v, v, MPFR_RNDN);
			mpfr_add(upv, u, v, MPFR_RNDN);
			mpfr_add(w, u, v, MPFR_RNDN);
			mpfr_sub(w, w, q, MPFR_RNDN);
			mpfr_div(w, w, v, MPFR_RNDN);
			mpfr_mul(w, w, mpfr_consts.e2_2, MPFR_RNDN);
			mpfr_sqr(k, w, MPFR_RNDN);
			mpfr_add(k, k, v, MPFR_RNDN);
			mpfr_add(k, k, u, MPFR_RNDN);
			mpfr_sqrt(k, k, MPFR_RNDN);
			mpfr_add(k, k, w, MPFR_RNDN);
			mpfr_div(k, upv, k, MPFR_RNDN);
			mpfr_add(l, k, mpfr_consts.e2, MPFR_RNDN);
			mpfr_div(l, k, l, MPFR_RNDN);
			mpfr_mul(D, l, t, MPFR_RNDN);
			mpfr_sqr(d, D, MPFR_RNDN);
			mpfr_add(d, d, z2, MPFR_RNDN);
			mpfr_sqrt(d, d, MPFR_RNDN);

			mpfr_sub(nva.alt, k, mpfr_consts.b2am2, MPFR_RNDN);
			mpfr_mul(nva.alt, nva.alt, d, MPFR_RNDN);
			mpfr_div(nva.alt, nva.alt, k, MPFR_RNDN);
			mpfr_ui_div(n, 1, d, MPFR_RNDN);
			mpfr_mul(nl, n, l, MPFR_RNDN);
			mpfr_mul(nva.n.x, nl, ecef.x, MPFR_RNDN);
			mpfr_mul(nva.n.y, nl, ecef.y, MPFR_RNDN);
			mpfr_mul(nva.n.z, n, ecef.z, MPFR_RNDN);

			mpfr_clear(s);
			mpfr_clear(up);
			mpfr_clear(sqrt_r38);
			mpfr_clear(sin_up);
			mpfr_clear(u);
			mpfr_clear(u2);
			mpfr_clear(v);
			mpfr_clear(upv);
			mpfr_clear(w);
			mpfr_clear(k);
			mpfr_clear(l);
			mpfr_clear(D);
			mpfr_clear(d);
			mpfr_clear(n);
			mpfr_clear(nl);
		} else {
			// On the singular disc (including center of earth).
			// Values are taken to have positive latitude.
			mpfr_t hp, sqrt_e4_p, n, bam2_n;
			mpfr_init2(hp, prec);
			mpfr_init2(sqrt_e4_p, prec);
			mpfr_init2(n, prec);
			mpfr_init2(bam2_n, prec);

			mpfr_sub(nva.alt, mpfr_consts.e2, p, MPFR_RNDN);
			mpfr_sqrt(nva.alt, nva.alt, MPFR_RNDN);
			mpfr_mul(nva.alt, nva.alt, mpfr_consts.neg_bem1, MPFR_RNDN);
			mpfr_sub(sqrt_e4_p, mpfr_consts.e4, p, MPFR_RNDN);
			mpfr_sqrt(sqrt_e4_p, sqrt_e4_p, MPFR_RNDN);

			mpfr_mul(n, nva.alt, mpfr_consts.neg_e2bm1, MPFR_RNDN);
			mpfr_ui_div(n, 1, n, MPFR_RNDN);
			mpfr_mul(bam2_n, mpfr_consts.bam2, n, MPFR_RNDN);
			mpfr_mul(nva.n.x, bam2_n, ecef.x, MPFR_RNDN);
			mpfr_mul(nva.n.y, bam2_n, ecef.y, MPFR_RNDN);
			mpfr_mul(nva.n.z, sqrt_e4_p, n, MPFR_RNDN);

			mpfr_clear(hp);
			mpfr_clear(sqrt_e4_p);
			mpfr_clear(n);
			mpfr_clear(bam2_n);
		}
		mpfr_clear(t2);
		mpfr_clear(t);
		mpfr_clear(z2);
		mpfr_clear(y2);
		mpfr_clear(q);
		mpfr_clear(p);
		mpfr_clear(r);
		mpfr_clear(r38);
		mpfr_clear(ev);
	}

	/** Computes latitude from spherical latitude and altitude.
	 *
	 * The computations are done according to Vermielle's 2011 method.
	 *
	 * @param lat_c spherical latitude.
	 * @param h_c spherical altitude.
	 * @param h_0 Reference for the spherical altitude.
	 * @param prec Precision of the arithmetics.
	 * @return the latitude.
	 */
	static inline mpreal f_c(const mpreal &lat_c, const mpreal &h_c, const mpreal &h_0,
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

		xyz_mpfr ecef;
		mpfr_init2(ecef.x, prec);
		mpfr_init2(ecef.y, prec);
		mpfr_init2(ecef.z, prec);

		mpfr_ptr r_ptr(r);
		mpfr_set(ecef.x, r_ptr, MPFR_RNDN);
		mpfr_set_d(ecef.y, 0.0, MPFR_RNDN);
		mpfr_ptr z_ptr(z);
		mpfr_set(ecef.z, z_ptr, MPFR_RNDN);

		// Using a transformation for only latitude/altitude made this code
		// run roughly 20% faster. It was judged not to be worse it.
		lla_mpfr la;
		mpfr_init2(la.lat, prec);
		mpfr_init2(la.alt, prec);
		mpfr_init2(la.lon, prec);
		util_mpgeo::ecef2lla(la, ecef);
		mpreal res(la.lat);
		mpfr_clear(la.lat);
		mpfr_clear(la.alt);

		mpfr_clear(ecef.x);
		mpfr_clear(ecef.y);
		mpfr_clear(ecef.z);
		return res;
	}

	/** Computes altitude from spherical latitude and altitude.
	 *
	 * The computations are done according to Vermielle's 2011 method.
	 *
	 * @param lat_c spherical latitude.
	 * @param h_c spherical altitude.
	 * @param h_0 Reference for the spherical altitude.
	 * @param prec Precision of the arithmetics.
	 * @return the altitude.
	 */
	static inline mpreal g_c(const mpreal &lat_c, const mpreal &h_c, const mpreal &h_0,
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

		xyz_mpfr ecef;
		mpfr_init2(ecef.x, prec);
		mpfr_init2(ecef.y, prec);
		mpfr_init2(ecef.z, prec);

		mpfr_ptr r_ptr(r);
		mpfr_set(ecef.x, r_ptr, MPFR_RNDN);
		mpfr_set_d(ecef.y, 0.0, MPFR_RNDN);
		mpfr_ptr z_ptr(z);
		mpfr_set(ecef.z, z_ptr, MPFR_RNDN);

		// Using a transformation for only latitude/altitude made this code
		// run roughly 20% faster. It was judged not to be worse it.
		lla_mpfr la;
		mpfr_init2(la.lat, prec);
		mpfr_init2(la.alt, prec);
		mpfr_init2(la.lon, prec);
		util_mpgeo::ecef2lla(la, ecef);
		mpreal res(la.alt);
		mpfr_clear(la.lat);
		mpfr_clear(la.alt);

		mpfr_clear(ecef.x);
		mpfr_clear(ecef.y);
		mpfr_clear(ecef.z);
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

#endif // ECEF2GEO_POLYNOMIAL_UTIL_MPGEO_HPP
