// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the BSD 2-Clause License.
/*
 * Computes and generates code for minimax multiplicative correction for ECEF
 * to geodetic coordinate transformations. Minimax approximations are computed
 * with Sollya and exported to header files.
 */

// Include these first since they are needed in subsequent includes.
// Included such that we can switch between mpreal.h and mplapack/mpreal.h.

// Defines to give an expected behaviour of mpfr.h for mpreal.h.
#define MPFR_USE_NO_MACRO
#define MPFR_USE_INTMAX_T
#include <cstdint>  // Required for mpreal.h. when mpfr is included first.
#include <mpfr.h>
#include <mpreal.h>

// ***** To handle inconsistencies between mplapack/mpreal.h and mpreal.h *****
void mpfr_set_mpreal(mpfr_ptr to, mpfr::mpreal &val) {
	mpfr_set(to, val.mpfr_ptr(), MPFR_RNDN);
}
inline mpfr::mpreal sum(const mpfr::mpreal tab[], unsigned long int n) {
	int status;
	return mpfr::sum(tab, n, status, mpfr::mpreal::get_default_rnd());
}
// ****************************************************************************

#include "settings.hpp"
#include "util_eval.hpp"
#include "remez_sollya.hpp"
#include "util_root.hpp"  // find_max
#include "util_mpgeo.hpp"  // f_c

#define _USE_MATH_DEFINES
#include <cmath>  // std::sqrt, std::cbrt, std::atan, std::sin, std::cos.
#include <vector>  // std::vector.
#include <iostream>  // std::cout.
#include <limits>  // std::numeric_limits.
#include <ostream>  // std::endl.
#include <string>  // std::basic_string

// Just to make the code more compact.
using str_vec = std::vector<std::string>;

// Bits of precision for Sollya amd mpfr computations.
constexpr mp_prec_t prec = 200;
// Mpfr consts with appropriate precision.
wgs84_mpfr_constants util_mpgeo::mpfr_consts = wgs84_mpfr_constants(prec);

static const char *sigma_preamble =
		"/*\n"
		" * Generated. DO NOT EDIT.\n"
		" */\n\n"
		"#ifndef ECEF2GEO_SIGMA_HPP\n"
		"#define ECEF2GEO_SIGMA_HPP\n\n"
		"namespace ecef2geo {\n\n";
static const char *sigma_proto =
		"namespace priv {\n\n"
		"/** Minimax approximations of sigma.\n"
		" *\n"
		" * For details see:\n"
		" * Nilsson, J.-O. Minimax polynomial ECEF to geodetic coordinate transformation\n"
		" * approximations. doi: https://doi.org/10.48550/arXiv.2405.06572\n"
		" *\n"
		" * @paramt L Total degree of the multiplicative polynomial correction.\n"
		" * @param d (phi-phi_c)^2.\n"
		" * @return Minimax approximation of sigma.\n"
		" */\n"
		"template<int L> inline double sigma(double d) = delete;  // Only allow provided specializations.\n";
static const char *sigma_signature = "template<> inline double sigma<%u>(double d) {";
static const char *sigma_postamble =
		"\n}  // namespace ecef2geo\n"
		"}  // namespace priv\n"
		"\n#endif // ECEF2GEO_SIGMA_HPP\n";

static const char *tau_preamble =
		"/*\n"
		" * Generated. DO NOT EDIT.\n"
		" */\n\n"
		"#ifndef ECEF2GEO_TAU_HPP\n"
		"#define ECEF2GEO_TAU_HPP\n\n"
		"namespace ecef2geo {\n\n";
static const char *tau_proto =
		"namespace priv {\n\n"
		"/** Minimax approximations of tau.\n"
		" *\n"
		" * Nilsson, J.-O. Minimax polynomial ECEF to geodetic coordinate transformation\n"
		" * approximations. doi: https://doi.org/10.48550/arXiv.2405.06572\n"
		" *\n"
		" * @paramt L Total degree of the multiplicative polynomial correction.\n"
		" * @param d (phi-phi_c)^2.\n"
		" * @return Minimax approximation of tau.\n"
		" */\n"
		"template<int L> inline double tau(double w, double d) = delete;  // Only allow provided specializations.\n";
static const char *tau_signature = "template<> inline double tau<%u>(double w, double d) {";
static const char *tau_postamble =
		"\n}  // namespace ecef2geo\n"
		"}  // namespace priv\n"
		"\n#endif // ECEF2GEO_TAU_HPP\n";

/** Divide one integer by another, rounding towards minus infinity.
 *
 * http://www.microhowto.info/howto/round_towards_minus_infinity_when_dividing_integers_in_c_or_c++.html
 *
 * @param x the dividend
 * @param y the divisor
 * @return the quotient, rounded towards minus infinity
 */
int div_floor(int x, int y) {
	int q = x / y;
	int r = x % y;
	if ((r != 0) && ((r < 0) != (y < 0))) --q;
	return q;
}

int main() {
	// Compute the value range of the approximations.
	// The minimum delta is by definition zero, at the poles.
	const double DELTA_MIN = 0.0;
	// In contrast, the maximum delta occurs at the lowest altitude.
	mpreal h_c = mpreal(ALT_LO_LIMIT, prec) + mpreal(util_mpgeo::mpfr_consts.b);
	mpreal h_0 = mpreal(ALT_REFERENCE, prec);
	// Numerical differentiation difference.
	// Selecting a proper "h" is a non-trivial problem.
	// See https://en.wikipedia.org/wiki/Numerical_differentiation.
	// Selecting a too low value will make the algorithm not converge.
	mpreal h;
	h.set_prec(prec);
	h = mul_2si(const_pi(prec) / 2, -(prec / 4), MPFR_RNDN);
	// Squared latitude to geocentric latitude difference.
	std::function<mpreal(const mpreal &)> f = [&h_c, &h_0](const mpreal &lat) {
		mpreal diff = util_mpgeo::f_c(lat, h_c, h_0, prec) - lat;
		return diff * diff;
	};
	mpreal DELTA_MAX_LAT = find_max(f, h, 0.0, const_pi(prec) * 2 - 0.1, prec);
	mpreal DELTA_MAX_VAL;
	DELTA_MAX_VAL.set_prec(prec);
	DELTA_MAX_VAL = f(DELTA_MAX_LAT);
	std::basic_string<char> DELTA_MAX_STR = DELTA_MAX_VAL.toString(17);
	const char* DELTA_MAX = DELTA_MAX_STR.c_str();

	// Initialize Sollya.
	sollya_lib_init();
	sollya_obj_t sollya_prec;
	sollya_prec = SOLLYA_CONST_SI64(prec);
	sollya_lib_set_prec(sollya_prec);

	// Open files and write preamble, ranges and protos for the approximations.
	FILE *fp_sig = fopen((CODEGEN_FOLDER + "sigma.hpp").c_str(), "w");
	fprintf(fp_sig, "%s", sigma_preamble);
	fprintf(fp_sig, "// Minimax approximation ranges.\n");
	fprintf(fp_sig, "constexpr double SIGMA_H_MIN = %.17e;\n", ALT_LO_LIMIT);
	fprintf(fp_sig, "constexpr double SIGMA_H_MAX = %.17e;\n", ALT_HI_LIMIT);
	fprintf(fp_sig, "constexpr double SIGMA_DELTA_MIN = %.17e;\n", DELTA_MIN);
	fprintf(fp_sig, "constexpr double SIGMA_DELTA_MAX = %s;\n\n", DELTA_MAX);
	fprintf(fp_sig, "%s", sigma_proto);
	FILE *fp_tau = fopen((CODEGEN_FOLDER + "tau.hpp").c_str(), "w");
	fprintf(fp_tau, "%s", tau_preamble);
	fprintf(fp_tau, "// Minimax approximation ranges.\n");
	fprintf(fp_tau, "constexpr double TAU_H_MIN = %.17e;\n", ALT_LO_LIMIT);
	fprintf(fp_tau, "constexpr double TAU_H_MAX = %.17e;\n", ALT_HI_LIMIT);
	fprintf(fp_tau, "constexpr double TAU_DELTA_MIN = %.17e;\n", DELTA_MIN);
	fprintf(fp_tau, "constexpr double TAU_DELTA_MAX = %s;\n\n", DELTA_MAX);
	fprintf(fp_tau, "%s", tau_proto);

	// Minimax approximation arguments.
	sollya_obj_t range, low, high, sigma, tau;
	low = sollya_lib_constant_from_double(DELTA_MIN);
	high = sollya_lib_parse_string(DELTA_MAX_VAL.toString().c_str());
	range = sollya_lib_range(low, high);
	// Truncated series approximations, well below double precision accuracy,
	// of sigma and tau. Note that the variable of sigma and tau is small
	// (below DELTA_MAX ~ 1.13*10^-5) to start with so the terms fall off
	// quickly.
	// @formatter:off (Disable reformate for CLion.)
	sigma = SOLLYA_ADD(
				SOLLYA_SUB(
					SOLLYA_ADD(
						SOLLYA_SUB(
							SOLLYA_ADD(
								SOLLYA_SUB(
									SOLLYA_ADD(
										SOLLYA_SUB(
											SOLLYA_CONST(1), // 1
											SOLLYA_MUL(SOLLYA_DIV(SOLLYA_CONST(1), SOLLYA_CONST(2)), SOLLYA_X_)),  // - 1/2*x
										SOLLYA_MUL(SOLLYA_DIV(SOLLYA_CONST(1), SOLLYA_CONST(24)), SOLLYA_POW(SOLLYA_X_, SOLLYA_CONST(2)))),  // + 1/24*x^2
									SOLLYA_MUL(SOLLYA_DIV(SOLLYA_CONST(1), SOLLYA_CONST(720)), SOLLYA_POW(SOLLYA_X_, SOLLYA_CONST(3)))),  // - 1/720*x^3
								SOLLYA_MUL(SOLLYA_DIV(SOLLYA_CONST(1), SOLLYA_CONST(40320)), SOLLYA_POW(SOLLYA_X_, SOLLYA_CONST(4)))),  // + 1/40320*x^4
							SOLLYA_MUL(SOLLYA_DIV(SOLLYA_CONST(1), SOLLYA_CONST(3628800)), SOLLYA_POW(SOLLYA_X_, SOLLYA_CONST(5)))),  // - 1/3628800*x^5
						SOLLYA_MUL(SOLLYA_DIV(SOLLYA_CONST(1), SOLLYA_CONST(479001600)), SOLLYA_POW(SOLLYA_X_, SOLLYA_CONST(6)))),  // + 1/479001600*x^6
					SOLLYA_MUL(SOLLYA_DIV(SOLLYA_CONST(1), SOLLYA_CONST(87178291200)), SOLLYA_POW(SOLLYA_X_, SOLLYA_CONST(7)))),  // - 1/87178291200*x^7
				SOLLYA_MUL(SOLLYA_DIV(SOLLYA_CONST(1), SOLLYA_CONST(20922789888000)), SOLLYA_POW(SOLLYA_X_, SOLLYA_CONST(8))));  // + 1/20922789888000*x^8
	tau = SOLLYA_ADD(
			SOLLYA_SUB(
				SOLLYA_ADD(
					SOLLYA_SUB(
						SOLLYA_ADD(
							SOLLYA_SUB(
								SOLLYA_ADD(
									SOLLYA_SUB(
										SOLLYA_CONST(1), // 1
										SOLLYA_MUL(SOLLYA_DIV(SOLLYA_CONST(1), SOLLYA_CONST(6)), SOLLYA_X_)),  // - 1/6*x
									SOLLYA_MUL(SOLLYA_DIV(SOLLYA_CONST(1), SOLLYA_CONST(120)), SOLLYA_POW(SOLLYA_X_, SOLLYA_CONST(2)))),  // + 1/120*x^2
								SOLLYA_MUL(SOLLYA_DIV(SOLLYA_CONST(1), SOLLYA_CONST(5040)), SOLLYA_POW(SOLLYA_X_, SOLLYA_CONST(3)))),  // - 1/5040*x^3
							SOLLYA_MUL(SOLLYA_DIV(SOLLYA_CONST(1), SOLLYA_CONST(362880)), SOLLYA_POW(SOLLYA_X_, SOLLYA_CONST(4)))),  // + 1/362880*x^4
						SOLLYA_MUL(SOLLYA_DIV(SOLLYA_CONST(1), SOLLYA_CONST(39916800)), SOLLYA_POW(SOLLYA_X_, SOLLYA_CONST(5)))),  // - 1/39916800*x^5
					SOLLYA_MUL(SOLLYA_DIV(SOLLYA_CONST(1), SOLLYA_CONST(6227020800)), SOLLYA_POW(SOLLYA_X_, SOLLYA_CONST(6)))),  // + 1/6227020800*x^6
				SOLLYA_MUL(SOLLYA_DIV(SOLLYA_CONST(1), SOLLYA_CONST(1307674368000)), SOLLYA_POW(SOLLYA_X_, SOLLYA_CONST(7)))),  // - 1/1307674368000*x^7
			SOLLYA_MUL(SOLLYA_DIV(SOLLYA_CONST(1), SOLLYA_CONST(355687428096000)), SOLLYA_POW(SOLLYA_X_, SOLLYA_CONST(8))));  // + 1/355687428096000*x^8
	// @formatter:on (Reenable reformate for CLion.)

	// Order of the polynomial approximation, i.e. how many of the terms to use.
	for (int L = 0; L < L_MAX; ++L) {
		// Polynomial terms of the sigma approximation.
		int nr_terms_sigma = L / 2 + 1;
		sollya_obj_t terms_sigma = sollya_lib_build_list(NULL);
		for (int i = 0; i < nr_terms_sigma; i++)
			terms_sigma = sollya_lib_append(terms_sigma, SOLLYA_CONST_UI64(i));

		// Polynomial terms of the tau approximation.
		int nr_terms_tau = div_floor(L - 1, 2) + 1;
		sollya_obj_t terms_tau = sollya_lib_build_list(NULL);
		for (int i = 0; i < nr_terms_tau; i++)
			terms_tau = sollya_lib_append(terms_tau, SOLLYA_CONST_UI64(i));

		// Run Remez exchange algorithm.
		str_vec coefs_sigma = remez(nr_terms_sigma, terms_sigma, sigma, range);
		str_vec coefs_tau = remez(nr_terms_tau, terms_tau, tau, range);

		// Print sigma approximation for L.
		fprintf(fp_sig, sigma_signature, L);
		std::vector<int> xs_sigma;
		std::string poly_sigma = print_poly(coefs_sigma, "d", xs_sigma);
		if (!xs_sigma.empty())
			fprintf(fp_sig, "%s", print_exp(xs_sigma, "d").c_str());
		fprintf(fp_sig, " return %s; }\n", poly_sigma.c_str());

		// Print tau approximation for L.
		fprintf(fp_tau, tau_signature, L);
		// Special handling for "no terms" avoiding ugly w * (0.0);
		if (nr_terms_tau == 0)
			fprintf(fp_tau, " return 0.0; }\n");
		else {
			std::vector<int> xs_tau;
			std::string poly_tau = print_poly(coefs_tau, "d", xs_tau);
			if (!xs_tau.empty())
				fprintf(fp_tau, "%s", print_exp(xs_tau, "d").c_str());
			if (nr_terms_tau > 1)
				fprintf(fp_tau, " return w * ( %s ); }\n", poly_tau.c_str());
			else
				fprintf(fp_tau, " return w * %s; }\n", poly_tau.c_str());
		}

		// Clean up loop-local objects.
		sollya_lib_clear_obj(terms_sigma);
		sollya_lib_clear_obj(terms_tau);
	}
	// Print postamble and close files.
	fprintf(fp_sig, "%s", sigma_postamble);
	fclose(fp_sig);
	fprintf(fp_tau, "%s", tau_postamble);
	fclose(fp_tau);

	// Clean up Sollya object and close library.
	sollya_lib_clear_obj(sigma);
	sollya_lib_clear_obj(tau);
	sollya_lib_clear_obj(range);
	sollya_lib_clear_obj(low);
	sollya_lib_clear_obj(high);
	sollya_lib_clear_obj(sollya_prec);
	sollya_lib_close();
}