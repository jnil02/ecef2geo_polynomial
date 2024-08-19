// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the GPLv3 License.

/*
 * Computes and generates code for minimax multiplicative correction for ECEF
 * to geodetic coordinate transformations. Minimax approximations are computed
 * with Sollya and exported to header files.
 */

#include "settings.hpp"
#include "util_eval.hpp"
#include "remez_sollya.hpp"

#define _USE_MATH_DEFINES

#include <cmath>  // std::sqrt, std::cbrt, std::atan, std::sin, std::cos.
#include <vector>  // std::vector.
#include <iostream>  // std::cout.
#include <iomanip> // std::setprecision.
#include <limits>  // std::numeric_limits.
#include <ostream>  // std::endl.

// Just to make the code more compact.
using str_vec = std::vector<std::string>;

// Bits of precision for Sollya computations.
static constexpr int BITS_OF_PRECISION = 200;

// Samples for delta_max computations.
const int ALT_SAMPLES = 10000;
const int PHI_C_SAMPLES = 10000;

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

// Primary WGS84 constants.
constexpr double a = 6378137.0;  // Semi-major axis / equatorial radius.
constexpr double f = 1.0 / 298.257223563;  // Flattening.
// Derived WGS84 constants.
constexpr double b = a - f * a;  // Semi-minor axis / polar radius.

struct xyz {
	double x;
	double y;
	double z;
};

/** Transformation of ECEF to geodetic latitude.
 *
 * Done with Vermeille's 2004 exact method.
 *
 * @param ecef ECEF coordinate.
 * @return Geodetic latitude.
 */
static double ecef2lat(xyz ecef) {
	constexpr double a2 = a * a;
	constexpr double e2 = 1.0 - b * b / (a * a);
	constexpr double ec2 = 1.0 - e2;
	constexpr double e4 = e2 * e2;
	constexpr double ec2_a2 = ec2 / a2;
	constexpr double a2_inv = 1.0 / a2;
	constexpr double c1 = 1.0 / 6.0;
	constexpr double c2 = e4 / 4.0;
	constexpr double c3 = e2 / 2.0;

	double x2y2 = ecef.x * ecef.x + ecef.y * ecef.y;
	double sqrtx2y2 = std::sqrt(x2y2);
	double z2 = ecef.z * ecef.z;

	double q = ec2_a2 * z2;
	double p = x2y2 * a2_inv;
	double r = (p + q - e4) * c1;
	double s = c2 * p * q / (r * r * r);
	double t = std::cbrt(1.0 + s + std::sqrt(s * (2.0 + s)));
	double u = r * (1.0 + t + 1.0 / t);
	double v = std::sqrt(u * u + e4 * q);
	double w = c3 * (u + v - q) / v;
	double k = std::sqrt(u + v + w * w) - w;
	double k_k_e2 = k / (k + e2);
	double d = k_k_e2 * sqrtx2y2;

	double tmp = std::sqrt(d * d + z2);

	return 2.0 * std::atan(ecef.z / (d + tmp));
}

// TODO(JO) This takes like 10s to compute. However, it should be doable in a
//  much shorter time.  Also, the current accuracy is pretty poor. The
//  function is smooth with a single maxima so some derivative free search
//  algorithm should probably be much better.

/** Compute maximum value of "delta".
 *
 * This is the range for multiplicative corrections.
 * Computes maximum value by brute force sampling over phi_c and the specified
 * altitude range.
 *
 * @param p_min Minimum spherical altitude.
 * @param p_max Maximum spherical altitude.
 * @return The maximum value of delta within the altitude range.
 */
static double compute_max_delta(double p_min, double p_max) {
	auto ps = std::vector<double>();
	for (int i = 0; i < ALT_SAMPLES - 1; ++i)
		ps.emplace_back(p_min + (p_max - p_min) * i / ALT_SAMPLES);
	ps.emplace_back(p_max);
	auto phi_cs = std::vector<double>();
	for (int i = 0; i < PHI_C_SAMPLES - 1; ++i)
		phi_cs.emplace_back(M_PI_2 * i / PHI_C_SAMPLES);
	phi_cs.emplace_back(M_PI_2);

	double max_delta = -std::numeric_limits<double>::infinity();
	for (auto p: ps) {
		for (auto phi_c: phi_cs) {
			double sin_phi_c = std::sin(phi_c);
			double cos_phi_c = std::cos(phi_c);
			xyz ecef = {cos_phi_c * p, 0.0, sin_phi_c * p};
			double phi = ecef2lat(ecef);
			double diff = phi - phi_c;
			double de = diff * diff;
			if (max_delta < de)
				max_delta = de;
		}
	}
	return max_delta;
}

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
	const double lo = ALT_LO_LIMIT + b;
	const double hi = ALT_HI_LIMIT + a;
	const double DELTA_MIN = 0.0;
	double DELTA_MAX = compute_max_delta(lo, hi);
	std::cout << std::setprecision(17) << DELTA_MAX << std::endl;

	// Initialize Sollya.
	sollya_lib_init();
	sollya_obj_t prec;
	prec = SOLLYA_CONST_SI64(BITS_OF_PRECISION);
	sollya_lib_set_prec(prec);

	// Open files and write preamble, ranges and protos for the approximations.
	FILE *fp_sig = fopen((CODEGEN_FOLDER + "sigma.hpp").c_str(), "w");
	fprintf(fp_sig, "%s", sigma_preamble);
	fprintf(fp_sig, "// Minimax approximation ranges.\n");
	fprintf(fp_sig, "constexpr double SIGMA_H_MIN = %.17e;\n", ALT_LO_LIMIT);
	fprintf(fp_sig, "constexpr double SIGMA_H_MAX = %.17e;\n", ALT_HI_LIMIT);
	fprintf(fp_sig, "constexpr double SIGMA_DELTA_MIN = %.17e;\n", DELTA_MIN);
	fprintf(fp_sig, "constexpr double SIGMA_DELTA_MAX = %.17e;\n\n", DELTA_MAX);
	fprintf(fp_sig, "%s", sigma_proto);
	FILE *fp_tau = fopen((CODEGEN_FOLDER + "tau.hpp").c_str(), "w");
	fprintf(fp_tau, "%s", tau_preamble);
	fprintf(fp_tau, "// Minimax approximation ranges.\n");
	fprintf(fp_tau, "constexpr double TAU_H_MIN = %.17e;\n", ALT_LO_LIMIT);
	fprintf(fp_tau, "constexpr double TAU_H_MAX = %.17e;\n", ALT_HI_LIMIT);
	fprintf(fp_tau, "constexpr double TAU_DELTA_MIN = %.17e;\n", DELTA_MIN);
	fprintf(fp_tau, "constexpr double TAU_DELTA_MAX = %.17e;\n\n", DELTA_MAX);
	fprintf(fp_tau, "%s", tau_proto);

	// Minimax approximation arguments.
	sollya_obj_t range, low, high, sigma, tau;
	low = SOLLYA_CONST_SI64(DELTA_MIN);
	high = SOLLYA_CONST(DELTA_MAX);
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
	sollya_lib_clear_obj(prec);
	sollya_lib_close();
}