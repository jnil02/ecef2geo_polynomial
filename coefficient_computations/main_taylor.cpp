// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the BSD 2-Clause License.
/*
 * Computes and generates code for minimax multiplicative correction for ECEF
 * to geodetic coordinate transformations. Minimax approximations are computed
 * with Sollya and exported to header files.
 */

#include "settings.hpp"
#include "util_eval.hpp"

#include <vector>  // std::vector.
#include <ostream>  // std::endl.
#include <string>  // std::basic_string
#include <mpreal.h>  // mpfr::mpreal type.

using mpfr::mpreal;

// Just to make the code more compact.
using str_vec = std::vector<std::string>;

// Bits of precision for Sollya amd mpfr computations.
constexpr mp_prec_t prec = 200;

static const char *sigma_preamble =
		"/*\n"
		" * Generated. DO NOT EDIT.\n"
		" */\n\n"
		"#ifndef ECEF2GEO_SIGMA2_HPP\n"
		"#define ECEF2GEO_SIGMA2_HPP\n\n"
		"namespace pteseries {\n\n";
static const char *sigma_proto =
		"namespace priv {\n\n"
		"/** Point-to-ellipse Fourier series expansions.\n"
		" *\n"
		" * Nilsson, J.-O. Point-to-ellipse Fourier series. doi:\n"
		" * https://doi.org/10.48550/arXiv.2507.08807\n"
		" *\n"
		" * @paramt L Total degree of the multiplicative polynomial correction.\n"
		" * @param d (phi-phi_c)^2.\n"
		" * @return Minimax approximation of sigma.\n"
		" */\n"
		"template<int L> inline double sigma(double d) = delete;  // Only allow provided specializations.\n";
static const char *sigma_signature = "template<> inline double sigma<%u>(double d) {";
static const char *sigma_postamble =
		"\n}  // namespace pteseries\n"
		"}  // namespace priv\n"
		"\n#endif // ECEF2GEO_SIGMA2_HPP\n";

static const char *tau_preamble =
		"/*\n"
		" * Generated. DO NOT EDIT.\n"
		" */\n\n"
		"#ifndef ECEF2GEO_TAU2_HPP\n"
		"#define ECEF2GEO_TAU2_HPP\n\n"
		"namespace pteseries {\n\n";
static const char *tau_proto =
		"namespace priv {\n\n"
		"/** Point-to-ellipse Fourier series expansions.\n"
		" *\n"
		" * Nilsson, J.-O. Point-to-ellipse Fourier series. doi:\n"
		" * https://doi.org/10.48550/arXiv.2507.08807\n"
		" *\n"
		" * @paramt L Total degree of the multiplicative polynomial correction.\n"
		" * @param d (phi-phi_c)^2.\n"
		" * @return Minimax approximation of tau.\n"
		" */\n"
		"template<int L> inline double tau(double w, double d) = delete;  // Only allow provided specializations.\n";
static const char *tau_signature = "template<> inline double tau<%u>(double w, double d) {";
static const char *tau_postamble =
		"\n}  // namespace pteseries\n"
		"}  // namespace priv\n"
		"\n#endif // ECEF2GEO_TAU2_HPP\n";

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

/** (-1)^n
 *
 * See
 * https://stackoverflow.com/questions/29110752/what-is-the-correct-way-to-obtain-1n
 *
 * @param n
 * @return
 */
inline long powm1(long n) {
	return 1L - ((n & 1L) << 1);
}

inline std::string to_scientific_string(const mpreal &x, int digits) {
	std::ostringstream ss;
	ss.precision(digits);
	ss << std::scientific << x;
	return ss.str();
}


int main() {
	mpfr::mpreal::set_default_prec(prec);

	int sig_digits = 17;  // Significant digits of double.

	// Open files and write preamble, ranges and protos for the approximations.
	FILE *fp_sig = fopen((CODEGEN_FOLDER + "ptepoly_sigma.hpp").c_str(), "w");
	fprintf(fp_sig, "%s", sigma_preamble);
	fprintf(fp_sig, "%s", sigma_proto);
	FILE *fp_tau = fopen((CODEGEN_FOLDER + "ptepoly_tau.hpp").c_str(), "w");
	fprintf(fp_tau, "%s", tau_preamble);
	fprintf(fp_tau, "%s", tau_proto);

	// Order of the polynomial approximation, i.e. how many of the terms to use.
	for (int L = 0; L < L_MAX; ++L) {
		// Polynomial terms of the sigma approximation.
		int nr_terms_sigma = L / 2 + 1;
		int nr_terms_tau = div_floor(L - 1, 2) + 1;

		// Generate coefficients.
		std::vector<std::string> coefs_sigma, coefs_tau;
		for (int l = 0; l < nr_terms_sigma; ++l) {
			mpreal tmp = mpfr::mpreal(powm1(l)) / mpfr::fac_ui(2 * l);
			coefs_sigma.emplace_back(to_scientific_string(tmp, sig_digits));
		}
		for (int l = 0; l < nr_terms_tau; ++l) {
			mpreal tmp = mpfr::mpreal(powm1(l)) / mpfr::fac_ui(2 * l + 1);
			coefs_tau.emplace_back(to_scientific_string(tmp, sig_digits));
		}

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
	}
	// Print postamble and close files.
	fprintf(fp_sig, "%s", sigma_postamble);
	fclose(fp_sig);
	fprintf(fp_tau, "%s", tau_postamble);
	fclose(fp_tau);
}