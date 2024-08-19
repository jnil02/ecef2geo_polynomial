// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the GPLv3 License.

/*
 * Computes and generates code for minimax inverse trigonometric function
 * approximations used for the ECEF to geodetic coordinate transformations.
 * Minimax approximations are computed with Sollya and exported to header files.
 */

#include "settings.hpp"
#include "util_eval.hpp"
#include "remez_sollya.hpp"

#include <string>  // std::string.
#include <vector>  // std::vector.

using str_vec = std::vector<std::string>;

// Bits of precision for Sollya computations.
static constexpr int BITS_OF_PRECISION = 100;
static constexpr int MAX_POLYNOMIAL_ORDER = 21;

static const char *preamble_chi =
		"/*\n"
		" * Generated. DO NOT EDIT.\n"
		" */\n\n"
		"#ifndef ECEF2GEO_CHI_HPP\n"
		"#define ECEF2GEO_CHI_HPP\n\n"
		"#define _USE_MATH_DEFINES\n"
		"#include <cmath>\n\n"
		"namespace ecef2geo {\n"
		"namespace priv {\n\n"
		"/** Minimax asin approximations valid on [0,1].\n"
		" *\n"
		" * @param N Degree of the polynomial.\n"
		" * @param x value in [0,1].\n"
		" * @return A minimax approximation of (pi/2 - asin(x))/sqrt(1.0 - x).\n"
		" */\n"
		"template<int N> inline double chi(double x) = delete;  // Only allow provided specializations.\n";
static const char *chi_signature = "template<> inline double chi<%u>(double x) {";
static const char *postamble_chi =
		"\n}  // namespace ecef2geo\n"
		"}  // namespace priv\n"
		"\n#endif // ECEF2GEO_CHI_HPP\n";

static const char *preamble_xi =
		"/*\n"
		" * Generated. DO NOT EDIT.\n"
		" */\n\n"
		"#ifndef ECEF2GEO_XI_HPP\n"
		"#define ECEF2GEO_XI_HPP\n\n"
		"namespace ecef2geo {\n"
		"namespace priv {\n\n"
		"/** Minimax atan approximations valid on [-1,1].\n"
		" *\n"
		" * @param N Degree of the polynomial.\n"
		" * @param x value in [-1,1].\n"
		" * @return An approximation of atan(x).\n"
		" */\n"
		"template<int N> inline double xi(double x) = delete;  // Only allow provided specializations.\n";
static const char *xi_signature = "template<> inline double xi<%u>(double x) {";
static const char *postamble_xi =
		"\n}  // namespace ecef2geo\n"
		"}  // namespace priv\n"
		"\n#endif // ECEF2GEO_XI_HPP\n";

int main() {
	// Initialize Sollya and set working precision.
	sollya_lib_init();
	sollya_obj_t prec;
	prec = SOLLYA_CONST_SI64(BITS_OF_PRECISION);
	sollya_lib_set_prec(prec);

	// Setup files for writing approximations to.
	FILE *fp_asin = fopen((CODEGEN_FOLDER + "chi.hpp").c_str(), "w");
	FILE *fp_atan = fopen((CODEGEN_FOLDER + "xi.hpp").c_str(), "w");
	fprintf(fp_asin, "%s", preamble_chi);
	fprintf(fp_atan, "%s", preamble_xi);

	// Remez minimax approximations arguments.
	sollya_obj_t f_asin, terms_asin, low_asin, high_asin, range_asin;
	sollya_obj_t f_atan, terms_atan, low_atan, high_atan, range_atan;
	// f=(Pi/2 - asin(x))/sqrt(1-x)
	f_asin = SOLLYA_DIV(
			SOLLYA_SUB(
					SOLLYA_DIV(SOLLYA_PI, SOLLYA_CONST_UI64(2)),
					SOLLYA_ASIN(SOLLYA_X_)),
			SOLLYA_SQRT(
					SOLLYA_SUB(SOLLYA_CONST_UI64(1), SOLLYA_X_))
	);
	low_asin = SOLLYA_CONST_SI64(0);
	high_asin = SOLLYA_CONST(0.9999999999999);
	range_asin = sollya_lib_range(low_asin, high_asin);
	f_atan = SOLLYA_ATAN(SOLLYA_X_);
	low_atan = SOLLYA_CONST_SI64(-1);
	high_atan = SOLLYA_CONST_SI64(1);
	range_atan = sollya_lib_range(low_atan, high_atan);

	// Order of the polynomial approximation, i.e. how many of the terms to use.
	for (int I = 1; I < MAX_POLYNOMIAL_ORDER; ++I) {

		terms_asin = sollya_lib_build_list(NULL);
		// Oth term needed to give a non-degenerate system. Should be ignored in
		// the results.
		terms_atan = sollya_lib_build_list(SOLLYA_CONST_UI64(0), NULL);

		for (int i = 0; i <= I; i++)
			terms_asin = sollya_lib_append(terms_asin, SOLLYA_CONST_UI64(i));
		int terms2_atan[I];
		for (int i = 0; i < I; i++) {
			terms2_atan[i] = i * 2 + 1;  // Separate list of coefs to retrieve.
			terms_atan = sollya_lib_append(terms_atan,
										   SOLLYA_CONST_UI64(i * 2 + 1));
		}

		str_vec cs_asin = remez(I + 1, terms_asin, f_asin, range_asin);
		str_vec cs_atan = remez(I, terms_atan, f_atan, range_atan, terms2_atan);

		fprintf(fp_asin, chi_signature, I);
		std::vector<int> xs_chi;
		std::string poly_chi = print_poly(cs_asin, "x", xs_chi);
		if (!xs_chi.empty())
			fprintf(fp_asin, "%s", print_exp(xs_chi, "x").c_str());
		fprintf(fp_asin, " return %s;", poly_chi.c_str());
		fprintf(fp_asin, " }\n");

		fprintf(fp_atan, xi_signature, I);
		// The polynomial is in "x2" so we first have to compute it.
		if (I > 1)
			fprintf(fp_atan, " double x2 = x * x;");
		std::vector<int> xs_xi;
		std::string poly_xi = print_poly(cs_atan, "x2", xs_xi);
		if (!xs_xi.empty())
			fprintf(fp_atan, "%s", print_exp(xs_xi, "x2").c_str());
		fprintf(fp_atan, " return x * (%s);", poly_xi.c_str());
		fprintf(fp_atan, " }\n");

		sollya_lib_clear_obj(terms_asin);
		sollya_lib_clear_obj(terms_atan);
	}
	fprintf(fp_asin, "%s", postamble_chi);
	fprintf(fp_atan, "%s", postamble_xi);
	fclose(fp_asin);
	fclose(fp_atan);

	sollya_lib_clear_obj(f_asin);
	sollya_lib_clear_obj(range_asin);
	sollya_lib_clear_obj(low_asin);
	sollya_lib_clear_obj(high_asin);
	sollya_lib_clear_obj(f_atan);
	sollya_lib_clear_obj(range_atan);
	sollya_lib_clear_obj(low_atan);
	sollya_lib_clear_obj(high_atan);
	sollya_lib_clear_obj(prec);
	sollya_lib_close();
}