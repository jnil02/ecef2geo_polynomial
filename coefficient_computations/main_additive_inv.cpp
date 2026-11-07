// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the BSD 2-Clause License.
/*
 * Computes and generates code for minimax polynomial approximations of the
 * point-to-ellipse Fourier series coefficients for ECEF to geodetic coordinate
 * transformations.
 *
 * This combines two ideas:
 *  - the point-to-ellipse Fourier series expansion of
 *    Nilsson, J.-O. Point-to-ellipse Fourier series.
 *    https://doi.org/10.48550/arXiv.2507.08807
 *  - the minimax polynomial approximation approach of
 *    Nilsson, J.-O. Minimax polynomial ECEF to geodetic coordinate
 *    transformation approximations.
 *    https://doi.org/10.48550/arXiv.2405.06572
 * i.e. the altitude/radial dependence of each Fourier series coefficient is
 * replaced by a minimax polynomial.
 *
 * Variable choice: the Fourier series of arXiv:2507.08807 expands in a / rho,
 * and the minimax polynomials of arXiv:2405.06572 are in rho. Here, however,
 * the coefficients are computed numerically, so the factor a need not be kept
 * explicit in the variable: the cleaner choice is to expand in 1 / rho and fold
 * a into the numerical coefficients (the a^i scaling applied by pot_scaling).
 * The generated functions therefore take u = 1 / rho, as used by ecef2nvaf2 in
 * ecef2geo_fourier.hpp.
 *
 * Minimax approximations are computed with Sollya and exported to header files.
 */

// Include these first since they are needed in subsequent includes.
// Included such that we can switch between mpreal.h and mplapack/mpreal.h.

// Defines to give an expected behaviour of mpfr.h for mpreal.h.
#define MPFR_USE_NO_MACRO
#define MPFR_USE_INTMAX_T
#include <cstdint>  // Required for mpreal.h. when mpfr is included first.
#include <mpfr.h>
#include <mpreal.h>

#include "util_four.hpp"
#include "util_sym.hpp"
#include <symengine/expression.h>
#include <symengine/real_mpfr.h>  // SymEngine::mpfr_class.

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
#include "util_mpgeo.hpp"  // f_c

#include "consts.hpp"


#define _USE_MATH_DEFINES
#include <cmath>  // std::sqrt, std::cbrt, std::atan, std::sin, std::cos.
#include <vector>  // std::vector.
#include <iostream>  // std::cout.
#include <limits>  // std::numeric_limits.
#include <ostream>  // std::endl.
#include <string>  // std::basic_string

#include <coefficients.hpp>


// Just to make the code more compact.
using vec_mpfr = std::vector<mpfr_ptr>;

// Bits of precision for Sollya amd mpfr computations.
constexpr mp_prec_t prec = 200;
// Mpfr consts with appropriate precision.
wgs84_mpfr_constants util_mpgeo::mpfr_consts = wgs84_mpfr_constants(prec);

static const char *phi_preamble =
		"/*\n"
		" * Generated. DO NOT EDIT.\n"
		" */\n\n"
		"#ifndef ECEF2GEO_PHI_HPP\n"
		"#define ECEF2GEO_PHI_HPP\n\n"
		"namespace ecef2geo {\n\n";
static const char *phi_proto =
		"namespace priv {\n\n"
		"/** Minimax approximations of phi.\n"
		" *\n"
		" * These approximate the point-to-ellipse Fourier series coefficients of\n"
		" *   Nilsson, J.-O. Point-to-ellipse Fourier series.\n"
		" *   https://doi.org/10.48550/arXiv.2507.08807\n"
		" * with minimax polynomials following the approach of\n"
		" *   Nilsson, J.-O. Minimax polynomial ECEF to geodetic coordinate\n"
		" *   transformation approximations.\n"
		" *   https://doi.org/10.48550/arXiv.2405.06572\n"
		" * i.e. a combination of the two.\n"
		" *\n"
		" * The former expands in a / rho and the latter uses polynomials in rho.\n"
		" * Since the coefficients here are computed numerically, the factor a\n"
		" * need not be kept explicit in the variable: the cleaner choice is to\n"
		" * expand in 1 / rho and fold a into the (numerical) coefficients. The\n"
		" * generated functions therefore take u = 1 / rho.\n"
		" *\n"
		" * @paramt N,M Orders of the polynomial approximation.\n"
		" * @param u 1 / rho.\n"
		" * @param v t^2.\n"
		" * @return Minimax approximation of phi.\n"
		" */\n"
		"template<int N, int M> inline double phi(double u, double v) = delete;  // Only allow provided specializations.\n";
static const char *phi_signature = "template<> inline double %s<%u,%u>(double u, double v) {";
static const char *phi_postamble =
		"\n}  // namespace priv\n"
		"}  // namespace ecef2geo\n"
		"\n#endif // ECEF2GEO_PHI_HPP\n";

static const char *h_preamble =
		"/*\n"
		" * Generated. DO NOT EDIT.\n"
		" */\n\n"
		"#ifndef ECEF2GEO_H_HPP\n"
		"#define ECEF2GEO_H_HPP\n\n"
		"namespace ecef2geo {\n\n";
static const char *h_proto =
		"namespace priv {\n\n"
		"/** Minimax approximations of h.\n"
		" *\n"
		" * These approximate the point-to-ellipse Fourier series coefficients of\n"
		" *   Nilsson, J.-O. Point-to-ellipse Fourier series.\n"
		" *   https://doi.org/10.48550/arXiv.2507.08807\n"
		" * with minimax polynomials following the approach of\n"
		" *   Nilsson, J.-O. Minimax polynomial ECEF to geodetic coordinate\n"
		" *   transformation approximations.\n"
		" *   https://doi.org/10.48550/arXiv.2405.06572\n"
		" * i.e. a combination of the two.\n"
		" *\n"
		" * The former expands in a / rho and the latter uses polynomials in rho.\n"
		" * Since the coefficients here are computed numerically, the factor a\n"
		" * need not be kept explicit in the variable: the cleaner choice is to\n"
		" * expand in 1 / rho and fold a into the (numerical) coefficients. The\n"
		" * generated functions therefore take u = 1 / rho.\n"
		" *\n"
		" * @paramt N,M Orders of the polynomial approximation.\n"
		" * @param u 1 / rho.\n"
		" * @param v t^2.\n"
		" * @return Minimax approximation of h.\n"
		" */\n"
		"template<int N, int M> inline double h(double u, double v) = delete;  // Only allow provided specializations.\n";
static const char *h_signature = "template<> inline double %s<%u,%u>(double u, double v) {";
static const char *h_postamble =
		"\n}  // namespace priv\n"
		"}  // namespace ecef2geo\n"
		"\n#endif // ECEF2GEO_H_HPP\n";


/* Plan för att implementera minimax
 *
 * - Lopa över funktioner som skall approximeras
 * 		- Lopa över N (Fourier koefficient)
 * 			- Generera Sollya-funktion med hjälp av lib.
 * 			- Polynom för Fourier-koefficient med tillräckligt stort antal varrho-termer.
 * 			- Lopa över M (antal varrho termer)
 * 				- Se till att Remez börjar vid rätt gradtal.
 * 				- Remez-approximation
 * 				- Generera polynom
 * 		- Stäng fil
 */

// FIXME(JO) Write some struct with
//  - File name, preamble, postamble, proto, signature
//  - Lamda for the reference.
//  - Lowest polynomial order
//  - Lowest sin-power order
//  - "offset"
//  Det är mer eller mindre identiskt med struct:en i main_fourier.cpp.
//  Därefter kan jag loopa över structen och få rätt gränser på

struct minimax_funcs {
	std::string filename;
	std::string preamble;
	std::string postamble;
	std::string proto;
	std::string signature;
	std::string name;
	int n_min;
	std::function<sollya_obj_t(int)> ref;
	mpfr::mpreal offset;
	mpfr::mpreal scale;
};


int main() {

//	mpreal VARRHO_MIN(0);
//	mpreal VARRHO_MIN(0.5);
	// The point-to-ellipse Fourier series expands in a / rho (arXiv:2507.08807)
	// and the generated polynomial is evaluated at a / rho (u = 1 / rho with a
	// folded into the coefficients). The minimax is therefore fitted over the
	// operating range of a / rho = [a / rho_max, a / rho_min], where
	//   rho_min = b + ALT_LO_LIMIT  (pole, lowest altitude)
	//   rho_max = a + ALT_HI_LIMIT  (equator, highest altitude).
	mpreal VARRHO_MIN = mpreal(util_mpgeo::mpfr_consts.a) / (mpreal(util_mpgeo::mpfr_consts.a) + ALT_HI_LIMIT);
	mpreal VARRHO_MAX = mpreal(util_mpgeo::mpfr_consts.a) / (mpreal(util_mpgeo::mpfr_consts.b) + ALT_LO_LIMIT);

	// Initialize Sollya.
	sollya_lib_init();
	sollya_obj_t sollya_prec;
	sollya_prec = SOLLYA_CONST_SI64(prec);
	sollya_lib_set_prec(sollya_prec);

	// Minimax approximation arguments.
	sollya_obj_t range, low, high;
	low = sollya_lib_parse_string(VARRHO_MIN.toString().c_str());
	high = sollya_lib_parse_string(VARRHO_MAX.toString().c_str());
	range = sollya_lib_range(low, high);

	// Generate reference functions for coefficients.
	sollya_obj_t phi_ref[N_MAX+1];
	sollya_obj_t h_ref[N_MAX+1];
	int phi_n_min = 0;
	for (int n = phi_n_min; n <= N_MAX; ++n) {
		phi_ref[n] = SOLLYA_CONST_SI64(0);
		for (int k = 1; k <= M_MAX + 4; ++k) {
			phi_ref[n] = SOLLYA_ADD(phi_ref[n],
									  SOLLYA_MUL(sollya_lib_parse_string(
														 d_phi(n, k).toString().c_str()),
												 SOLLYA_POW(SOLLYA_X_,
															SOLLYA_CONST(k))));
		}
	}
	int h_n_min = 1;
	for (int n = h_n_min; n <= N_MAX; ++n) {
		h_ref[n] = SOLLYA_CONST_SI64(0);
		for (int k = 0; k <= M_MAX + 4; ++k) {
			h_ref[n] = SOLLYA_ADD(h_ref[n],
									SOLLYA_MUL(sollya_lib_parse_string(d_h(n, k).toString().c_str()),
											   SOLLYA_POW(SOLLYA_X_, SOLLYA_CONST(k))));
		}
	}

	std::vector<minimax_funcs> funcs;
	funcs.push_back({
			.filename   = "phi.hpp",
			.preamble   = phi_preamble,
			.postamble  = phi_postamble,
			.proto      = phi_proto,
			.signature  = phi_signature,
			.name       = "phi",
			.n_min      = phi_n_min,
			.ref        = [&phi_ref](int n) { return phi_ref[n]; },
			.offset     = mpfr::mpreal(0),
			.scale      = mpfr::mpreal(1)
	});
	funcs.push_back({
			.filename   = "h.hpp",
			.preamble   = h_preamble,
			.postamble  = h_postamble,
			.proto      = h_proto,
			.signature  = h_signature,
			.name       = "h",
			.n_min      = h_n_min,
			.ref        = [&h_ref](int n) { return h_ref[n]; },
			.offset     = -mp_a(),
			.scale      = mp_a()
	});

	SymEngine::Expression u("u"), v("v");
	for(minimax_funcs &func : funcs) {

		// Open files and write preamble, ranges and protos for the approximations.
		FILE *fp = fopen((CODEGEN_FOLDER + func.filename).c_str(), "w");
		fprintf(fp, "%s", func.preamble.c_str());
		fprintf(fp, "// Minimax approximation ranges.\n");
		fprintf(fp, "//constexpr double H_MIN = %.17e;\n\n", ALT_LO_LIMIT);
		fprintf(fp, "//constexpr double E2 = %.17e;\n\n", mp_e2().toDouble());
		fprintf(fp, "//constexpr double A = %.17e;\n\n", mp_a().toDouble());
		fprintf(fp, "%s", func.proto.c_str());

		// Order of the polynomial approximation, i.e. how many of the terms to use.
		std::vector<std::vector<vec_mpfr>> minimax_all;
		for (int M = 0; M <= M_MAX; ++M) {
			std::vector<vec_mpfr> minimax_coefs{};
			// Loop over sin-power series coefficients.
			for (int n = func.n_min; n <= N_MAX; ++n) {
				// Polynomial terms of the approximation.
				// Note, the zeroth term has to be included even if it is not
				// part of the series expansions since otherwise the system
				// becomes degenerate. Further, it cannot be dropped later on
				// since we do not have a symmetric case when it should be zero.
				sollya_obj_t terms = sollya_lib_build_list(NULL);
				for (int i = 0; i <= M; i++)
					terms = sollya_lib_append(terms,
													SOLLYA_CONST_UI64(i));

				// Run Remez exchange algorithm.
				sollya_obj_t f = func.ref(n);
				vec_mpfr result = remez(M + 1, terms, f, range);

				pot_scaling(0, result, util_mpgeo::mpfr_consts.a);
				minimax_coefs.push_back(result);
				// Clean up loop-local objects.
				sollya_lib_clear_obj(terms);

			}
			minimax_all.push_back(minimax_coefs);
		}
		for (int M = 0; M <= M_MAX; ++M) {
			// Coefficient of degree M in varrho.
			for (int n = func.n_min; n <= N_MAX; ++n) {
				// With "n" coefficients sin-power coefficients.

				// Sum an expression of up to n:th sin-power terms.
				// And M:th degree varrho terms.
				SymEngine::Expression d{0};
				for (int _n = func.n_min; _n <= n; ++_n) {
					SymEngine::Expression tmp(0);
					for (int m = 0; m <= M; m++) {
						SymEngine::mpfr_class mc(minimax_all[M][_n - func.n_min][m]);
						SymEngine::Expression a = SymEngine::Expression(SymEngine::real_mpfr(mc));
						tmp = tmp + a * SymEngine::pow(u, m);
					}
					d = d + tmp * SymEngine::pow(v, _n);
				}
				SymEngine::mpfr_class mc3(func.scale.mpfr_srcptr());
				d *= SymEngine::Expression(SymEngine::real_mpfr(mc3));
				SymEngine::mpfr_class mc2(func.offset.mpfr_srcptr());
				d += SymEngine::Expression(SymEngine::real_mpfr(mc2));
				std::string expr = gen_2var_poly(d, u, v, [](auto a) { return a; });
				std::ostringstream om_func = expr2cfunc(n, M, expr, func.signature.c_str(), func.name.c_str());  // Generate a complete "C-function".
				fprintf(fp, "%s", om_func.str().c_str());
			}
		}
		for (int M = 0; M <= M_MAX; ++M) {
			// Clear Remez results of M degree.
			for (int n = func.n_min; n <= N_MAX; ++n) {
				clear_mpfrs(minimax_all[M][n - func.n_min]);
			}
		}
		// Print postamble and close files.
		fprintf(fp, "%s", func.postamble.c_str());
		fclose(fp);
	}

	// Clean up Sollya object and close library.
	for (int n = phi_n_min; n <= N_MAX; ++n)
		sollya_lib_clear_obj(phi_ref[n]);
	for (int n = h_n_min; n <= N_MAX; ++n)
		sollya_lib_clear_obj(h_ref[n]);
	sollya_lib_clear_obj(range);
	sollya_lib_clear_obj(low);
	sollya_lib_clear_obj(high);
	sollya_lib_clear_obj(sollya_prec);
	sollya_lib_close();
}