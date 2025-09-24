// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the BSD 2-Clause License.

/*
 * Computes and generates code for minimax additive corrections for ECEF to
 * geodetic coordinate transformations. Minimax approximations are computed with
 * a basic multiprecision Remez algorithm implementation and exported to header
 * files.
 */

// Include these first since they are needed in subsequent includes.
// Included such that we can switch between mpreal.h and mplapack/mpreal.h.

#include <mpfr.h>
#include <mplapack/mpreal.h>
// ***** To handle inconsistencies between mplapack/mpreal.h and mpreal.h *****
void mpfr_set_mpreal(mpfr_ptr to, mpfr::mpreal &val) {
		mpfr_ptr r_ptr(val);
		mpfr_set(to, r_ptr, MPFR_RNDN);
}
inline bool isnan    (const mpfr::mpreal& v){    return mpfr::_isnan(v);    }
inline bool isinf    (const mpfr::mpreal& v){    return mpfr::_isinf(v);    }
inline mpfr::mpreal sum(const mpfr::mpreal tab[], unsigned long int n) {
	return mpfr::sum(tab, n, MPFR_RNDN);
}
// ****************************************************************************

#include "settings.hpp"
#include "util_mpgeo.hpp"
#include "util_sym.hpp"
#include "remez_mp.hpp"

#include <cstdio> // std::sprintf
#include <iostream>  // std::cout.
#include <functional>  // std::function.
#include <string>  // std::string, std::to_string.
#include <utility>  // std::move.
#include <algorithm>  // For std::function etc.
#include <vector>  // std::vector, std::begin, std::end.
#include <ostream>  // std::ostringstream, std::endl.
#include <execution>  // For std::execution::par, i.e. parallel for loops.
#include <chrono>  // For measuring elapsed time.

#include <symengine/eval.h>
#include <symengine/real_mpfr.h>

// Just get rid of some long type names.
using time_point = std::chrono::steady_clock::time_point;
using ostrstream = std::ostringstream;

// TODO(JO) It is unclear why such a large number of bits is required.
//  Something is clearly ill conditioned.
// Number of bits of precision in the computations.
// 200 is required to get stable results for high order coefficients. Reasonable
// results can be attained with 100. 53 corresponds to double. (Multiplying with
// 0.3 gives roughly the decimal precision.)
constexpr mp_prec_t prec = 200;
// Mpfr consts with appropriate precision.
wgs84_mpfr_constants util_mpgeo::mpfr_consts = wgs84_mpfr_constants(prec);
// Max number of Remez algorithm iterations.
constexpr int ITER_MAX = 500;
// Terminate Remez algorithm when relative difference between max errors is less
// than 1 + REL_ERR_MAX.
constexpr double REL_ERR_MAX = 0.00001;
// Number of integral samples for computing Fourier coefficients.
// 2048 samples versus 64 samples gives a difference in 22:nd digit.
// 256 gives 22:nd digit.
// 128 gives 23:rd digit.
// 46 samples give difference in the 14th digit.
// 32 samples give differences in the first digit.
constexpr int NR_INTEGRAL_SAMPLES = 128;

static const char *omega_proto =
		"namespace priv {\n\n"
		"/** Polynomial additive latitude correction\n"
		" *\n"
		" * Nilsson, J.-O. Minimax polynomial ECEF to geodetic coordinate transformation\n"
		" * approximations. doi: https://doi.org/10.48550/arXiv.2405.06572\n"
		" *\n"
		" * @paramt N Number of Fourier coefficients in the approximation.\n"
		" * @paramt M Number of coefficient for the altitude dependency.\n"
		" * @param u h_c.\n"
		" * @param v t^2.\n"
		" * @return Minimax approximation of omega.\n"
		" */\n"
		"template<int N,int M> inline double omega(double u, double v) = delete;  // Only allow provided specializations.\n";
static const char *mu_proto =
		"namespace priv {\n\n"
		"/** Polynomial additive altitude correction\n"
		" *\n"
		" * Nilsson, J.-O. Minimax polynomial ECEF to geodetic coordinate transformation\n"
		" * approximations. doi: https://doi.org/10.48550/arXiv.2405.06572\n"
		" *\n"
		" * @paramt N Number of Fourier coefficients in the approximation.\n"
		" * @paramt M Number of coefficient for the altitude dependency.\n"
		" * @param u h_c.\n"
		" * @param v t^2.\n"
		" * @return Minimax approximation of mu.\n"
		" */\n"
		"template<int N,int M> inline double mu(double u, double v) = delete;  // Only allow provided specializations.\n";
static const char *signature = "template<> inline double %s<%u,%u>(double u, double v) {";
static const char *preamble =
		"/*\n"
		" * Generated. DO NOT EDIT!\n"
		" */\n\n"
		"#ifndef ECEF2GEO_%s_HPP\n"
		"#define ECEF2GEO_%s_HPP\n\n"
		"namespace ecef2geo {\n\n";
static const char *postamble =
		"\n}  // namespace\n"
		"}  // namespace priv\n"
		"\n#endif // ECEF2GEO_%s_HPP\n";

/** Struct for input parameters and results of Remez optimization.
 *
 * The parameters and the results need to be gathered in a common object such
 * that they can be placed in a standard container such that language support
 * for parallel execution over containers can be used.
 */
struct coef_struct {
	int n;  // Fourier coefficient index.
	int M;  // Number of polynomial coefficients of the approximation.
	char c;  // Coefficient "name", e.g. "b".
	std::function<mpfr::mpreal(const mpfr::mpreal &)> f;
	mpfr::mpreal *coefs;  // Pointer to the coefficients.

	// Just for initializing the array.
	coef_struct() : n(-1), M(-1), c('x'), f(nullptr), coefs(nullptr) {}

	coef_struct(int n, int M, char c,
				std::function<mpfr::mpreal(const mpfr::mpreal &)> f,
				mpfr_prec_t prec) : n(n), M(M), c(c), f(std::move(f)) {
		coefs = new mpfr::mpreal[M + 1];
		for (int i = 0; i <= M; ++i)
			coefs[i].set_prec(prec);
	}

	// Delete copy constructor.
	coef_struct(const coef_struct &coefStruct) = delete;

	// Move constructors such that memory is not freed.
	coef_struct(coef_struct &&other) noexcept {
		n = other.n;
		M = other.M;
		c = other.c;
		f = std::move(other.f);
		coefs = other.coefs;
		other.coefs = nullptr;  // Such that destructor does not free memory.
	}

	coef_struct &operator=(coef_struct &&other) noexcept {
		n = other.n;
		M = other.M;
		c = other.c;
		if (other.f)
			f = std::move(other.f);
		delete[] coefs;
		coefs = other.coefs;
		other.coefs = nullptr;  // Such that destructor does not free memory.
		return *this;
	}

	~coef_struct() {
		delete[] coefs;
	}
};


int main() {
	// TODO(JO) This should not be required. The precision of the arithmetics is
	//  specified everywhere but in mpreal.h on line 1897, in
	//  "inline mpreal::mpreal(const mpreal &u)" the precision is taken from
	//  default_prec. This seems like a bug or at least like a weird policy.
	//  This behaviour is at least in contrast to that of gmpxx.h. However, to
	//  handle this case, the default precision is set here.
	mpreal::default_prec = prec;

	// Reference radius.
	mpreal h_0 = mpreal(ALT_REFERENCE, prec);

	// Upper and lower limits for the minimax approximation.
	// Note that this is relative h_0 and not geodetic altitude.
	mpreal lo = mpreal(ALT_LO_LIMIT, prec);
	lo = lo + mpreal(util_mpgeo::mpfr_consts.b) - h_0;
	mpreal hi = mpreal(ALT_HI_LIMIT, prec);
	hi = hi + mpreal(util_mpgeo::mpfr_consts.a) - h_0;

	// Numerical differentiation difference.
	// Selecting a proper "h" is a non-trivial problem.
	// See https://en.wikipedia.org/wiki/Numerical_differentiation.
	// Selecting a too low value will make the algorithm not converge.
	mpreal h;
	h.set_prec(prec);
	h = mul_2si(hi - lo, -(prec / 4), MPFR_RNDN);
	// TODO(JO) 4 seems to be a large number. 2 does not work.
	//  The coefficient-function seems noisy. Could be related to why an
	//  excessive number of bits of precision is required to make the matrices
	//  in the Remez algorithm non-singular.

	// Fill in Remez "job" info.
	std::vector<coef_struct> coefsIO;
	for (int n = 0; n <= N_MAX; n++) {
		// Count down in order to start with the most expensive jobs.
		for (int M = M_MAX; M >= 0; M--) {
			// Zeroth latitude error coefficient is zero so don't compute.
			if (n != 0) {
				// Compute the n:th Fourier coefficient of the latitude error.
				std::function<mpreal(const mpreal &)> f = [&h_0, n](
						const mpreal &h_c) {
					auto integrand = [&h_c, &h_0, n](const mpreal &lat) {
						return (util_mpgeo::f_c(lat, h_c, h_0, prec) - lat) *
							   sin(lat * 2 * n) / const_pi(prec) * 2;
					};
					// TODO(JO) Use something better than composite Simpson's.
					//  For oscillatory functions Clenshaw-Curtis integration is
					//  preferable. See
					//  https://www.gnu.org/software/gsl/doc/html/integration.html#qawo-adaptive-integration-for-oscillatory-functions
					//  This could lead to smaller number of samples required
					//  and hence an overall speedup.
					return util_mpgeo::simpson_integration(
							integrand, const_pi(prec) / -2, const_pi(prec) / 2,
							NR_INTEGRAL_SAMPLES, prec);
				};
				coefsIO.emplace_back(n, M, 'b', f, prec);
			}
			// Compute the n:th Fourier coefficient of the altitude error.
			std::function<mpreal(const mpreal &)> f = [&h_0, n](
					const mpreal &h_c) {
				auto integrand = [&h_c, &h_0, n](const mpreal &lat) {
					return (util_mpgeo::g_c(lat, h_c, h_0, prec) - h_c) *
						   cos(lat * 2 * n) / const_pi(prec) * 2;
				};
				return util_mpgeo::simpson_integration(
						integrand, const_pi(prec) / -2, const_pi(prec) / 2,
						NR_INTEGRAL_SAMPLES, prec);
			};
			coefsIO.emplace_back(n, M, 'c', f, prec);
		}
	}

	// Just for displaying how long computations took.
	time_point begin = std::chrono::steady_clock::now();

	// Parallel execution of Remez exchange algorithm.
	// Std::execution::par gives the actual parallelism. Change to
	// std::execution::seq for debugging.
	std::for_each(
			std::execution::par,
			std::begin(coefsIO), std::end(coefsIO),
			[&h, &lo, &hi](const coef_struct &val) {
				remez(val.coefs, val.M + 1, val.f, h, lo, hi, REL_ERR_MAX,
					  ITER_MAX, prec);
			});

	// Display how long the computations took.
	time_point end = std::chrono::steady_clock::now();
	std::cout << "Coefficient computations took: "
			  << std::chrono::duration_cast<std::chrono::seconds>(
					  end - begin).count() << "[s]" << std::endl;

	// Build arrays of coefficient structs to enable access with indices.
	coef_struct coefs_c[M_MAX + 1][N_MAX + 1];
	coef_struct coefs_b[M_MAX + 1][N_MAX + 1];
	for (auto &c: coefsIO) {
		if (c.c == 'c')
			coefs_c[c.M][c.n] = std::move(c);
		else if (c.c == 'b')
			coefs_b[c.M][c.n] = std::move(c);
	}

	// Generate coefficient symbols to use in the symbolic polynomials.
	SymEngine::Expression syms[N_MAX + 1][M_MAX + 1];
	for (int n = 0; n <= N_MAX; ++n)
		for (int m = 0; m <= M_MAX; ++m)
			syms[n][m] = SymEngine::Expression(SymEngine::symbol(
					"c" + std::to_string(n) + std::to_string(m)));

	// Open files to which generated approximations/code is to be written.
	FILE *fp_omega = fopen((CODEGEN_FOLDER + "omega.hpp").c_str(), "w");
	FILE *fp_mu = fopen((CODEGEN_FOLDER + "mu.hpp").c_str(), "w");

	// Write some defines to the generated files used for compile-time asserts.
	fprintf(fp_omega, preamble, "OMEGA", "OMEGA");
	fprintf(fp_omega, "// Minimax approximation ranges.\n");
	fprintf(fp_omega, "constexpr double OMEGA_H_MIN = %.17e;\n", ALT_LO_LIMIT);
	fprintf(fp_omega, "constexpr double OMEGA_H_MAX = %.17e;\n", ALT_HI_LIMIT);
	fprintf(fp_omega, "constexpr double OMEGA_H_0 = %.17e;\n", ALT_REFERENCE);
	if (ALT_REFERENCE == 0.0) {
		fprintf(fp_omega, "// For ensuring compile-time removal of h_0.\n");
		fprintf(fp_omega, "#define H_0_IS_ZERO\n");
	}
	fprintf(fp_mu, preamble, "MU", "MU");
	fprintf(fp_mu, "// Minimax approximation ranges.\n");
	fprintf(fp_mu, "constexpr double MU_H_MIN = %.17e;\n", ALT_LO_LIMIT);
	fprintf(fp_mu, "constexpr double MU_H_MAX = %.17e;\n", ALT_HI_LIMIT);
	fprintf(fp_mu, "constexpr double MU_H_0 = %.17e;\n", ALT_REFERENCE);

	// Write documentation and template function signatures.
	fprintf(fp_omega, "\n%s", omega_proto);
	fprintf(fp_mu, "\n%s", mu_proto);

	// Write the actual generated approximations.
	SymEngine::Expression u("u"), v("v");
	for (int N = 0; N <= N_MAX; ++N)
		for (int M = 0; M <= M_MAX; ++M) {
			SymEngine::Expression om = SymEngine::expand(omega_ex(N, M, u, v, syms));
			// Lambda for substituting symbolic coefficients in expression.
			auto subs = [N, M, &coefs_b, &syms](const SymEngine::Expression& u_coef) -> SymEngine::Expression {
				std::map<SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCPBasicKeyLess> sub_map;
				for (int n = 1; n <= N; ++n)
					for (int m = 0; m <= M; ++m) {
						// Construct a SymEngine object directly from the mpreal object.
						SymEngine::mpfr_class mc(coefs_b[M][n].coefs[m].operator mpfr_ptr());
						// M, n and m above designate the m:th coefficient of the
						// M-coefficient approximation of the n:th sin/cos
						// coefficient.
						SymEngine::RCP<const SymEngine::Basic> mp_basic = SymEngine::real_mpfr(mc);
						sub_map[syms[n][m].get_basic()] = mp_basic;
					}
				return u_coef.subs(sub_map);
			};
			std::string expr = gen_2var_poly(om, u, v, subs);
			ostrstream om_func = expr2cfunc(N, M, expr, signature, "omega");
			fprintf(fp_omega, "%s", om_func.str().c_str());
			SymEngine::Expression mu = SymEngine::expand(mu_ex(N, M, u, v, syms));
			// Lambda for substituting symbolic coefficients in expression.
			auto subs2 = [N, M, &coefs_c, &syms](const SymEngine::Expression& u_coef) -> SymEngine::Expression {
				std::map<SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCPBasicKeyLess> sub_map;
				for (int n = 0; n <= N; ++n)
					for (int m = 0; m <= M; ++m) {
						// Construct a SymEngine object directly from the mpreal object.
						SymEngine::mpfr_class mc(coefs_c[M][n].coefs[m].operator mpfr_ptr());
						// M, n and m above designate the m:th coefficient of the
						// M-coefficient approximation of the n:th sin/cos
						// coefficient.
						SymEngine::RCP<const SymEngine::Basic> mp_basic = SymEngine::real_mpfr(mc);
						sub_map[syms[n][m].get_basic()] = mp_basic;
					}
				return u_coef.subs(sub_map);
			};
			std::string mu_str = gen_2var_poly(mu, u, v, subs2);
			ostrstream mu_func = expr2cfunc(N, M, mu_str, signature, "mu");
			fprintf(fp_mu, "%s", mu_func.str().c_str());
		}

	// Clean up.
	fprintf(fp_omega, postamble, "OMEGA");
	fclose(fp_omega);
	fprintf(fp_mu, postamble, "MU");
	fclose(fp_mu);
}