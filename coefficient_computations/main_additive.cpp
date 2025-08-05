// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the GPLv3 License.

/*
 * Computes and generates code for minimax additive corrections for ECEF to
 * geodetic coordinate transformations. Minimax approximations are computed with
 * a basic multiprecision Remez algorithm implementation and exported to header
 * files.
 */

#include "settings.hpp"
#include "util_mpgeo.hpp"
#include "util_eval.hpp"
#include "remez_mp.hpp"

#include <cstdio> // std::sprintf
#include <iostream>  // std::cout.
#include <functional>  // std::function.
#include <string>  // std::string, std::to_string.
#include <utility>  // std::move.
#include <algorithm>  // For std::function etc.
#include <vector>  // std::vector, std::begin, std::end.
#include <ostream>  // std::ostringstream, std::endl.
#include <iomanip>  // std::setprecision.
#include <stdexcept>  // std::invalid_argument.
#include <execution>  // For std::execution::par, i.e. parallel for loops.
#include <chrono>  // For measuring elapsed time.
#include <ginac/ginac.h>  // For final polynomial manipulation.

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

static const char *omega_preamble =
		"/*\n"
		" * Generated. DO NOT EDIT!\n"
		" */\n\n"
		"#ifndef ECEF2GEO_OMEGA_HPP\n"
		"#define ECEF2GEO_OMEGA_HPP\n\n"
		"namespace ecef2geo {\n\n";
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
static const char *omega_signature = "template<> inline double omega<%u,%u>(double u, double v) {";
static const char *omega_postamble =
		"\n}  // namespace\n"
		"}  // namespace priv\n"
		"\n#endif // ECEF2GEO_OMEGA_HPP\n";


static const char *mu_preamble =
		"/*\n"
		" * Generated. DO NOT EDIT!\n"
		" */\n\n"
		"#ifndef ECEF2GEO_MU_HPP\n"
		"#define ECEF2GEO_MU_HPP\n\n"
		"namespace ecef2geo {\n\n";
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
static const char *mu_signature = "template<> inline double mu<%u,%u>(double u, double v) {";
static const char *mu_postamble =
		"\n}  // namespace\n"
		"}  // namespace priv\n"
		"\n#endif // ECEF2GEO_MU_HPP\n";

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

/** Construct GiNaC symbolic expression for mu.
 *
 * @param N Summation index limit N.
 * @param M Summation index limit M.
 * @param u Polynomial variable u.
 * @param v Polynomial variable v.
 * @param syms Symbols to use for the coefficients arranged in a nested array.
 * @return A GiNaC expression for mu given the input parameters.
 */
GiNaC::ex mu_ex(int N, int M, GiNaC::symbol &u, GiNaC::symbol &v,
				GiNaC::symbol syms[N_MAX + 1][M_MAX + 1]) {
	GiNaC::ex mu_NM = u;
	for (int m = 0; m <= M; ++m)
		mu_NM += GiNaC::numeric(1, 2) * syms[0][m] * GiNaC::pow(u, m);
	for (int n = 1; n <= N; ++n) {
		GiNaC::ex sum;
		for (int m = 0; m <= M; ++m)
			sum += syms[n][m] * GiNaC::pow(u, m);
		for (int k = 0; k <= n; ++k)
			mu_NM += GiNaC::pow(GiNaC::numeric(-1), k) *
					 GiNaC::binomial(2 * n, 2 * k) * GiNaC::pow(v, k) *
					 GiNaC::pow(GiNaC::numeric(1) - v, n - k) * sum;
	}
	return mu_NM;
}

/** Construct GiNaC symbolic expression for omega.
 *
 * @param N Summation index limit N.
 * @param M Summation index limit M.
 * @param u Polynomial variable u.
 * @param v Polynomial variable v.
 * @param syms Symbols to use for the coefficients arranged in a nested array.
 * @return A GiNaC expression for omega given the input parameters.
 */
GiNaC::ex omega_ex(int N, int M, GiNaC::symbol &u, GiNaC::symbol &v,
				   GiNaC::symbol syms[N_MAX + 1][M_MAX + 1]) {
	GiNaC::ex omega_NM;
	for (int n = 1; n <= N; ++n) {
		GiNaC::ex sum;
		for (int m = 0; m <= M; ++m)
			sum += syms[n][m] * GiNaC::pow(u, m);
		for (int k = 0; k < n; ++k)
			omega_NM += GiNaC::pow(GiNaC::numeric(-1), k) *
						GiNaC::binomial(2 * n, 2 * k + 1) * GiNaC::pow(v, k) *
						GiNaC::pow(GiNaC::numeric(1) - v, n - k - 1) * sum;
	}
	return omega_NM;
}

/** Generate a C-code evaluation schema from 2-var polynomial.
 *
 * @param expr The polynomial expression.
 * @param u Inner polynomial variable of expression.
 * @param v Outer polynomial variable of expression.
 * @param N Order of N.
 * @param M Order of U.
 * @param coefs Coefficient values.
 * @param syms Coefficient symbols.
 * @param n_min first n-index.
 * @return A C-code string representation for evaluating expr.
 */
std::string gen_2var_poly(const GiNaC::ex &expr, const GiNaC::symbol &u,
						  const GiNaC::symbol &v, int N, int M,
						  coef_struct coefs[M_MAX + 1][N_MAX + 1],
						  GiNaC::symbol syms[N_MAX + 1][M_MAX + 1],
						  int n_min = 0) {
	// Retrieve "v-coefficients".
	GiNaC::lst v_coefs;
	for (int i = expr.ldegree(v); i <= expr.degree(v); ++i)
		v_coefs.append(expr.coeff(v, i));
	// Loop over the "v-coefficients".
	std::vector<std::string> v_coefs_str;
	std::vector<int> us;
	for (int i = 0; i < v_coefs.nops(); ++i) {
		// Retrieve "u-coefficients".
		std::vector<std::string> u_coefs_str;
		for (int j = v_coefs[i].ldegree(u); j <= v_coefs[i].degree(u); ++j) {
			GiNaC::ex u_coef = v_coefs[i].coeff(u, j);
			// Substituting all coefficients with values.
			for (int n = n_min; n <= N; ++n)
				for (int m = 0; m <= M; ++m) {
					// to_string gives full precision.
					GiNaC::numeric b_nm = GiNaC::numeric(
							coefs[M][n].coefs[m].to_string().c_str());
					u_coef = u_coef.subs(syms[n][m] == b_nm);
				}
			std::ostringstream s;
			s << std::setprecision(17);
			if (GiNaC::is_a<GiNaC::numeric>(u_coef))
				// TODO(JO) It would be nice to use std::format but support is
				//  so far limited.
				//  https://stackoverflow.com/questions/554063/how-do-i-print-a-double-value-with-full-precision-using-cout
				// s << std::format("{}", GiNaC::ex_to<GiNaC::numeric>(omega_u_coef).to_double());
				s << GiNaC::ex_to<GiNaC::numeric>(u_coef).to_double();
			else {
				// To see the expressions resulting in the exception below,
				// uncomment below and comment out exception.
//                s << u_coef;
				throw std::invalid_argument("Non-numeric coefficient.");
			}
			u_coefs_str.emplace_back(s.str());
		}
		std::string inner_coef = print_poly(u_coefs_str, "u", us);
		v_coefs_str.emplace_back(inner_coef);
	}
	std::vector<int> vs;
	std::string poly_str = print_poly(v_coefs_str, "v", vs);
	if (poly_str == "0")
		poly_str = "0.0";
	return print_exp(us, "u") + print_exp(vs, "v") + " return " + poly_str;
}

/** Construct the actual C-function for a 2-var polynomial.
 *
 * @param N Summation index limit N.
 * @param M Summation index limit M.
 * @param expr_str Expression string for the polynomial.
 * @param signature C-function signature.
 * @return A complete C-code function definition for the polynomial.
 */
std::ostringstream
expr2cfunc(int N, int M, const std::string &expr_str, const char *signature) {
	std::ostringstream func;
	// TODO(JO) Using std::format would be preferable but support is so far
	//  limited.
	//  https://stackoverflow.com/questions/554063/how-do-i-print-a-double-value-with-full-precision-using-cout
	char signature_str[512];  // 512 characters should always be enough.
	std::sprintf(signature_str, signature, N, M);
	func << std::string(signature_str);
	func << expr_str << "; }" << std::endl;
	return func;
}


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
	GiNaC::symbol syms[N_MAX + 1][M_MAX + 1];
	for (int n = 0; n <= N_MAX; ++n)
		for (int m = 0; m <= M_MAX; ++m)
			syms[n][m] = GiNaC::symbol(
					"c" + std::to_string(n) + std::to_string(m));

	// Open files to which generated approximations/code is to be written.
	FILE *fp_omega = fopen((CODEGEN_FOLDER + "omega.hpp").c_str(), "w");
	FILE *fp_mu = fopen((CODEGEN_FOLDER + "mu.hpp").c_str(), "w");

	// Write some defines to the generated files used for compile-time asserts.
	fprintf(fp_omega, "%s", omega_preamble);
	fprintf(fp_omega, "// Minimax approximation ranges.\n");
	fprintf(fp_omega, "constexpr double OMEGA_H_MIN = %.17e;\n", ALT_LO_LIMIT);
	fprintf(fp_omega, "constexpr double OMEGA_H_MAX = %.17e;\n", ALT_HI_LIMIT);
	fprintf(fp_omega, "constexpr double OMEGA_H_0 = %.17e;\n", ALT_REFERENCE);
	if (ALT_REFERENCE == 0.0) {
		fprintf(fp_omega, "// For ensuring compile-time removal of h_0.\n");
		fprintf(fp_omega, "#define H_0_IS_ZERO\n");
	}
	fprintf(fp_mu, "%s", mu_preamble);
	fprintf(fp_mu, "// Minimax approximation ranges.\n");
	fprintf(fp_mu, "constexpr double MU_H_MIN = %.17e;\n", ALT_LO_LIMIT);
	fprintf(fp_mu, "constexpr double MU_H_MAX = %.17e;\n", ALT_HI_LIMIT);
	fprintf(fp_mu, "constexpr double MU_H_0 = %.17e;\n", ALT_REFERENCE);

	// Write documentation and template function signatures.
	fprintf(fp_omega, "\n%s", omega_proto);
	fprintf(fp_mu, "\n%s", mu_proto);

	// Write the actual generated approximations.
	GiNaC::symbol u("u"), v("v");
	for (int N = 0; N <= N_MAX; ++N)
		for (int M = 0; M <= M_MAX; ++M) {
			GiNaC::ex om = omega_ex(N, M, u, v, syms).expand().eval();
			std::string expr = gen_2var_poly(om, u, v, N, M, coefs_b, syms, 1);
			ostrstream om_func = expr2cfunc(N, M, expr, omega_signature);
			fprintf(fp_omega, "%s", om_func.str().c_str());
			GiNaC::ex mu = mu_ex(N, M, u, v, syms).expand().eval();
			std::string mu_str = gen_2var_poly(mu, u, v, N, M, coefs_c, syms);
			ostrstream mu_func = expr2cfunc(N, M, mu_str, mu_signature);
			fprintf(fp_mu, "%s", mu_func.str().c_str());
		}

	// Clean up.
	fprintf(fp_omega, "%s", omega_postamble);
	fclose(fp_omega);
	fprintf(fp_mu, "%s", mu_postamble);
	fclose(fp_mu);
}