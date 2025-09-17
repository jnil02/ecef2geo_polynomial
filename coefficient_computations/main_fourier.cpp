// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the BSD 2-Clause License.

/*
 * Computes and generates code for truncated series expansions for ECEF to
 * geodetic coordinate transformations. Coefficients are computed with the
 * external lib ptecoef.
 */

#include "settings.hpp"
#include "util_sym.hpp"
#include "util_four.hpp"

#include <cstdio> // std::sprintf.
#include <functional>  // std::function.
#include <string>  // std::string, std::to_string.
#include <utility>  // std::move.
#include <vector>  // std::vector, std::begin, std::end.
#include <chrono>  // For measuring elapsed time.
#include <mpreal.h>  // mpfr::mpreal type.
#include <symengine/real_mpfr.h>  // SymEngine::mpfr_class.

// Just get rid of some long type names.
using time_point = std::chrono::steady_clock::time_point;
using ostrstream = std::ostringstream;

// Number of bits of precision in the computations.
constexpr mp_prec_t prec = 200;

static const char *preamble =
		"/*\n"
		" * Generated. DO NOT EDIT!\n"
		" */\n\n"
		"#ifndef POINT_TO_ELLIPSE_FOURIER_%s_HPP\n"
		"#define POINT_TO_ELLIPSE_FOURIER_%s_HPP\n\n"
		"namespace ecef2geo {\n\n";
static const char *proto =
		"namespace priv {\n\n"
		"/** Point-to-ellipse Fourier series expansions.\n"
		" *\n"
		" * Nilsson, J.-O. Point-to-ellipse Fourier series. doi:\n"
		" * https://doi.org/10.48550/arXiv.2507.08807\n"
		" *\n"
		" * @paramt N Number of Fourier coefficients in the approximation.\n"
		" * @paramt M Number of coefficient for the altitude dependency.\n"
		" * @param u a / rho.\n"
		" * @param v t^2.\n"
		" * @return Minimax approximation of omega.\n"
		" */\n"
		"template<int N,int M> inline double %s(double u, double v) = delete;  // Only allow provided specializations.\n";
static const char *signature = "template<> inline double %s<%u,%u>(double u, double v) {";
static const char *postamble =
		"\n}  // namespace\n"
		"}  // namespace priv\n"
		"\n#endif // POINT_TO_ELLIPSE_FOURIER_%s_HPP\n";

/** Struct for input parameters and results of Remez optimization.
 *
 * The parameters and the results need to be gathered in a common object such
 * that they can be placed in a standard container such that language support
 * for parallel execution over containers can be used.
 */
struct coef_structs {
	const char* name;
	int n_min;  // Lower fourier coefficient index.
	int k_min;  // Lower varrho coefficient index.
	std::function<mpfr::mpreal(int, int)> coef;  // Get coefficient.

	// Just for initializing the array.
	coef_structs() : name(nullptr), n_min(-1), k_min(-1), coef(nullptr) {}

	coef_structs(const char* name, int n, int M,
				 std::function<mpfr::mpreal(int, int)> coef)
				 : name(name), n_min(n), k_min(M), coef(std::move(coef)) {
	}

	// Delete copy constructor.
	coef_structs(const coef_structs &coefStruct) = delete;

	// Move constructors such that memory is not freed.
	coef_structs(coef_structs &&other) noexcept {
		name = other.name;
		n_min = other.n_min;
		k_min = other.k_min;
		coef = std::move(other.coef);
	}

	coef_structs &operator=(coef_structs &&other) noexcept {
		name = other.name;
		n_min = other.n_min;
		k_min = other.k_min;
		coef = std::move(other.coef);
		return *this;
	}

	~coef_structs() {
	}
};



/* PLAN
 * - Städa upp och commit:a!!!!!!!! (bara polynomgenereringen)
 * - Skapa en helt separat main som genererar omega och tau (inte minimax).
 * - Mha dessa, gör faktiska nva-approximationer.
 * - Göra alla approximationer som jag har för minimax-varianten.
 * - Gör minimax av dessa!
 * - Ta över världen.
 */

int main() {
	set_precision_bits(prec);

	// Generate coefficient symbols to use in the symbolic polynomials.
	SymEngine::Expression syms[N_MAX + 1][M_MAX + 1];
	for (int n = 0; n <= N_MAX; ++n)
		for (int m = 0; m <= M_MAX; ++m)
			syms[n][m] = SymEngine::Expression(SymEngine::symbol(
					"c" + std::to_string(n) + std::to_string(m)));

	// Definitions of the polynomials to generate.
	std::vector<coef_structs> funcs;
	funcs.emplace_back("ptepoly_phi",
					   0,  // n_min
					   1,  // k_min
					   [](int n, int k) { return d_phi(n,k); }
					   );
	funcs.emplace_back("ptepoly_h",
					   1,  // n_min
					   0,  // k_min
					   [](int n, int k) { return d_h(n,k); }
					   );
	funcs.emplace_back("ptepoly_sin",
					   0,  // n_min
					   1,  // k_min
					   [](int n, int k) { return d_sin(n,k); }
					   );
	funcs.emplace_back("ptepoly_cos",
					   1,  // n_min
					   1,  // k_min
					   [](int n, int k) { return d_cos(n,k); }
					   );

	// Note, the coefficient functions in funcs are not thread safe so we cannot
	// run the jobs in parallel.
	for (coef_structs &coef :funcs) {
		// Open files to which generated approximations/code is to be written.
		FILE *fp_omega = fopen((CODEGEN_FOLDER + coef.name + ".hpp").c_str(), "w");
		fprintf(fp_omega, preamble, coef.name, coef.name);

		// Write documentation and template function signatures.
		fprintf(fp_omega, proto, coef.name);

		// Write the actual generated approximations.
		SymEngine::Expression u("u"), v("v");
		for (int N = 0; N <= N_MAX; ++N)
			for (int M = 0; M <= M_MAX; ++M) {
				SymEngine::Expression om = SymEngine::expand(double_power_series(N, M, u, v, coef.n_min, coef.k_min, syms));  // Plain symbolic double sum.
				// Lambda for substituting symbolic coefficients in expression.
				auto subs = [N, M, &coef, &syms](const SymEngine::Expression& u_coef) -> SymEngine::Expression {
					std::map<SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCPBasicKeyLess> sub_map;
					for (int n = coef.n_min; n <= N; ++n)
						for (int m = coef.k_min; m <= M; ++m) {
							// Construct a SymEngine object directly from the mpreal coefficient.
							SymEngine::mpfr_class mc(coef.coef(n,m).mpfr_ptr());
							SymEngine::RCP<const SymEngine::Basic> mp_basic = SymEngine::real_mpfr(mc);
							sub_map[syms[n][m].get_basic()] = mp_basic;
						}
					return u_coef.subs(sub_map);
				};
				// Replace coefficient symbols with values and generate literal expression. This includes the expressions for the exponents.
				std::string expr = gen_2var_poly(om, u, v, subs);
				ostrstream om_func = expr2cfunc(N, M, expr, signature, coef.name);  // Generate a complete "C-function".
				fprintf(fp_omega, "%s", om_func.str().c_str());
			}

		// Clean up.
		fprintf(fp_omega, postamble, coef.name);
		fclose(fp_omega);
	}
}