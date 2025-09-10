// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the GPLv3 License.

/*
 * Utility functions for symbolic calculations.
 */

#ifndef ECEF2GEO_POLYNOMIAL_UTIL_SYM_HPP
#define ECEF2GEO_POLYNOMIAL_UTIL_SYM_HPP

#include "util_eval.hpp"

#include <ginac/ginac.h>  // For final polynomial manipulation.

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

/** Construct GiNaC symbolic expression for omega.
 *
 * @param N Summation index limit N.
 * @param M Summation index limit M.
 * @param u Polynomial variable u.
 * @param v Polynomial variable v.
 * @param syms Symbols to use for the coefficients arranged in a nested array.
 * @return A GiNaC expression for omega given the input parameters.
 */
GiNaC::ex phi_ex(int N, int M, GiNaC::symbol &u, GiNaC::symbol &v, int n_min, int k_min,
				 GiNaC::symbol syms[N_MAX + 1][M_MAX + 1]) {
	GiNaC::ex omega_NM;
	for (int n = n_min; n <= N; ++n) {
		GiNaC::ex sum;
		for (int k = k_min; k <= M; ++k)
			sum += syms[n][k] * GiNaC::pow(u, k);
		omega_NM += sum * GiNaC::pow(v, n);
	}
	return omega_NM;
}

/** Generate a C-code evaluation schema from 2-var polynomial.
 *
 * @param expr The polynomial expression.
 * @param u Inner polynomial variable of expression.
 * @param v Outer polynomial variable of expression.
 * @return A C-code string representation for evaluating expr.
 */
std::string gen_2var_poly(const GiNaC::ex &expr, const GiNaC::symbol &u,
						  const GiNaC::symbol &v,
						  const std::function<GiNaC::ex(GiNaC::ex)>& subs
) {
	// Retrieve "v-coefficients".
	GiNaC::lst v_coefs;
	for (int i = 0; i <= expr.degree(v); ++i)
		v_coefs.append(expr.coeff(v, i));
	// Loop over the "v-coefficients".
	std::vector<std::string> v_coefs_str;
	std::vector<int> us;
	for (int i = 0; i < v_coefs.nops(); ++i) {
		// Retrieve "u-coefficients".
		std::vector<std::string> u_coefs_str;
		for (int j = 0; j <= v_coefs[i].degree(u); ++j) {
			// Retrieve coefficient and substituting all symbols with values.
			GiNaC::ex u_coef = subs(v_coefs[i].coeff(u, j));
			// Convert u coefficient to a literal string like 0.345...
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
static std::ostringstream
expr2cfunc(int N, int M, const std::string &expr_str, const char *signature, const char *name) {
	std::ostringstream func;
	// TODO(JO) Using std::format would be preferable but support is so far
	//  limited.
	//  https://stackoverflow.com/questions/554063/how-do-i-print-a-double-value-with-full-precision-using-cout
	char signature_str[512];  // 512 characters should always be enough.
	std::sprintf(signature_str, signature, name, N, M);
	func << std::string(signature_str);
	func << expr_str << "; }" << std::endl;
	return func;
}

#endif //ECEF2GEO_POLYNOMIAL_UTIL_SYM_HPP
