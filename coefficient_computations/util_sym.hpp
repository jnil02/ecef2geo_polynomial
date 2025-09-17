// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the BSD 2-Clause License.

/*
 * Utility functions for symbolic calculations.
 */

#ifndef ECEF2GEO_POLYNOMIAL_UTIL_SYM_HPP
#define ECEF2GEO_POLYNOMIAL_UTIL_SYM_HPP

#include "util_eval.hpp"

#include <iomanip>  // std::setprecision.
#include <stdexcept>  // std::invalid_argument.
#include <utility>
#include <symengine/expression.h>
#include <symengine/functions.h>

using SymEngine::Expression;
using SymEngine::integer;
using SymEngine::rational;
using SymEngine::pow;
using SymEngine::binomial;
using SymEngine::Mul;
using SymEngine::Add;
using SymEngine::Pow;
using SymEngine::Symbol;
using SymEngine::is_a_Number;

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

/** Construct SymEngine symbolic expression for mu.
 *
 * @param N Summation index limit N.
 * @param M Summation index limit M.
 * @param u Polynomial variable u.
 * @param v Polynomial variable v.
 * @param syms Symbols to use for the coefficients arranged in a nested array.
 * @return A GiNaC expression for mu given the input parameters.
 */
Expression mu_ex(int N, int M, Expression &u, Expression &v,
				 Expression syms[N_MAX + 1][M_MAX + 1]) {
	Expression mu_NM = u;
	for (int m = 0; m <= M; ++m)
		mu_NM += rational(1, 2) * syms[0][m] * pow(u, m);
	for (int n = 1; n <= N; ++n) {
		Expression sum(0);
		for (int m = 0; m <= M; ++m)
			sum += syms[n][m] * pow(u, m);
		for (int k = 0; k <= n; ++k)
			mu_NM += Expression(powm1(k))
					 * binomial(SymEngine::Integer(2 * n), 2 * k)
					 * pow(v, k)
					 * pow(integer(1) - v, n - k)
					 * sum;
	}
	return mu_NM;
}

/** Construct SymEngine symbolic expression for omega.
 *
 * @param N Summation index limit N.
 * @param M Summation index limit M.
 * @param u Polynomial variable u.
 * @param v Polynomial variable v.
 * @param syms Symbols to use for the coefficients arranged in a nested array.
 * @return A GiNaC expression for omega given the input parameters.
 */
Expression omega_ex(int N, int M, Expression &u, Expression &v,
					Expression syms[N_MAX + 1][M_MAX + 1]) {
	Expression omega_NM(0);
	for (int n = 1; n <= N; ++n) {
		Expression sum(0);
		for (int m = 0; m <= M; ++m)
			sum += syms[n][m] * pow(u, m);
		for (int k = 0; k < n; ++k)
			omega_NM += Expression(powm1(k))
						* binomial(SymEngine::Integer(2 * n), 2 * k + 1)
						* pow(v, k)
						* pow(integer(1) - v, n - k - 1)
						* sum;
	}
	return omega_NM;
}

/** Construct SymEngine double power series.
 *
 * @param N Summation index limit N.
 * @param M Summation index limit M.
 * @param u Polynomial variable u.
 * @param v Polynomial variable v.
 * @param syms Symbols to use for the coefficients arranged in a nested array.
 * @return A GiNaC expression for omega given the input parameters.
 */
Expression double_power_series(int N, int M, Expression &u, Expression &v,
							   int n_min, int k_min,
							   Expression syms[N_MAX + 1][M_MAX + 1]) {
	Expression omega_NM(0);
	for (int n = n_min; n <= N; ++n) {
		Expression sum(0);
		for (int k = k_min; k <= M; ++k)
			sum += syms[n][k] * pow(u, k);
		omega_NM += sum * pow(v, n);
	}
	return omega_NM;
}

/**
 * Class for iterating over the coefficients of a polynomial.
 */
class CoefRange {
public:
	/** Construct the coefficient range.
	 *
	 * @param expr The expression with coefficients to iterate over.
	 * @param x The variable with respect to which the coefficients are defined.
	 */
	CoefRange(const Expression &expr, Expression x) : x(std::move(x)) {

		// Extract all terms.
		Expression expanded = SymEngine::expand(expr);
		if (is_a<Add>(*expanded.get_basic())) {
			for (const auto &term : expanded.get_basic()->get_args())
				process_term(Expression(term));
		} else
			process_term(expanded);

		// Combine all terms into a complete list of coefficients.
		coeffs.resize(max_exp + 1, Expression(0));
		for (auto &kv : terms)
			coeffs[kv.first] += kv.second;
	}

	// Enable range-based for loops.
	[[nodiscard]] auto begin() const { return coeffs.begin(); }
	[[nodiscard]] auto end()   const { return coeffs.end(); }

private:
	Expression x; // The variable.
	std::vector<Expression> coeffs; // The coefficients.
	std::map<int, Expression> terms; // All terms powers and coefficients.
	int max_exp = 0; // The highest power.

	// Extract the coefficient and the powers for this term.
	void process_term(const Expression &term) {
		int power = 0;
		Expression rest(1);
		const auto &b = term.get_basic();

		if (is_a<Mul>(*b)) {
			for (const auto &f : b->get_args()) {
				if (is_a<Pow>(*f)) {
					auto pw = rcp_static_cast<const Pow>(f);
					if (Expression(pw->get_base()) == x) {
						auto &iexp = SymEngine::down_cast<const SymEngine::Integer &>(*pw->get_exp());
						power += (int)iexp.as_int();
					} else {
						rest *= Expression(f);
					}
				} else if (is_a<Symbol>(*f) && Expression(f) == x) {
					power += 1;
				} else {
					rest *= Expression(f);
				}
			}
		} else if (is_a<Pow>(*b)) {
			auto pw = rcp_static_cast<const Pow>(b);
			if (Expression(pw->get_base()) == x) {
				auto &iexp = SymEngine::down_cast<const SymEngine::Integer &>(*pw->get_exp());
				power = (int)iexp.as_int();
			} else {
				rest = term;
			}
		} else if (is_a<Symbol>(*b) && Expression(b) == x) {
			power = 1;
		} else {
			rest = term;
		}

		terms[power] += rest;
		if (power > max_exp)
			max_exp = power;
	}
};

// Free function to build CoefRange.
inline CoefRange coeffs(const Expression &expr, const Expression &sym) {
	return {expr, sym};
}

/** Generate a C-code evaluation schema from 2-var polynomial.
 *
 * @param expr The polynomial expression.
 * @param u Inner polynomial variable of expression.
 * @param v Outer polynomial variable of expression.
 * @return A C-code string representation for evaluating expr.
 */
std::string gen_2var_poly(const Expression &expr, const Expression &u,
						  const Expression &v,
						  const std::function<Expression(Expression)> &subs) {
	std::vector<std::string> v_coefs_str;
	std::vector<int> us;

	// Loop over v-coefficients.
	for (const auto& v_coef : coeffs(expr, v)) {
		std::vector<std::string> u_coefs_str;

		// Loop over u-coefficients.
		for (const auto& uc : coeffs(v_coef, u)) {
			Expression u_coef = subs(uc);

			std::ostringstream s;
			s << std::setprecision(17);
			if (is_a_Number(*u_coef.get_basic())) {
				// TODO(JO) It would be nice to use std::format but support is
				//  so far limited.
				//  https://stackoverflow.com/questions/554063/how-do-i-print-a-double-value-with-full-precision-using-cout
				// s << std::format("{}", GiNaC::ex_to<GiNaC::numeric>(omega_u_coef).to_double());
				double val = eval_double(u_coef);
				// FIXME(JO) This does not ensure we get a decimal dot, i.e. a double constant.
				s << val;
			} else {
				// Uncomment if non-numeric coefficients should be allowed.
				// s << val;
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
