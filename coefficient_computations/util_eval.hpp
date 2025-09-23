// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the BSD 2-Clause License.

#ifndef ECEF2GEO_POLYNOMIAL_UTIL_EVAL_HPP
#define ECEF2GEO_POLYNOMIAL_UTIL_EVAL_HPP

/*
 * Utility functions for polynomial evaluation code generation.
 */

#include "settings.hpp"

#include <string>  // std::string, std::to_string.
#include <vector>  // std::vector.
#include <map>  // std::map.
#include <iterator>  // std::distance.
#include <algorithm>  // std::min_element.

// Just to avoid long names.
using vec_map = std::vector<std::map<int, int>>;

/** Determine all solutions to the make-change problem.
 *
 * https://en.wikipedia.org/wiki/Change-making_problem
 *
 * @param total Value to which the "coins" should add upp.
 * @param coins "Coins" to use to add up to the total value.
 * @return A vector with all solutions, map of <coin value, number of coins>, to
 *         count to total.
 */
vec_map
make_change_all(int total, const std::vector<int> &coins) {
	// [[Counter()], [], [], [], ...]. "Counter()" is the way to count to 0.
	// ways[i] is all the ways to count to "i".
	std::vector<vec_map> ways(total + 1, vec_map());
	// An empty map is the way to count to zero.
	ways[0].emplace_back(std::map<int, int>());
//    ways = [[Counter()]] + [[] for _ in range(total)];
	for (int coin: coins) {
		// Start from "coin" since this is the lowest value one can count to
		// with coin.
		for (int i = coin; i <= total; i++) {
			// "way" is a Counter() of ways to count to "i-coin".
			for (std::map<int, int> way: ways[i - coin]) {
				// We can count to "i" by adding one coin to it.
				way[coin] += 1;
				// Add (list concatenation) this way to count to "i".
				ways[i].emplace_back(way);
			}
		}
	}
	// Ways contains all the ways to count to all values from "0" to "total".
	// Specifically pick out the ways for "total".
	return ways[total];
}

/** Solve the minimum make-change problem for the provided coins.
 *
 * https://en.wikipedia.org/wiki/Change-making_problem
 *
 * @param total Value to which the "coins" should add upp.
 * @param coins "Coins" to use to add up to the total value.
 * @return A map of <coin value, number of coins> which is the minimum
 *         make-change solution.
 */
std::map<int, int> make_change(int total, const std::vector<int> &coins) {
	vec_map ways = make_change_all(total, coins);
	// Count the number of coins for each possible way.
	std::vector<int> nr_coins(ways.size());
	for (int i = 0; i < ways.size(); ++i)
		for (const auto &_coins: ways[i])
			nr_coins[i] += _coins.second;
	// Select the way with the least number of coins.
	long minWay = std::distance(nr_coins.begin(),
								std::min_element(nr_coins.begin(),
												 nr_coins.end()));
	return ways[minWay];
}

/** Print C-code for evaluation of the polynomial.
 *
 * Recursively splits the polynomial in a lower part of degree n - f(n) and
 * higher part of degree f(n), i.e.
 * p(x) = b(x) + a(x) * x^f(n)
 *
 * @param c Polynomial coefficients
 * @param n Polynomial degree.
 * @param N Total polynomial degree.
 * @param var Polynomial variable.
 * @param f Polynomial evaluation schema (splitting). The first argument is the
 *          order of the (sub-) polynomial to split and the second argument is
 *          total order of the whole polynomial.
 * @param xs Output of the exponents contained in the polynomial, i.e. [1, 3]
 *           for x(1+x^3). This is such that code can be generated for computing
 *           them.
 * @return A C-code string representation for computing the polynomial.
 */
static std::string
print_poly(std::vector<std::string>::iterator c, int n, int N,
		   const std::string &var, int (*f)(int, int), std::vector<int> &xs) {
	if (n == 1)
		return *c;
	int k = f(n - 1, N);
	if (k == 0)
		return print_poly(c, n, N, var, evaluation_schema::direct_schema, xs);
	std::string a = print_poly(c + k, n - k, N, var, f, xs);
	std::string b = print_poly(c, k, N, var, f, xs);
	// FIXME(JO) Make this number comparison more robust.
	if (a == "0" && b == "0")
		return "0";
	if (a == "0")
		return b;
	if (k > 1)
		xs.emplace_back(k);
	std::string xk = k == 1 ? var : var + std::to_string(k);
	if (b == "0") {
		if (a.find(" + ") < a.length())
			return "(" + a + ")" + " * " + xk;
		else {
			if (a == "1")  // To not print "1 * x"
				return xk;
			return a + " * " + xk;
		}
	}
	// If there is a "+" in the expression, add parenthesis.
	if (a.find(" + ") < a.length())
		return b + " + " + "(" + a + ")" + " * " + xk;
	else {
		if (a == "1")  // To not print "1 * x"
			return b + " + " + xk;
		return b + " + " + a + " * " + xk;
	}
	// Gave significantly slower code.
//    if (a == "1")  // To not print "1 * x"
//        return b + " + " + xk;
//    return "std::fma(" + a + " ," + xk + " ," + b + ")";
}

/** Print C-code for the polynomial according to evaluation schema in settings.
 *
 * @param c polynomial coefficients.
 * @param var polynomial variable.
 * @param xs The exponents contained in the polynomial, i.e. [0, 1, 3] for
 *           x(1+x^3). This is such that code can be generated for computing
 *           them.
 * @return A C-code string representation for computing the polynomial.
 */
static std::string
print_poly(std::vector<std::string> c, const std::string &var,
		   std::vector<int> &xs) {
	return print_poly(c.begin(), c.size(), c.size(), var, POLY_EVAL_SCHEMA, xs);
}

/** Generate code for computing the exponents of xs.
 *
 * From lower to higher x^y, starts with x^1 and solves the change-making
 * problem of constructing the x^y from previously constructed x^y.
 *
 * @param xs The exponents to compute
 * @param var The variable name.
 * @return C-code for computing the exponents.
 */
static std::string print_exp(std::vector<int> &xs, const std::string &var) {

	// Internal struct just to improve naming and to avoid har to read
	// std::pair<const int, int>.
	struct base_and_exp {
		base_and_exp(std::pair<const int, int> p) : x(p.first), p(p.second) {}
		int x;  // base.
		int p;  // exponent.
	};

	if (xs.empty())
		return "";
	// Sort and remove duplicates.
	std::sort(xs.begin(), xs.end());
	xs.erase(unique(xs.begin(), xs.end()), xs.end());
	std::vector<int> comp;  // Computed exponents.
	std::string code;
	// Solve the make-change of computing x^y from previously computed x^y.
	std::vector<int> available(1, 1);  // "1" if for starting with "x^1".
	for (int y: xs) {
		std::map<int, int> way = make_change(y, available);
		code += " double " + var + std::to_string(y) + " = ";
		bool flag = false;
		for (const base_and_exp xp: way) {
			if (flag)
				code += " * ";
			flag = true;
			code += var + (xp.x != 1 ? std::to_string(xp.x) : "");
			for (int i = 1; i < xp.p; i++)
				code += " * " + var + (xp.x != 1 ? std::to_string(xp.x) : "");
		}
		code += ";";
		available.emplace_back(y);
	}
	return code;
}

#endif // ECEF2GEO_POLYNOMIAL_UTIL_EVAL_HPP
