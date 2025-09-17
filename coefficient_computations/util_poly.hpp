// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the BSD 2-Clause License.

#ifndef ECEF2GEO_POLYNOMIAL_UTIL_POLY_HPP
#define ECEF2GEO_POLYNOMIAL_UTIL_POLY_HPP

/*
 * Utility functions for manipulating monomial polynomials represented by mpfr
 * vectors.
 */

#include <mpfr.h>
#include <vector>  // std::vector

static void printVec(std::vector<mpfr_ptr> vec) {
	printf("[");
	for (mpfr_ptr e: vec)
		mpfr_printf("%.50Re, ", e);
	printf("]\n");
}

/** Add two polynomial represented by their coefficients.
 *
 * The function take ownership of arguments.
 * Precision is taken from the first argument and the second argument for those
 * elements for which there is none in the first argument.
 *
 * @param x First vector of polynomial coefficients.
 * @param y Second vector of polynomial coefficients.
 * @return
 */
static std::vector<mpfr_ptr>
poly_add(const std::vector<mpfr_ptr> &x, const std::vector<mpfr_ptr> &y) {
	size_t n = x.size() - 1;
	size_t m = y.size() - 1;
	std::vector<mpfr_ptr> res((n > m ? n : m) + 1);
	for (size_t i = 0, size = res.size(); i < size; ++i) {
		if (i > n)
			res[i] = y[i];
		else if (i > m)
			res[i] = x[i];
		else {
			res[i] = x[i];  // Reuse mpfr_t.
			mpfr_add(res[i], x[i], y[i], MPFR_RNDN);
			mpfr_clear(y[i]);
		}
	}
	return res;
}

/** Multiply two polynomial represented by their coefficients.
 *
 * The function take ownership of arguments.
 * Precision is taken from the first argument component of each element of the
 * result.
 *
 * @param x First vector of polynomial coefficients.
 * @param y Second vector of polynomial coefficients.
 * @param take_x Whether to take ownership of x.
 * @param take_y Whether to take ownership of y.
 * @return
 */
static std::vector<mpfr_ptr>
poly_mul(const std::vector<mpfr_ptr> &x, const std::vector<mpfr_ptr> &y,
		 bool take_x = true, bool take_y = true) {
	size_t n = x.size() - 1;
	size_t m = y.size() - 1;
	std::vector<mpfr_ptr> res(n + m + 1);
	for (size_t i = 0, size = n + m + 1; i < size; ++i) {
		size_t jfrom = i > m ? i - m : 0;
		size_t jto = std::min(n, i);
		mpfr_ptr zero;
		zero = new mpfr_t;
		mpfr_init2(zero, mpfr_get_prec(x[jfrom]));
		mpfr_set_d(zero, 0.0, MPFR_RNDN);
		res[i] = zero;
		for (size_t j = jfrom; j <= jto; ++j) {
			mpfr_t prod;
			mpfr_init2(prod, mpfr_get_prec(x[j]));
			mpfr_mul(prod, x[j], y[i - j], MPFR_RNDN);
			mpfr_add(res[i], res[i], prod, MPFR_RNDN);
			mpfr_clear(prod);
		}
	}
	// Release argument mpfr_t for which we have taken ownership.
	if (take_x)
		for (int i = 0; i <= n; ++i) {
			mpfr_clear(x[i]);
			delete[] x[i];
		}
	if (take_y)
		for (int i = 0; i <= m; ++i) {
			mpfr_clear(y[i]);
			delete[] y[i];
		}
	return res;
}

#endif // ECEF2GEO_POLYNOMIAL_UTIL_POLY_HPP
