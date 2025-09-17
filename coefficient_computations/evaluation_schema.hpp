// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the BSD 2-Clause License.

#ifndef ECEF2GEO_POLYNOMIAL_EVALUATION_SCHEMA_HPP
#define ECEF2GEO_POLYNOMIAL_EVALUATION_SCHEMA_HPP

/*
 * Polynomial evaluation schemas such as Horner, Estrin etc. These are meant to
 * be provided as arguments to the print_poly function. However, in order to
 * enable static settings without including the functionality itself, these
 * reside in this separate file.
 *
 * The polynomials are split by calling the evaluation schema on the polynomial
 * and then recursively on the sub-polynomials. Denote the schema with f(k,N).
 * Then each call splits the (sub-) polynomial p(x) in a lower part b(x) of
 * degree k - f(k) and higher part a(x) of degree f(k,N), i.e.
 *
 * p(x) = b(x) + a(x) * x^f(k,N)
 *
 * The first argument k of the evaluation schema is the order of the (sub-)
 * polynomial to split and the second argument N is total order of the whole
 * polynomial. This means that the returned value must be smaller than k.
 */

#include <cmath> // for std::floor and std::log2.

namespace evaluation_schema {

/** Horner(-1) evaluation schema.
 *
 * @param k Order of the sub-polynomial to be split.
 * @param N Total order of the polynomial.
 * @return The split of the sub-polynomial.
 */
int horner_schema(int k, int N) {
	return 1;
}

/** Horner-2 evaluation schema with a 2nd degree innermost part.
 *
 * @param k Order of the sub-polynomial to be split.
 * @param N Total order of the polynomial.
 * @return The split of the sub-polynomial.
 */
int horner2_inner_schema(int k, int N) {
	// This will determine the "phase" of the split.
	// This gives a 2nd degree innermost part. Without it, this gives a 2nd
	// degree outermost part.
	if (!(N & 1) && k == N - 1)
		return 1;
	if (k == 1)
		return 1;
	return 2;
}

/** Horner-2 evaluation schema with a 2nd degree outermost part.
 *
 * @param k Order of the sub-polynomial to be split.
 * @param N Total order of the polynomial.
 * @return The split of the sub-polynomial.
 */
int horner2_outer_schema(int k, int N) {
	if (k == 1)
		return 1;
	return 2;
}

/** Direct (monomial) evaluation schema.
 *
 * @param k Order of the sub-polynomial to be split.
 * @param N Total order of the polynomial.
 * @return The split of the sub-polynomial.
 */
int direct_schema(int k, int N) {
	return k;
}

// Exponentiation by squaring.
template<class T>
inline T _ipow(T base, int exp) {
	T result = 1;
	for (;;) {
		if (exp & 1)
			result *= base;
		exp >>= 1;
		if (!exp)
			break;
		base *= base;
	}
	return result;
}

/** Estrin evaluation schema.
 *
 * @param k Order of the sub-polynomial to be split.
 * @param N Total order of the polynomial.
 * @return The split of the sub-polynomial.
 */
int estrin_schema(int k, int N) {
	return _ipow(2, (int) std::floor(std::log2(k)));
}

/** "Balanced" evaluation schema.
 *
 * https://arxiv.org/abs/1307.5655v2
 *
 * @param k Order of the sub-polynomial to be split.
 * @param N Total order of the polynomial.
 * @return The split of the sub-polynomial.
 */
int balanced_schema(int k, int N) {
	return k / 2;
}

}  // namespace evaluation_schema

#endif // ECEF2GEO_POLYNOMIAL_EVALUATION_SCHEMA_HPP