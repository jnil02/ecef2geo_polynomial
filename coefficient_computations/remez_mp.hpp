// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the GPLv3 License.

#ifndef ECEF2GEO_POLYNOMIAL_REMEZ_MP_HPP
#define ECEF2GEO_POLYNOMIAL_REMEZ_MP_HPP

/*
 * Basic multiprecision Remez exchange algorithm implementation and supporting
 * functions.
 */

#include <functional>  // std::function.

#include "util_root.hpp"
#include "util_cheb.hpp"

#include "mplapack/mpblas_mpfr.h"
#include "mplapack/mplapack_mpfr.h"

/** Naive multiprecision Remez exchange algorithm implementation.
 *
 * Only provides absolute error minimization.
 * Only provides minimization of a full polynomial up to a certain degree.
 *
 * @param c Output of the N polynomial coefficients of the approximation.
 * @param N Number of coefficients (degree of the polynomial minus 1).
 * @param f The function to approximate. The function should return a value of
 *          the same precision as the provided precision. May be evaluated over
 *          the range [lo,hi].
 * @param h Distance over which the derivative of the function is smooth. This
 *          is used to compute a numerical derivatives. Should have the
 *          provided precision. Selecting a proper h is a non-trivial problem,
 *          see https://en.wikipedia.org/wiki/Numerical_differentiation and may
 *          require some experimentation.
 * @param lo Lower bound of function interval over which the function should be
 *           approximated.
 * @param hi Upper bound of function interval over which the function should be
 *           approximated.
 * @param x_init Whether x contains initial points.
 * @param rel_err_limit Minimax relative error between residual error extrema
 *                      for which the algorithm terminates.
 * @param max_iter Maximum number of iteration until termination.
 * @param prec Precision of the arithmetics.
 */
void remez(mpreal *c, int N, std::function<mpreal(const mpreal &)> f,
		   const mpreal &h, const mpreal &lo, const mpreal &hi,
		   const mpreal &rel_err_limit, int max_iter, mp_prec_t prec) {

	// Set precision of output variables.
	for (int i = 0; i < N; ++i)
		c[i].set_prec(prec);

	mpreal *x = new mpreal[N + 1];
	for (int i = 0; i < N + 1; ++i)
		x[i].set_prec(prec);
	// Column-major (Fortran) matrix.
	mpreal *A = new mpreal[(N + 1) * (N + 1)];
	for (int i = 0; i < (N + 1) * (N + 1); ++i)
		A[i].set_prec(prec);
	mpreal *vander = new mpreal[(N + 1) * N];
	for (int i = 0; i < (N + 1) * N; ++i)
		vander[i].set_prec(prec);
	mpreal *b = new mpreal[N + 1];
	for (int i = 0; i < N + 1; ++i)
		b[i].set_prec(prec);
	mpreal *zeros = new mpreal[N + 2];
	for (int j = 0; j < N + 2; ++j)
		zeros[j].set_prec(prec);
	mpreal *extrema = new mpreal[N + 1];
	for (int j = 0; j < N + 1; ++j)
		extrema[j].set_prec(prec);
	mpreal *res = new mpreal[N + 1];
	for (int j = 0; j < N + 1; ++j)
		res[j].set_prec(prec);
	mplapackint *ipiv = new mplapackint[N + 1];

	// Initialize the node points with Chebyshev nodes.
	chebyshev_nodes(x, N + 1, lo, hi);

	zeros[0] = hi;
	zeros[N + 1] = lo;

	extrema[0] = hi;
	extrema[N] = lo;

	mplapackint _n = N + 1;
	mplapackint info;

	mpreal mae;
	mae.set_prec(prec);
	mpreal max_err;
	max_err.set_prec(prec);
	mpreal max_rel_err;
	max_rel_err.set_prec(prec);

	// Remez exchange algorithm iterations.
	for (int i = 0; i < max_iter; ++i) {

		// These values have to be written every time since A is used for
		// output LD-decomposition in dgesv.
		for (int j = 0; j < N + 1; ++j)
			A[j + N * (N + 1)] = j % 2 ? 1 : -1;

		// Build the rest of the A matrix.
		cheb_vander(vander, x, N - 1, N + 1);
		for (int j = 0; j < N + 1; ++j)  // Loop over rows.
			for (int k = 0; k < N; ++k)  // Loop over column 1 to
				A[j + (N + 1) * k] = vander[j + (N + 1) * k];

		// Compute the function values for the current points of interest.
		for (int j = 0; j < N + 1; ++j)
			b[j] = f(x[j]);

		// Solve the linear equation system.
		Rgesv(_n, (mplapackint) 1, A, _n, ipiv, b, _n, info);

		// Copy the new coefficients, ignoring the last "E".
		for (int j = 0; j < N; ++j)
			c[j] = b[j];

		// Lambda for the residual error.
		std::function<mpreal(const mpreal &)> re = [&f, &c, &N](
				const mpreal &x) {
			return f(x) - chebval(x, c, N);
		};

		// Find residual error roots of the updated coefficients.
		for (int j = 1; j < N + 1; ++j)
//            zeros[j] = bisection_search(r_i, x[j - 1], x[j], prec);
			zeros[j] = itp_mpreal(re, x[j - 1], x[j], prec, 0.0);

		// Find residual error extrema of the updated coefficients.
		for (int j = 1; j < N; ++j)
			extrema[j] = find_max(re, h, zeros[j], zeros[j + 1], prec);

		// Compute the absolute residual errors.
		for (int j = 0; j < N + 1; ++j)
			res[j] = abs(re(extrema[j]), MPFR_RNDN);

		mae = sum(res, N + 1, MPFR_RNDN) / (N + 1);
		max_err = 0.0;
		max_err.set_prec(prec);
		max_rel_err = 0.0;
		max_rel_err.set_prec(prec);
		for (int j = 0; j < N + 1; ++j) {
			max_err = max(max_err, res[j]);
			max_rel_err = max(max_rel_err,
							  abs((res[j] - mae) / mae, MPFR_RNDN));
		}
		if (max_rel_err < rel_err_limit)
			break;

		// Set the new points of interest.
		for (int j = 0; j < N + 1; ++j)
			x[j] = extrema[j];
	}

	// Convert back to monomial polynomial basis.
	cheb2mono(c, c, N, prec);

	delete[] x;
	delete[] A;
	delete[] vander;
	delete[] b;
	delete[] zeros;
	delete[] extrema;
	delete[] res;
	delete[] ipiv;
}

#endif // ECEF2GEO_POLYNOMIAL_REMEZ_MP_HPP
