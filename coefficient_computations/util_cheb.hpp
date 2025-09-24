// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the BSD 2-Clause License.

#ifndef ECEF2GEO_POLYNOMIAL_UTIL_CHEB_HPP
#define ECEF2GEO_POLYNOMIAL_UTIL_CHEB_HPP

/*
 * Utility functions for, primarily, Chebyshev polynomials manipulation. The
 * polynomials coefficients are represented by mpreal numbers.
 */

// Note this file assumes either of below has been included.
//#include <mplapack/mpreal.h>
//#include <mpreal.h>

using namespace mpfr;

void print_mpreal_array(mpreal *a, int size) {
	printf("[");
	for (int i = 0; i < size; ++i)
		printf("%s,", a[i].to_string().c_str());
	printf("]\n");
}

// Prints column-major matrix
void print_mpreal_matrix(mpreal *a, int rows, int cols) {
	printf("[\n");
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j)
			printf("%s,", a[i + j * rows].to_string().c_str());
		printf("],\n");
	}
	printf("]\n");
}

/** Polynomial Chebyshev to monomial basis.
 *
 * @param c_mono Output polynomial coefficient in monomial basis.
 * @param c_cheb Input polynomial coefficients in Chebyshev basis.
 * @param n Polynomial degree.
 * @param prec
 */
void cheb2mono(mpreal *c_mono, mpreal *c_cheb, int n, mp_prec_t prec) {
	if (n < 3) {
		for (int i = 0; i < n; ++i)
			c_mono[i] = c_cheb[i];
	} else {
		mpreal *c0 = new mpreal[n];
		mpreal *c1 = new mpreal[n];
		mpreal *c0_tmp = new mpreal[n];
		mpreal *c1_tmp = new mpreal[n];
		for (int i = 0; i < n; ++i) {
			c0[i].set_prec(prec);
			c1[i].set_prec(prec);
			c0_tmp[i].set_prec(prec);
			c1_tmp[i].set_prec(prec);
		}

		c0[0] = c_cheb[n - 2];
		c1[0] = c_cheb[n - 1];
		int size_c0 = 1;
		int size_c1 = 1;
		for (int i = n - 1; i > 1; --i) {
			// Store c0 and c1.
			for (int j = 0; j < size_c0; ++j)
				c0_tmp[j] = c0[j];
			for (int j = 0; j < size_c1; ++j)
				c1_tmp[j] = c1[j];
			// Polynomial subtraction c0 = c[i-2] - c1.
			c0[0] = c_cheb[i - 2] - c1_tmp[0];
			for (int j = 1; j < size_c1; ++j)
				c0[j] = -c1_tmp[j];
			// Polynomial addition c1 = tmp + x * c1 * 2.
			c1[0] = c0_tmp[0];
			for (int j = 1; j < size_c0; ++j)
				c1[j] = c0_tmp[j] + c1_tmp[j - 1] * 2;
			for (int j = size_c0; j < size_c1; ++j)
				c1[j] = c1_tmp[j - 1] * 2;
			c1[size_c1] = c1_tmp[size_c1 - 1] * 2;

			size_c0 = size_c1;
			size_c1++;
		}

		// c_poly = c0 + x * c1.
		c_mono[0] = c0[0];
		for (int j = 1; j < size_c0; ++j)
			c_mono[j] = c0[j] + c1[j - 1];
		for (int j = size_c0; j < size_c1; ++j)
			c_mono[j] = c1[j - 1];
		c_mono[size_c1] = c1[size_c1 - 1];
	}
}

/** Generate a set of Chebyshev nodes in the range [lower, upper].
 *
 * https://en.wikipedia.org/wiki/Chebyshev_nodes
 * The points are computed to the precision of the first variable in the output
 * array.
 *
 * @param x Output array of size N. Most be initialized.
 * @param n Number of points.
 * @param lo Lower range limit.
 * @param hi Upper range limit.
 * @param prec Precision to which the points are computed.
 */
void chebyshev_nodes(mpreal *x, int n, const mpreal &lo, const mpreal &hi) {
	mp_prec_t prec = x->get_prec();
	mpreal signed_range = hi - lo;
	mpreal range = fabs(signed_range, MPFR_RNDN);
	for (int i = 0; i < n; ++i) {
		mpreal _x;
		_x.set_prec(prec);
		_x = const_pi(prec, MPFR_RNDN) * (2 * i + 1) / (2 * n);
		x[i] = (cos(_x) + 1.0) * 0.5 * range + lo;
	}
}

/** Evaluate polynomial in Chebyshev basis.
 *
 * The precision of the arithmetic is the same as that of the provided point.
 *
 * @param x Point to evaluate polynomial in.
 * @param c Chebyshev basis coefficients.
 * @param size_c Number of coefficients.
 * @return polynomial value in provided point.
 */
mpreal chebval(const mpreal &x, mpreal *c, int size_c) {
	mp_prec_t prec = x.get_prec();
	mpreal c0;
	c0.set_prec(prec);
	mpreal c1;
	c1.set_prec(prec);
	if (size_c == 1) {
		c0 = c[0];
		c1 = 0.0;
	} else if (size_c == 2) {
		c0 = c[0];
		c1 = c[1];
	} else {
		mpreal x2 = mpreal(2 * x);
		c0 = c[size_c - 2];
		c1 = c[size_c - 1];
		for (int i = 3; i < size_c + 1; ++i) {
			mpreal tmp = mpreal(c0);
			c0 = c[size_c - i] - c1;
			c1 = tmp + c1 * x2;
		}
	}
	return c0 + c1 * x;
}

/** Pseudo-Vandermonde matrix of given degree.
 *
 * @param v Column major matrix of size (size) x (ideg+1)
 * @param x points of polynomial os size (size).
 * @param ideg degree of polynomial (cols-1)
 * @param size size of x (rows).
 */
void cheb_vander(mpreal *v, mpreal *x, int ideg, int size) {
	for (int i = 0; i < size; ++i)
		v[i] = mpreal(1.0);  // Fill in first column.
	if (ideg > 0) {
		mpreal *x2 = new mpreal[size];
		for (int i = 0; i < size; ++i) {
			v[i + size] = x[i];  // Fill in second column.
			x2[i] = 2. * x[i];
		}
		// Forward recursion for remaining columns.
		for (int j = 2; j < ideg + 1; ++j) {  // Loop over cols.
			for (int i = 0; i < size; ++i) { // Loop over rows.
				v[i + j * size] =
						v[i + (j - 1) * size] * x2[i] - v[i + (j - 2) * size];
			}
		}
		delete[] x2;
	}
}


/** Multiply a polynomial in Chebyshev basis with x.
 *
 * The input and output arrays can be the same.
 *
 * @param c_out n+1 long array of resulting Chebyshev coefficients. Should be
 *              initialized the desired arithmetic precision.
 * @param c_in n long array of polynomial Chebyshev coefficients.
 * @param n The number of coefficients in the original polynomial.
 */
int chebmulx(mpreal *c_out, mpreal *c_in, int n) {
	// The zero series need special treatment.
	if (n == 1 && c_in[0] == 0.0) {
		c_out[0] = c_in[0];
		return 1;
	}
	if (n > 1) {
		mpreal tmp[n - 1];
		for (int i = 0; i < n - 1; ++i)
			tmp[i] = c_in[i + 1] / 2;
		c_out[1] = c_in[0];
		c_out[0] = 0.0;
		for (int i = 0; i < n - 1; ++i) {
			c_out[i + 2] = tmp[i];
			c_out[i] += tmp[i];
		}
	} else {
		c_out[1] = c_in[0];
		c_out[0] = 0.0;
	}
	return n + 1;
	// c_out[0] =           c_in[1]/2
	// c_out[1] = c_in[0] + c_in[2]/2
	// c_out[2] = c_in[1]/2 + c_in[3]/2
	// ...
	// c_out[n-2] = c_in[n-3]/2 + c_in[n-1]/2
}

/** Add two Chebyshev polynomials, i.e. just adding the coefficients.
 *
 * @param c_sum Polynomial coefficients of polynomial sum.
 * @param c1 Coefficient of first polynomial to add.
 * @param c2 Coefficient of second polynomial to add.
 * @param n1 Degree of first polynomial.
 * @param n2 Degree of second polynomial.
 */
void chebadd(mpreal *c_sum, mpreal *c1, mpreal *c2, int n1, int n2) {
	int n_min = std::min(n1, n2);
	for (int i = 0; i < n_min; ++i)
		c_sum[i] = c1[i] + c2[i];
	for (int i = std::min(n_min, n1); i < n1; ++i)
		c_sum[i] = c1[i];
	for (int i = std::min(n_min, n2); i < n2; ++i)
		c_sum[i] = c2[i];
}

/** Polynomial monomial to Chebyshev basis.
 *
 * @param c_mono Input polynomial coefficient in monomial basis.
 * @param c_cheb Output polynomial coefficients in Chebyshev basis.
 * @param n Polynomial degree.
 * @param prec
 */
void mono2cheb(mpreal *c_cheb, mpreal *c_mono, int n, mp_prec_t prec) {
	int n_cheb = 1;
	c_cheb[0] = 0.0;
	for (int i = n - 1; i >= 0; --i) {
		n_cheb = chebmulx(c_cheb, c_cheb, n_cheb);
		chebadd(c_cheb, c_cheb, c_mono + i, n_cheb, 1);
	}
}

/** Evaluate a Chebyshev polynomial.
 *
 * @param x Point to evaluate the polynomial in.
 * @param n The degree of the polynomial.
 * @return polynomial value in provided point.
 */
mpreal chebval(const mpreal &x, int n) {
	mpreal Tnm2 = mpreal(1);
	if (n == 0)
		return Tnm2;
	mpreal Tnm1 = x;
	if (n == 1)
		return Tnm1;
	for (int i = 2; i <= n; ++i) {
		mpreal Tn = Tnm1 * x * 2 - Tnm2;
		Tnm2 = Tnm1;
		Tnm1 = Tn;
	}
	return Tnm1;
}

// Map x from range [min, max] to [min_to, max_to]
mpreal range_map(const mpreal &x, const mpreal &min_x, const mpreal &max_x,
				 const mpreal &min_to, const mpreal &max_to) {
	return (x - min_x) / (max_x - min_x) * (max_to - min_to) + min_to;
}

/** Computes the n first Chebyshev coefficients of a function.
 *
 * @param f The function.
 * @param n The number of coefficients to compute.
 * @param min The lower limits of the f input.
 * @param max The upper limit of the f input.
 * @param coef Allocated output array of coefficients of size at least n.
 */
void chebcoef(mpreal (*f)(mpreal), int n, const mpreal &min, const mpreal &max,
			  mpreal *coef) {
	for (int i = 0; i < n; ++i)
		coef[i] = mpreal(0);
	for (int i = 0; i < n; i++) {
		double val =
				f(range_map(cos(M_PI * (i + .5) / n), -1, 1, min, max)) * 2 / n;
		for (int j = 0; j < n; j++)
			coef[j] += val * cos(M_PI * j * (i + .5) / n);
	}
}

mpreal chebcoef_single(std::function<mpreal(const mpreal &)> f, int n,
					   const mpreal &min, const mpreal &max) {
	mpreal coef = mpreal(0);
	for (int i = 0; i < n; i++) {
		mpreal val =
				f(range_map(cos(M_PI * (i + .5) / n), -1, 1, min, max)) * 2 / n;
		coef += val * cos(M_PI * (n - 1) * (i + .5) / n);
	}
	return coef;
}

#endif // ECEF2GEO_POLYNOMIAL_UTIL_CHEB_HPP
