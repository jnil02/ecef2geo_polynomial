// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the GPLv3 License.

#ifndef ECEF2GEO_POLYNOMIAL_UTIL_ROOT_HPP
#define ECEF2GEO_POLYNOMIAL_UTIL_ROOT_HPP

/*
 * Utility functions for root-finding.
 */

#include <functional>  // std::function.
#include <cmath>  // std::sqrt, std::pow.

#include "mplapack/mpreal.h"

using namespace mpfr;

/** ITP root search of f on the interval [a,b].
 *
 * If f(a) and f(b) have equal sign, the extremum closest to zero is returned.
 *
 * @param f Function to find the root of.
 * @param a Lower root interval limit.
 * @param b Upper root interval limit.
 * @param ya Function value at a.
 * @param yb Function value at b.
 * @param prec Bits of precision of arithmetics.
 * @param epsilon Precision below which the root should be found. If <=0.0, the
 *                maximum ulp of the interval [a,b] is used.
 * @param k1 R+
 * @param k2 in [1,1+phi] where phi = 1/2(1+sqrt(5)) is the golden ratio. k2 is
 *           the asymptotic order of convergence.
 * @param n0
 * @return The value x in [a,b] for which f(x) == 0 or NAN if no such value
 *         could be found.
 */
inline mpreal
itp_mpreal(std::function<mpreal(const mpreal &)> &f, mpreal a, mpreal b,
		   mpreal ya, mpreal yb, mp_prec_t prec,
		   const mpreal &epsilon = 1e-10,
		   const mpreal &k1 = 0.0,
		   const mpreal &k2 = 2.0, const mpreal &n0 = 1.0) {
	// Ensure prec.
	a.set_prec(prec);
	b.set_prec(prec);
	// Check that a < b
	if (a >= b) {
//        fprintf(stderr, "Warning. ITP root interval limits are swapped.");
		mpfr::swap(a, b);
		mpfr::swap(ya, yb);
	}
	// If k1 <=0 then the default value of k1, based on a and b.
	mpreal k1set;
	k1set.set_prec(prec);
	if (k1 <= 0.0) {
		k1set = 0.2 / (b - a);
	} else {
		k1set = k1;
	}
	if (k2 < 1.0 || k2 >= 1.0 + (1.0 + std::sqrt(5.0)) / 2.0) {
		fprintf(stderr, "k2 must be in [ 1, 1 + (1 + sqrt(5)) / 2 )");
		return NAN;
	}
	if (n0 < 0.0) {
		fprintf(stderr, "n0 must be non-negative");
		return NAN;
	}
	// If epsilon <= 0.0, set to largest ulp in interval.
	mpreal _epsilon;
	_epsilon.set_prec(prec);
	if (epsilon <= 0.0) {
		mpreal next1;
		mpreal next2;
		next1.set_prec(prec);
		next1 = mpfr::nextabove(a);
		next2.set_prec(prec);
		next2 = mpfr::nextbelow(b);
		// The largest ulp is in either limit.
		// IEEE-754 standard guarantee that two non-identical numbers does not
		// produce zero. I hope this is the case with mpfr?
		_epsilon = mpfr::max(next1 - a, b - next2);
	} else {
		_epsilon = epsilon;
	}

	// Check that f(a) and f(b) are finite.
	if (_isinf(ya) || _isnan(ya) || _isinf(yb) || _isnan(yb)) {
		fprintf(stderr, "f(a) and f(b) must be finite.");
		return NAN;
	}
	mpreal root, froot, estimprec;
	root.set_prec(prec);
	froot.set_prec(prec);
	estimprec.set_prec(prec);
	// Check whether (a, b) already satisfies the convergence criterion.
	if (b - a <= 2 * _epsilon) {
		root = (a + b) * 0.5;
		froot = f(root);//, pars) ;
		estimprec = (b - a) * 0.5;
		return root;
	}
	// Check whether the root lies on a limit of the input interval.
	if (ya == 0.0)
		return a;
	if (yb == 0.0)
		return b;
	// Check that f(a) and f(b) have opposite signs. Otherwise, return the
	// extremum closest to zero.
	double signya = (ya > 0.0) - (ya < 0.0);
	double signyb = (yb > 0.0) - (yb < 0.0);
	if (signya * signyb > 0.0)
		return mpfr::fabs(ya) <= mpfr::fabs(yb) ? a : b;
	// Set for_rk
	// M_LOG2E is 1/log_e(2). log_2(x) = log_e(x)/log_e(2) = log_e(x) * M_LOG2E
	mpreal log2e;
	log2e.set_prec(prec);
	log2e = log(_epsilon) * M_LOG2E;
	mpreal log2bma;
	log2bma.set_prec(prec);
	log2bma = log(b - a) * M_LOG2E;
	mpreal for_rk;
	for_rk.set_prec(prec);
	for_rk = std::pow(2.0, n0 - 1 + log2e + ceil(log2bma - log2e));
	// Initialise the counter k and implement the loop.
	int k = 0;
	mpreal xf, x12, delta, rk, xt, xITP, yITP, interval, tmp;
	xf.set_prec(prec);
	x12.set_prec(prec);
	delta.set_prec(prec);
	rk.set_prec(prec);
	xt.set_prec(prec);
	xITP.set_prec(prec);
	yITP.set_prec(prec);
	interval.set_prec(prec);
	tmp.set_prec(prec);
	interval = b - a;
	while (interval >= 2.0 * _epsilon) {
		// Interpolation. Regular falsi, equation (5)
		xf = (yb * a - ya * b) / (yb - ya);
		// Truncation
		x12 = (a + b) * 0.5;
		// Equation (13)
		int sigma = (x12 > xf) - (x12 < xf);
		delta = k1set * pow(b - a, k2);
		// Equation (14)
		if (delta <= fabs(x12 - xf)) {
			xt = xf + sigma * delta;
		} else {
			xt = x12;
		}
		// Projection, equation (15)
		rk = for_rk - (b - a) * 0.5;
		if (fabs(xt - x12) <= rk) {
			xITP = xt;
		} else {
			xITP = x12 - sigma * rk;
		}
		// Update (a, b)
		yITP = f(xITP);//, pars) ;
		if (yITP * signyb > 0.0) {
			b = xITP;
			yb = yITP;
		} else if (yITP * signyb < 0.0) {
			a = xITP;
			ya = yITP;
		} else {
			a = xITP;
			ya = yITP;
			b = xITP;
			yb = yITP;
		}
		root = (a + b) * 0.5;
		// Ensure that if the interval cannot become smaller, break.
		tmp = b - a;
		if (interval == tmp)
			break;
		interval = tmp;
		// Update the first term of rk
		for_rk *= 0.5;
		k += 1;
	}
	return root;
}

/** ITP root search of f on the interval [a,b].
 *
 * If f(a) and f(b) have equal sign, the extremum closest to zero is returned.
 *
 * @param f Function to find the root of.
 * @param a Lower root interval limit.
 * @param b Upper root interval limit.
 * @param prec Bits of precision of arithmetics.
 * @param epsilon Precision below which the root should be found. If <=0.0, the
 *                maximum ulp of the interval [a,b] is used.
 * @param k1 R+
 * @param k2 in [1,1+phi] where phi = 1/2(1+sqrt(5)) is the golden ratio. k2 is
 *           the asymptotic order of convergence.
 * @param n0
 * @return The value x in [a,b] for which f(x) == 0 or NAN if no such value
 *         could be found.
 */
inline mpreal
itp_mpreal(std::function<mpreal(const mpreal &)> &f, const mpreal &a,
		   const mpreal &b, mp_prec_t prec,
		   const mpreal &epsilon = 1e-10,
		   const mpreal &k1 = 0.0,
		   const mpreal &k2 = 2.0, const mpreal &n0 = 1.0) {
	mpreal ya, yb;
	ya.set_prec(prec);
	yb.set_prec(prec);
	ya = f(a);
	yb = f(b);
	return itp_mpreal(f, a, b, ya, yb, prec, epsilon, k1, k2, n0);
}



/** Bisection root search of f on the open interval (lo,hi).
 *
 * Values given to the function is of specified precision.
 * Computes the root to given precision.
 * Only evaluates f in [lo,hi].
 *
 * @param f Function to compute the root of. Assumed monotone over the interval
 *          (lo,hi).
 * @param lo Lower interval limit.
 * @param hi Upper interval limit.
 * @param f_lo Function value at lower limit.
 * @param f_hi Function value at upper limit.
 * @param prec Precision of the arithmetics.
 * @return The root of f on (lo,hi) or the f(lo) or f(hi) closest to 0.
 */
mpreal
bisection_search(std::function<mpreal(const mpreal &)> &f, mpreal lo, mpreal hi,
				 const mpreal &f_lo, const mpreal &f_hi, mp_prec_t prec) {
	if (f_hi < f_lo)
		mpfr::swap(lo, hi);
	// Initialize mid point value.
	mpreal mid;
	mid.set_prec(prec);
	mid = (lo + hi) / 2;
	// Each loop adds one bit of precision.
	for (int i = 0; i < prec + 1; ++i) {
		mpreal val(f(mid));
		// Bracket up.
		if (val < 0)
			lo = mid;
			// Bracket down.
		else
			hi = mid;
		// Update mid point.
		mid = (lo + hi) / 2;
	}
	return mid;
}

/** Bisection root search of f on the open interval (lo,hi).
 *
 * Values given to the function is of specified precision.
 * Computes the root to given precision.
 * Only evaluates f in [lo,hi].
 *
 * @param f Function to compute the root of. Assumed monotone over the interval
 *          (lo,hi).
 * @param lo Lower interval limit.
 * @param hi Upper interval limit.
 * @param prec Precision of the arithmetics.
 * @return The root of f on (lo,hi) or the f(lo) or f(hi) closest to 0.
 */
mpreal
bisection_search(std::function<mpreal(const mpreal &)> &f, const mpreal &lo,
				 const mpreal &hi, mp_prec_t prec) {
	mpreal f_lo, f_hi;
	f_lo.set_prec(prec);
	f_hi.set_prec(prec);
	f_lo = f(lo);
	f_hi = f(hi);
	return bisection_search(f, lo, hi, f_lo, f_hi, prec);
}

/** Find extremum value for a concave function over the interval (lo,hi).
 *
 * @param f The function for which the extremum is to be found.
 * @param lo Lower range limit.
 * @param hi Upper range limit.
 * @param h Distance over which the function has a smooth derivative.
 * @param prec Precision of the arithmetics and the values given to f.
 * @return The extremum value of f.
 */
mpreal find_max(std::function<mpreal(const mpreal &)> &f, const mpreal &h,
				mpreal lo, mpreal hi, mp_prec_t prec) {
	// Finite differences.
	std::function<mpreal(const mpreal &)> forward_difference = [&f, &h](
			const mpreal &x) {
		return ((f(x + h) - f(x)) / h);
	};
	std::function<mpreal(const mpreal &)> backward_difference = [&f, &h](
			const mpreal &x) {
		return ((f(x) - f(x - h)) / h);
	};
	std::function<mpreal(const mpreal &)> midpoint_difference = [&f, &h](
			const mpreal &x) {
		return ((f(x + h) - f(x - h)) / (h * 2));
	};

	lo.set_prec(prec);
	hi.set_prec(prec);
	mpreal res, f_lo, f_hi;
	res.set_prec(prec);
	f_lo.set_prec(prec);
	f_hi.set_prec(prec);
	f_lo = forward_difference(lo);
	f_hi = backward_difference(hi);
//    res = bisection_search(midpoint_difference, lo, hi, f_lo, f_hi, prec);
	res = itp_mpreal(midpoint_difference, lo, hi, f_lo, f_hi, prec);
	return res;
}

#endif // ECEF2GEO_POLYNOMIAL_UTIL_ROOT_HPP
