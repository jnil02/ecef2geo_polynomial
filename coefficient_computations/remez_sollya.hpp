// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the GPLv3 License.

#ifndef ECEF2GEO_POLYNOMIAL_REMEZ_SOLLYA_HPP
#define ECEF2GEO_POLYNOMIAL_REMEZ_SOLLYA_HPP

/*
 * Functions for interfacing Sollya's minimax implementation and related mpfr
 * and polynomial manipulation.
 */

#include "util_poly.hpp"

#include <vector>  // std::vector
#include <sollya.h>
#include <mpfr.h>
#include <stdexcept>  // for std::invalid_argument exception.

/** Retrieve mpfr_t constant from sollya_obj_t.
 *
 * sollya_lib_get_head_function(&type, f) should give
 * SOLLYA_BASE_FUNC_CONSTANT. The function will return an error if this is not
 * the case.
 *
 * @param c constant output variable.
 * @param f the Sollya object containing a constant.
 * @return Non-zero if successful.
 */
static int get_mpfr_constant(mpfr_t c, sollya_obj_t &f) {
	mp_prec_t prec;
	int res = sollya_lib_get_prec_of_constant(&prec, f);
	if (res) {
		mpfr_init2(c, prec);
		sollya_lib_get_constant(c, f); // Exact conversion.
	} else {
		mpfr_init2(c, 200); // Initialization at some default precision.
		res = sollya_lib_get_constant(c, f);
		if (!res)
			sollya_lib_printf("Error: %b is not a constant expression\n", f);
	}
	return res;
}

/** Recursive retrieval of polynomial coefficients from a Sollya arithmetic tree.
 *
 * The caller should call both mpfr_clear and delete[] on all elements of the
 * returned vector.
 *
 * @param f Sollya arithmetic tree consisting of +-*^.
 * @return vector of ordered polynomial coefficients including zeros for
 *         lacking terms.
 */
static std::vector<mpfr_ptr> get_poly_coefs(sollya_obj_t &f) {
	if (!f)
		return {};
	sollya_base_function_t type;
	sollya_lib_get_head_function(&type, f);
	switch (type) {
		case SOLLYA_BASE_FUNC_CONSTANT: {
			mpfr_ptr a;
			a = new mpfr_t;
			get_mpfr_constant(a, f);
			std::vector<mpfr_ptr> result(1);
			result[0] = a;
			return result;
		}
		case SOLLYA_BASE_FUNC_FREE_VARIABLE: {
			std::vector<mpfr_ptr> result(2);
			mpfr_ptr zero, one;
			zero = new mpfr_t;
			one = new mpfr_t;
			mpfr_init2(zero, mpfr_get_default_prec());
			mpfr_init2(one, mpfr_get_default_prec());
			mpfr_set_d(zero, 0.0, MPFR_RNDN);
			mpfr_set_d(one, 1.0, MPFR_RNDN);
			result[0] = zero;
			result[1] = one;
			return result;
		}
		case SOLLYA_BASE_FUNC_POW: {
			sollya_obj_t g_1 = NULL;
			sollya_obj_t g_2 = NULL;
			int n;
			sollya_lib_get_subfunctions(f, &n, &g_1, &g_2, NULL);
			int exp;
			sollya_lib_get_constant_as_int(&exp, g_2);
			std::vector<mpfr_ptr> base = get_poly_coefs(g_1);
			std::vector<mpfr_ptr> result(base.size());
			// Allocate and copy a new set base values.
			for (int i = 0; i < result.size(); ++i) {
				mpfr_ptr e;
				e = new mpfr_t;
				mpfr_init_set(e, base[i], MPFR_RNDN);
				result[i] = e;
			}
			for (int i = 1; i < exp; ++i)
				result = poly_mul(result, base, true, false);
			// Free base values which we have the ownership of.
			for (int i = 0; i < base.size(); ++i) {
				mpfr_clear(base[i]);
				delete[] base[i];
			}
			return result;
		}
		case SOLLYA_BASE_FUNC_ADD: {
			sollya_obj_t g_1 = NULL;
			sollya_obj_t g_2 = NULL;
			int n;
			sollya_lib_get_subfunctions(f, &n, &g_1, &g_2, NULL);
			std::vector<mpfr_ptr> g_1_coefs = get_poly_coefs(g_1);
			std::vector<mpfr_ptr> g_2_coefs = get_poly_coefs(g_2);
			return poly_add(g_1_coefs, g_2_coefs);
		}
		case SOLLYA_BASE_FUNC_MUL: {
			sollya_obj_t g_1 = NULL;
			sollya_obj_t g_2 = NULL;
			int n;
			sollya_lib_get_subfunctions(f, &n, &g_1, &g_2, NULL);
			std::vector<mpfr_ptr> g_1_coefs = get_poly_coefs(g_1);
			std::vector<mpfr_ptr> g_2_coefs = get_poly_coefs(g_2);
			return poly_mul(g_1_coefs, g_2_coefs);
		}
		default:
			throw std::invalid_argument(
					"Arithmetic tree contains unhandled operation.");
	}
}

/** Compute minimax approximation with Sollya.
 *
 * @param nr_terms
 * @param terms Terms to consider in the approximation.
 * @param f Function to approximate.
 * @param range of the approximation.
 * @param subterms Subset of terms to retrieve.
 * @return A vector of coefficients literals.
 */
static std::vector<std::string>
remez(int nr_terms, sollya_obj_t &terms, sollya_obj_t &f, sollya_obj_t &range,
	  int *subterms = nullptr) {
	// Run Remez optimization.
	sollya_obj_t poly = nr_terms ? sollya_lib_remez(f, terms, range, NULL)
								 : NULL;

	// Retrieve polynomial coefficients.
	std::vector<mpfr_ptr> result_mpfr = get_poly_coefs(poly);

	// Export string representations for the coefficients.
	std::vector<std::string> coefs;
	char tmp[512]; // 512 digits should always be enough.
	if (subterms) {
		for (int i = 0; i < nr_terms; i++) {
			mpfr_sprintf(tmp, "%.16Re", result_mpfr[subterms[i]]);
			coefs.emplace_back(std::string(tmp));
		}
	} else {
		for (auto c: result_mpfr) {
			mpfr_sprintf(tmp, "%.16Re", c);
			coefs.emplace_back(std::string(tmp));
		}
	}

	// Clear all Sollya and mpfr objects.
	for (auto &i: result_mpfr) {
		mpfr_clear(i);
		delete[] i;
	}
	sollya_lib_clear_obj(poly);

	return coefs;
}

#endif // ECEF2GEO_POLYNOMIAL_REMEZ_SOLLYA_HPP
