// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the BSD 2-Clause License.

#ifndef ECEF2GEO_POLYNOMIAL_CONSTS_HPP
#define ECEF2GEO_POLYNOMIAL_CONSTS_HPP

#include <mpreal.h>

// Set default precision once for your app (like mp.dps but in bits).
inline void set_precision_bits(mpfr_prec_t bits) {
	mpfr::mpreal::set_default_prec(bits);
}

// WGS84 constants. Lazy evaluation such that set precision is respected.
// The literal values are defined by the standard.
inline const mpfr::mpreal& mp_a()  { static const mpfr::mpreal v = mpfr::mpreal("6378137.0");                               return v; } // Semi-major axis.
inline const mpfr::mpreal& mp_f()  { static const mpfr::mpreal v = mpfr::mpreal(1) / mpfr::mpreal("298.257223563");         return v; } // Flattening.
inline const mpfr::mpreal& mp_b()  { static const mpfr::mpreal v = mp_a() - mp_f() * mp_a();                                return v; } // Semi-minor axis.
inline const mpfr::mpreal& mp_e2() { static const mpfr::mpreal v = mpfr::mpreal(1) - (mp_b() * mp_b()) / (mp_a() * mp_a()); return v; } // Eccentricity squared.

#endif //ECEF2GEO_POLYNOMIAL_CONSTS_HPP
