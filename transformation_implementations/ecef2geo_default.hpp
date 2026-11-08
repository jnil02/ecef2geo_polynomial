// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the BSD 2-Clause License.

/*
 * Convenience header that bundles the checked-in, common-range coefficient set
 * with the transformation headers, so the repository can be used without
 * running the coefficient optimizer.
 *
 * Note, this is NOT a universal approximation. The checked-in coefficients are
 * optimized for a specific "common" altitude range (currently [-5000, 100000]
 * m). Using them outside that range degrades accuracy. For any other range,
 * generate a dedicated coefficient set (see coefficient_computations and the
 * generate_coefficients CMake target) and wire it up the same way. For an
 * example, see the neighboring ecef2geo_benchmark project and specifically the
 * file nilsson.hpp. Choosing this default set is therefore a deliberate,
 * explicit choice. That is why the transformation headers do not bundle it.
 */

#ifndef ECEF2GEO_DEFAULT_HPP_
#define ECEF2GEO_DEFAULT_HPP_

// Range-specific generated coefficients (polynomial + trigonometric).
#include "generated/mu.hpp"
#include "generated/omega.hpp"
#include "generated/sigma.hpp"
#include "generated/tau.hpp"
#include "generated/xi.hpp"
#include "generated/chi.hpp"

// Generated coefficients for the Fourier (point-to-ellipse) approximations.
#include "generated/ptepoly_phi.hpp"
#include "generated/ptepoly_h.hpp"
#include "generated/ptepoly_cos.hpp"
#include "generated/ptepoly_sin.hpp"
#include "generated/ptepoly_sigma.hpp"
#include "generated/ptepoly_tau.hpp"
#include "generated/phi.hpp"
#include "generated/h.hpp"

// Transformations (coefficient-agnostic, requiring the coefficients above to
// be included first, which they are).
#include "ecef2geo_trigonometric.hpp"
#include "ecef2geo_polynomial.hpp"
#include "ecef2geo_fourier.hpp"

#endif  // ECEF2GEO_DEFAULT_HPP_
