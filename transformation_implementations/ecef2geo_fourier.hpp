// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the BSD 2-Clause License.

/*
 * Polynomial ECEF coordinate (x,y,z) to geodetic coordinate
 * (latitude, longitude, altitude) or (n-vector, altitude) transformation
 * approximations.
 */

#ifndef ECEF2GEO_FOURIER_HPP_
#define ECEF2GEO_FOURIER_HPP_

// Input and output (coordinate) structs.
#include "coordinate_structs.hpp"

// Generated polynomial approximations components.
#include "generated/ptepoly_phi.hpp"
#include "generated/ptepoly_h.hpp"
#include "generated/ptepoly_cos.hpp"
#include "generated/ptepoly_sin.hpp"
#include "generated/ptepoly_sigma.hpp"
#include "generated/ptepoly_tau.hpp"

// To get WGS84 constants.
#include "ecef2geo_common.hpp"

#include <cmath>  // std::copysign, std::sqrt

namespace ecef2geo {

template<int L, int N, int M, int Nh, int Mh>
static nva ecef2nvaf(xyz ecef) {
	double &x = ecef.x;
	double &y = ecef.y;
	double &z = ecef.z;

	double p = std::sqrt(x * x + y * y + z * z);
	double rho_inv = 1. / p;
	double t = z * rho_inv;
	double t2 = t * t;
	double h = p + ecef2geo::priv::ptepoly_h<Nh,Mh>(rho_inv, t2);
	double w = ecef2geo::priv::ptepoly_phi<N,M>(rho_inv, t2);
	double de = t2 * (1.0 - t2) * w * w;
	double si = pteseries::priv::sigma<L>(de);
	double ta = pteseries::priv::tau<L>(w, de);
	double inv_p_rho = (si - t2 * ta) * rho_inv;
	double t_eta = t * (si + (1.0 - t2) * ta);
	return {{x * inv_p_rho, y * inv_p_rho, t_eta}, h};
}

}  // namespace

#endif // ECEF2GEO_FOURIER_HPP_
