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

#include "generated/phi.hpp"
#include "generated/h.hpp"
#include "generated/sigma.hpp"
#include "generated/tau.hpp"

// To get WGS84 constants.
#include "ecef2geo_common.hpp"

#include <cmath>  // std::copysign, std::sqrt

namespace ecef2geo {

/** ECEF to n-vector and altitude via the point-to-ellipse Fourier series.
 *
 * Transforms ECEF coordinates (x, y, z) to the geodetic n-vector and
 * ellipsoidal altitude. The point-to-ellipse Fourier series is evaluated
 * directly as truncated power-series polynomials (the ptepoly_* components):
 *
 * Nilsson, J.-O. Point-to-ellipse Fourier series. doi:
 * https://doi.org/10.48550/arXiv.2507.08807
 *
 * The geocentric radius p = |ecef| gives the altitude (h = p + a height
 * correction), and t = z / p (the sine of the geocentric latitude) drives the
 * latitude-dependent series. The template orders trade accuracy for speed.
 *
 * @paramt L Degree of the sigma/tau correction polynomials.
 * @paramt N Number of Fourier coefficients of the latitude (phi) series.
 * @paramt M Order of the radial (altitude) dependency of the phi series.
 * @paramt Nh Number of Fourier coefficients of the height (h) series.
 * @paramt Mh Order of the radial dependency of the height series.
 * @param ecef {x, y, z}
 * @return {n-vector, altitude}
 */
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

/** ECEF to n-vector and altitude via minimax point-to-ellipse approximations.
 *
 * Same transformation and algorithm as ecef2nvaf(), but the point-to-ellipse
 * Fourier series coefficients and the sigma/tau corrections are replaced by
 * minimax polynomial approximations (the phi/h and sigma/tau components),
 * combining the point-to-ellipse Fourier series
 *
 * Nilsson, J.-O. Point-to-ellipse Fourier series. doi:
 * https://doi.org/10.48550/arXiv.2507.08807
 *
 * with the minimax polynomial approach of
 *
 * Nilsson, J.-O. Minimax polynomial ECEF to geodetic coordinate transformation
 * approximations. doi: https://doi.org/10.48550/arXiv.2405.06572
 *
 * For the same template orders this generally gives a smaller worst-case error
 * over the supported altitude range than the truncated series of ecef2nvaf().
 *
 * @paramt L Degree of the sigma/tau minimax correction polynomials.
 * @paramt N Number of Fourier coefficients of the latitude (phi) series.
 * @paramt M Order of the radial (altitude) dependency of the phi series.
 * @paramt Nh Number of Fourier coefficients of the height (h) series.
 * @paramt Mh Order of the radial dependency of the height series.
 * @param ecef {x, y, z}
 * @return {n-vector, altitude}
 */
template<int L, int N, int M, int Nh, int Mh>
static nva ecef2nvaf2(xyz ecef) {
	double &x = ecef.x;
	double &y = ecef.y;
	double &z = ecef.z;

	double p = std::sqrt(x * x + y * y + z * z);
	double rho_inv = 1. / p;
	double t = z * rho_inv;
	double t2 = t * t;
	double h = p + ecef2geo::priv::h<Nh,Mh>(rho_inv, t2);
	double w = ecef2geo::priv::phi<N,M>(rho_inv, t2);
	double de = t2 * (1.0 - t2) * w * w;
	double si = priv::sigma<L>(de);  // Minimax versions.
	double ta = priv::tau<L>(w, de);  // Minimax versions.
	double inv_p_rho = (si - t2 * ta) * rho_inv;
	double t_eta = t * (si + (1.0 - t2) * ta);
	return {{x * inv_p_rho, y * inv_p_rho, t_eta}, h};
}

}  // namespace

#endif // ECEF2GEO_FOURIER_HPP_
