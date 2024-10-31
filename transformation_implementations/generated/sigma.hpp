/*
 * Generated. DO NOT EDIT.
 */

#ifndef ECEF2GEO_SIGMA_HPP
#define ECEF2GEO_SIGMA_HPP

namespace ecef2geo {

// Minimax approximation ranges.
constexpr double SIGMA_H_MIN = -5.00000000000000000e+03;
constexpr double SIGMA_H_MAX = 1.00000000000000000e+05;
constexpr double SIGMA_DELTA_MIN = 0.00000000000000000e+00;
constexpr double SIGMA_DELTA_MAX = 1.1334967451904287e-5;

namespace priv {

/** Minimax approximations of sigma.
 *
 * For details see:
 * Nilsson, J.-O. Minimax polynomial ECEF to geodetic coordinate transformation
 * approximations. doi: https://doi.org/10.48550/arXiv.2405.06572
 *
 * @paramt L Total degree of the multiplicative polynomial correction.
 * @param d (phi-phi_c)^2.
 * @return Minimax approximation of sigma.
 */
template<int L> inline double sigma(double d) = delete;  // Only allow provided specializations.
template<> inline double sigma<0>(double d) { return 9.9999716626081372e-01; }
template<> inline double sigma<1>(double d) { return 9.9999716626081372e-01; }
template<> inline double sigma<2>(double d) { return 9.9999999999933083e-01 + -4.9999952770986795e-01 * d; }
template<> inline double sigma<3>(double d) { return 9.9999999999933083e-01 + -4.9999952770986795e-01 * d; }
template<> inline double sigma<4>(double d) { double d2 = d * d; return 1.0000000000000000e+00 + -4.9999999999989962e-01 * d + 4.1666643052156918e-02 * d2; }
template<> inline double sigma<5>(double d) { double d2 = d * d; return 1.0000000000000000e+00 + -4.9999999999989962e-01 * d + 4.1666643052156918e-02 * d2; }
template<> inline double sigma<6>(double d) { double d2 = d * d; return 1.0000000000000000e+00 + -5.0000000000000000e-01 * d + (4.1666666666662683e-02 + -1.3888883266386210e-03 * d) * d2; }
template<> inline double sigma<7>(double d) { double d2 = d * d; return 1.0000000000000000e+00 + -5.0000000000000000e-01 * d + (4.1666666666662683e-02 + -1.3888883266386210e-03 * d) * d2; }

}  // namespace ecef2geo
}  // namespace priv

#endif // ECEF2GEO_SIGMA_HPP
