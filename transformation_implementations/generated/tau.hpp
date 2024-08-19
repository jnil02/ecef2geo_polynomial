/*
 * Generated. DO NOT EDIT.
 */

#ifndef ECEF2GEO_TAU_HPP
#define ECEF2GEO_TAU_HPP

namespace ecef2geo {

// Minimax approximation ranges.
constexpr double TAU_H_MIN = -5.00000000000000000e+03;
constexpr double TAU_H_MAX = 1.00000000000000000e+05;
constexpr double TAU_DELTA_MIN = 0.00000000000000000e+00;
constexpr double TAU_DELTA_MAX = 1.13349674431948401e-05;

namespace priv {

/** Minimax approximations of tau.
 *
 * Nilsson, J.-O. Minimax polynomial ECEF to geodetic coordinate transformation
 * approximations. doi: https://doi.org/10.48550/arXiv.2405.06572
 *
 * @paramt L Total degree of the multiplicative polynomial correction.
 * @param d (phi-phi_c)^2.
 * @return Minimax approximation of tau.
 */
template<int L> inline double tau(double w, double d) = delete;  // Only allow provided specializations.
template<> inline double tau<0>(double w, double d) { return 0.0; }
template<> inline double tau<1>(double w, double d) { return w * 9.9999905541991507e-01; }
template<> inline double tau<2>(double w, double d) { return w * 9.9999905541991507e-01; }
template<> inline double tau<3>(double w, double d) { return w * ( 9.9999999999986617e-01 + -1.6666657220863013e-01 * d ); }
template<> inline double tau<4>(double w, double d) { return w * ( 9.9999999999986617e-01 + -1.6666657220863013e-01 * d ); }
template<> inline double tau<5>(double w, double d) { double d2 = d * d; return w * ( 1.0000000000000000e+00 + -1.6666666666665233e-01 * d + 8.3333299598317599e-03 * d2 ); }
template<> inline double tau<6>(double w, double d) { double d2 = d * d; return w * ( 1.0000000000000000e+00 + -1.6666666666665233e-01 * d + 8.3333299598317599e-03 * d2 ); }
template<> inline double tau<7>(double w, double d) { double d2 = d * d; return w * ( 1.0000000000000000e+00 + -1.6666666666666667e-01 * d + (8.3333333333328908e-03 + -1.9841263594044442e-04 * d) * d2 ); }

}  // namespace ecef2geo
}  // namespace priv

#endif // ECEF2GEO_TAU_HPP
