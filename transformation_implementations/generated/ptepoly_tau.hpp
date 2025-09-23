/*
 * Generated. DO NOT EDIT.
 */

#ifndef ECEF2GEO_TAU2_HPP
#define ECEF2GEO_TAU2_HPP

namespace pteseries {

namespace priv {

/** Point-to-ellipse Fourier series expansions.
 *
 * Nilsson, J.-O. Point-to-ellipse Fourier series. doi:
 * https://doi.org/10.48550/arXiv.2507.08807
 *
 * @paramt L Total degree of the multiplicative polynomial correction.
 * @param d (phi-phi_c)^2.
 * @return Minimax approximation of tau.
 */
template<int L> inline double tau(double w, double d) = delete;  // Only allow provided specializations.
template<> inline double tau<0>(double w, double d) { return 0.0; }
template<> inline double tau<1>(double w, double d) { return w * 1.00000000000000000e+00; }
template<> inline double tau<2>(double w, double d) { return w * 1.00000000000000000e+00; }
template<> inline double tau<3>(double w, double d) { return w * ( 1.00000000000000000e+00 + -1.66666666666666667e-01 * d ); }
template<> inline double tau<4>(double w, double d) { return w * ( 1.00000000000000000e+00 + -1.66666666666666667e-01 * d ); }
template<> inline double tau<5>(double w, double d) { double d2 = d * d; return w * ( 1.00000000000000000e+00 + -1.66666666666666667e-01 * d + 8.33333333333333333e-03 * d2 ); }
template<> inline double tau<6>(double w, double d) { double d2 = d * d; return w * ( 1.00000000000000000e+00 + -1.66666666666666667e-01 * d + 8.33333333333333333e-03 * d2 ); }
template<> inline double tau<7>(double w, double d) { double d2 = d * d; return w * ( 1.00000000000000000e+00 + -1.66666666666666667e-01 * d + (8.33333333333333333e-03 + -1.98412698412698413e-04 * d) * d2 ); }

}  // namespace pteseries
}  // namespace priv

#endif // ECEF2GEO_TAU2_HPP
