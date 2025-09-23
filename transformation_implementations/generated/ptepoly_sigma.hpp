/*
 * Generated. DO NOT EDIT.
 */

#ifndef ECEF2GEO_SIGMA2_HPP
#define ECEF2GEO_SIGMA2_HPP

namespace pteseries {

namespace priv {

/** Point-to-ellipse Fourier series expansions.
 *
 * Nilsson, J.-O. Point-to-ellipse Fourier series. doi:
 * https://doi.org/10.48550/arXiv.2507.08807
 *
 * @paramt L Total degree of the multiplicative polynomial correction.
 * @param d (phi-phi_c)^2.
 * @return Minimax approximation of sigma.
 */
template<int L> inline double sigma(double d) = delete;  // Only allow provided specializations.
template<> inline double sigma<0>(double d) { return 1.00000000000000000e+00; }
template<> inline double sigma<1>(double d) { return 1.00000000000000000e+00; }
template<> inline double sigma<2>(double d) { return 1.00000000000000000e+00 + -5.00000000000000000e-01 * d; }
template<> inline double sigma<3>(double d) { return 1.00000000000000000e+00 + -5.00000000000000000e-01 * d; }
template<> inline double sigma<4>(double d) { double d2 = d * d; return 1.00000000000000000e+00 + -5.00000000000000000e-01 * d + 4.16666666666666667e-02 * d2; }
template<> inline double sigma<5>(double d) { double d2 = d * d; return 1.00000000000000000e+00 + -5.00000000000000000e-01 * d + 4.16666666666666667e-02 * d2; }
template<> inline double sigma<6>(double d) { double d2 = d * d; return 1.00000000000000000e+00 + -5.00000000000000000e-01 * d + (4.16666666666666667e-02 + -1.38888888888888889e-03 * d) * d2; }
template<> inline double sigma<7>(double d) { double d2 = d * d; return 1.00000000000000000e+00 + -5.00000000000000000e-01 * d + (4.16666666666666667e-02 + -1.38888888888888889e-03 * d) * d2; }

}  // namespace pteseries
}  // namespace priv

#endif // ECEF2GEO_SIGMA2_HPP
