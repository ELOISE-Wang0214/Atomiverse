/// @file core/gto.h
/// @brief Gaussian-type orbital normalization functions.
///
/// Pure math, no backend dependency. Used by both CPU and CUDA code.
/// CUDA-free header (uses only the C++ `<cmath>` header).
#pragma once

#include <cmath>

namespace ATV {

/// @brief Physical GTO normalization: N = (2*alpha/pi)^{3/4} for s-type (l=0).
/// @details The standard textbook normalization ensuring <g|g> = 1 for a
///   single primitive Gaussian exp(-alpha*r^2).
/// @param alpha Gaussian exponent.
/// @return Normalization factor.
inline double phys_norm_s(double alpha) {
    return std::pow(2.0 * alpha / M_PI, 0.75);
}

/// @brief GTO normalization factor (libcint convention, matches CINTgto_norm).
/// @details Includes the 2*sqrt(pi) spherical harmonic factor on top of phys_norm.
///   For validation against libcint only.
/// @param l Angular momentum quantum number.
/// @param alpha Gaussian exponent.
/// @return Normalization factor.
inline double gto_norm(int l, double alpha) {
    double n1 = (2*l + 2 + 1) * 0.5;
    double gi = std::exp(std::lgamma(n1)) / (2.0 * std::pow(2.0 * alpha, n1));
    return 1.0 / std::sqrt(gi);
}

/// @brief Boys function F_0(x) = integral_0^1 exp(-x*t^2) dt.
/// @param x Non-negative argument.
/// @return F_0(x).
inline double boys_f0(double x) {
    if (x < 1e-12) return 1.0;
    double sx = std::sqrt(x);
    return 0.5 * std::sqrt(M_PI / x) * std::erf(sx);
}

/// @brief General Boys function F_n(x) = integral_0^1 t^{2n} exp(-x*t^2) dt.
/// @details Uses upward recurrence from F_0:
///   F_n(x) = [(2n-1)*F_{n-1}(x) - exp(-x)] / (2x)  for x > threshold.
///   F_n(0) = 1/(2n+1).
/// @param n Order (n >= 0).
/// @param x Non-negative argument.
/// @return F_n(x).
inline double boys_fn(int n, double x) {
    if (n == 0) return boys_f0(x);
    if (x < 1e-12) return 1.0 / (2 * n + 1);
    double ex = std::exp(-x);
    double inv2x = 0.5 / x;
    double fn = boys_f0(x);
    for (int k = 1; k <= n; ++k)
        fn = ((2 * k - 1) * fn - ex) * inv2x;
    return fn;
}

} // namespace ATV
