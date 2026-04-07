/// @file cpu/scf/xc_functional.h
/// @brief Exchange-correlation functional evaluator for UKS.
///
/// Supports LSDA (Slater exchange + VWN5 correlation) and PBE GGA.
/// Evaluates XC energy density and potential derivatives at a single grid point.
///
/// References:
///   - Slater: Dirac, Proc. Cambridge Phil. Soc. 26, 376 (1930)
///   - VWN5: Vosko, Wilk, Nusair, Can. J. Phys. 58, 1200 (1980), parameterization V
///   - PBE: Perdew, Burke, Ernzerhof, PRL 77, 3865 (1996)
#pragma once

namespace ATV {
namespace cpu {
namespace scf {

/// @brief XC functional type selector.
enum class XCType { LSDA, PBE };

/// @brief Output from XC evaluation at a single grid point.
struct XCOutput {
    double exc;          ///< XC energy density per electron f(rho_alpha, rho_beta) [Hartree].
    double vrho_alpha;   ///< df/drho_alpha (XC potential for alpha).
    double vrho_beta;    ///< df/drho_beta (XC potential for beta).
    // GGA only:
    double vsigma_aa;    ///< df/dsigma_alpha_alpha.
    double vsigma_ab;    ///< df/dsigma_alpha_beta.
    double vsigma_bb;    ///< df/dsigma_beta_beta.
};

/// @brief Evaluate XC functional at a single grid point.
///
/// Returns the XC energy density f such that E_xc = integral f(rho) dtau,
/// along with potential derivatives for V_xc matrix construction.
///
/// Applies density cutoff: if rho_total < 1e-12, returns zero for all outputs.
///
/// @param type  Functional type (LSDA or PBE).
/// @param rho_a Alpha density at grid point [electrons/Bohr^3].
/// @param rho_b Beta density at grid point [electrons/Bohr^3].
/// @param sigma_aa |grad rho_alpha|^2 (GGA only, ignored for LSDA).
/// @param sigma_ab grad rho_alpha . grad rho_beta (GGA only).
/// @param sigma_bb |grad rho_beta|^2 (GGA only).
/// @return XCOutput with energy density and potential derivatives.
XCOutput evaluate_xc(XCType type, double rho_a, double rho_b,
                     double sigma_aa = 0, double sigma_ab = 0, double sigma_bb = 0);

/// @brief Batched XC evaluation over an array of grid points.
///
/// Processes `n` grid points at once, enabling SIMD auto-vectorization of
/// transcendental functions (pow, log, exp).
///
/// @param type     Functional type (LSDA or PBE).
/// @param n        Number of grid points.
/// @param rho_a    Alpha densities [n].
/// @param rho_b    Beta densities [n].
/// @param sigma_aa |grad rho_alpha|^2 array [n] (GGA only, may be nullptr for LDA).
/// @param sigma_ab grad rho_alpha . grad rho_beta array [n] (GGA only).
/// @param sigma_bb |grad rho_beta|^2 array [n] (GGA only).
/// @param[out] exc     XC energy density per electron [n].
/// @param[out] vrho_a  df/drho_alpha [n].
/// @param[out] vrho_b  df/drho_beta [n].
/// @param[out] vsigma_aa df/dsigma_alpha_alpha [n] (GGA only, may be nullptr for LDA).
/// @param[out] vsigma_ab df/dsigma_alpha_beta [n] (GGA only).
/// @param[out] vsigma_bb df/dsigma_beta_beta [n] (GGA only).
void evaluate_xc_batch(XCType type, int n,
                       const double* rho_a, const double* rho_b,
                       const double* sigma_aa, const double* sigma_ab,
                       const double* sigma_bb,
                       double* exc, double* vrho_a, double* vrho_b,
                       double* vsigma_aa, double* vsigma_ab,
                       double* vsigma_bb);

} // namespace scf
} // namespace cpu
} // namespace ATV
