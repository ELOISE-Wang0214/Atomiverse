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

#include "atomiverse/core/real.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief XC functional type selector.
enum class XCType { LSDA, PBE };

/// @brief Output from XC evaluation at a single grid point.
struct XCOutput {
    Real exc;          ///< XC energy density per electron f(ρ_α, ρ_β) [Hartree].
    Real vrho_alpha;   ///< ∂f/∂ρ_α (XC potential for alpha).
    Real vrho_beta;    ///< ∂f/∂ρ_β (XC potential for beta).
    // GGA only:
    Real vsigma_aa;    ///< ∂f/∂σ_αα.
    Real vsigma_ab;    ///< ∂f/∂σ_αβ.
    Real vsigma_bb;    ///< ∂f/∂σ_ββ.
};

/// @brief Evaluate XC functional at a single grid point.
///
/// Returns the XC energy density f such that E_xc = ∫ f(ρ) dτ,
/// along with potential derivatives for V_xc matrix construction.
///
/// Applies density cutoff: if ρ_total < 1e-12, returns zero for all outputs.
///
/// @param type  Functional type (LSDA or PBE).
/// @param rho_a Alpha density at grid point [electrons/Bohr³].
/// @param rho_b Beta density at grid point [electrons/Bohr³].
/// @param sigma_aa |∇ρ_α|² (GGA only, ignored for LSDA).
/// @param sigma_ab ∇ρ_α · ∇ρ_β (GGA only).
/// @param sigma_bb |∇ρ_β|² (GGA only).
/// @return XCOutput with energy density and potential derivatives.
XCOutput evaluate_xc(XCType type, Real rho_a, Real rho_b,
                     Real sigma_aa = 0, Real sigma_ab = 0, Real sigma_bb = 0);

/// @brief Batched XC evaluation over an array of grid points.
///
/// Processes `n` grid points at once, enabling SIMD auto-vectorization of
/// transcendental functions (pow, log, exp).
///
/// @param type     Functional type (LSDA or PBE).
/// @param n        Number of grid points.
/// @param rho_a    Alpha densities [n].
/// @param rho_b    Beta densities [n].
/// @param sigma_aa |∇ρ_α|² array [n] (GGA only, may be nullptr for LDA).
/// @param sigma_ab ∇ρ_α · ∇ρ_β array [n] (GGA only).
/// @param sigma_bb |∇ρ_β|² array [n] (GGA only).
/// @param[out] exc     XC energy density per electron [n].
/// @param[out] vrho_a  ∂f/∂ρ_α [n].
/// @param[out] vrho_b  ∂f/∂ρ_β [n].
/// @param[out] vsigma_aa ∂f/∂σ_αα [n] (GGA only, may be nullptr for LDA).
/// @param[out] vsigma_ab ∂f/∂σ_αβ [n] (GGA only).
/// @param[out] vsigma_bb ∂f/∂σ_ββ [n] (GGA only).
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
