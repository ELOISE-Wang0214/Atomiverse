/// @file cpu/scf/xc_libxc.h
/// @brief libxc-based XC evaluator — ground truth reference for native validation.
///
/// Provides the same interface as evaluate_xc() but delegates to libxc.
/// Gated by ATV_HAS_LIBXC. When enabled, native XC output must match
/// libxc to 1e-12 at every test point — any discrepancy is a bug in native code.
///
/// LSDA: XC_LDA_X (Slater) + XC_LDA_C_VWN (VWN5, functional id 7)
/// PBE:  XC_GGA_X_PBE + XC_GGA_C_PBE
#pragma once

#ifdef ATV_HAS_LIBXC

#include "atomiverse/cpu/scf/xc_functional.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Evaluate XC functional via libxc at a single grid point.
///
/// Same interface as evaluate_xc(). Used as ground truth for
/// validating the native implementation.
///
/// @param type  Functional type (LSDA or PBE).
/// @param rho_a Alpha density [electrons/Bohr^3].
/// @param rho_b Beta density [electrons/Bohr^3].
/// @param sigma_aa |grad rho_alpha|^2 (GGA only).
/// @param sigma_ab grad rho_alpha . grad rho_beta (GGA only).
/// @param sigma_bb |grad rho_beta|^2 (GGA only).
/// @return XCOutput with energy density and potential derivatives.
XCOutput evaluate_xc_libxc(XCType type, Real rho_a, Real rho_b,
                           Real sigma_aa = 0, Real sigma_ab = 0, Real sigma_bb = 0);

} // namespace scf
} // namespace cpu
} // namespace ATV

#endif // ATV_HAS_LIBXC
