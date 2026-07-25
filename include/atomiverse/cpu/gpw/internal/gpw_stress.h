/// @file cpu/gpw/internal/gpw_stress.h
/// @brief GPW cell stress via central FD of SCF energy under strain (v1).
#pragma once

#include <array>
#include <vector>

#include "atomiverse/cpu/gpw/gpw_solver.h"

namespace ATV {
namespace cpu {
namespace gpw {
namespace internal {

/// @brief Strain-FD Voigt stress [Hartree/Bohr³], order xx,yy,zz,xy,xz,yz.
///
/// Orthorhombic GPW only: diagonal components from central FD of \(E\) under
/// axis-aligned strain at fixed fractional coordinates; off-diagonal are 0.
///
/// @param atoms Geometry [Bohr] + basis/PP (same as `GpwDft::run`).
/// @param method SCF configuration (`cell` must match atoms' box).
/// @param h Dimensionless strain step (default 2e-3).
/// @param axes Per-axis mask (x,y,z); a masked axis skips both FD SCFs and
///        reports exactly 0 (variable-cell relax of constrained axes).
/// @param guess_C Optional converged unstrained MO coefficients — warm start
///        for every strained SCF (cost only; values unchanged within SCF tol).
/// @return Voigt stress σ [Hartree/Bohr³].
std::array<double, 6>
gpw_stress_impl(const std::vector<GpwAtom>& atoms, const GpwDft& method,
                double h = 2e-3,
                const std::array<bool, 3>& axes = {true, true, true},
                const linalg::RealMatrix* guess_C = nullptr);

} // namespace internal
} // namespace gpw
} // namespace cpu
} // namespace ATV
