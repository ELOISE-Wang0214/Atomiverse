/// @file cpu/gpw/internal/periodic.h
/// @brief Orthorhombic Γ-point lattice-image enumeration shared by the analytic
///        one-electron / GTH integrals and the grid collocation.
///
/// At the Γ point each contracted basis function is a periodic image sum
/// χ_μ(r) = Σ_T φ_μ(r − T) over lattice vectors T = (t_x a, t_y b, t_z c). The
/// analytic overlap/kinetic/pseudopotential matrices and the real-space density
/// collocation must use the *same* set of images for the invariant
/// ∫ρ dr = Tr(P S) to hold. This header centralizes the per-axis image count so
/// both sides derive it from one shared cutoff radius (set by EPS_DEFAULT).
#pragma once

#include <cmath>

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief Per-axis half-count of lattice images covering a radius `rcut`.
///
/// For an orthorhombic cell an image t_d a_d can contribute to an integral whose
/// integrand reaches `rcut` only if |t_d| a_d <= rcut, i.e. |t_d| <= rcut / a_d.
/// @param L Box lengths a,b,c [Bohr].
/// @param periodic When false, all counts are 0 (isolated / Martyna-Tuckerman).
/// @param rcut Reach of the integrand [Bohr] (e.g. twice the basis radius for a
///        Gaussian *pair*).
/// @param m Output: half-count per axis; images run t_d ∈ [-m[d], m[d]].
inline void image_counts(const double L[3], bool periodic, double rcut, int m[3]) {
    for (int d = 0; d < 3; ++d)
        m[d] = (periodic && L[d] > 0.0)
                   ? static_cast<int>(std::ceil(rcut / L[d] - 1e-12))
                   : 0;
}

} // namespace gpw
} // namespace cpu
} // namespace ATV
