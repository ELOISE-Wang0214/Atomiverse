/// @file cpu/gpw/internal/vdw_df.h
/// @brief Nonlocal vdW-DF correlation (DRSLL / DF1) for GPW grids.
///
/// Román–Pérez–Soler FFT evaluation of the Dion et al. nonlocal correlation
///
///   E_c^{nl} = ½ ∬ ρ(r) φ(r,r') ρ(r') dr dr'
///
/// with q₀(ρ,|∇ρ|²) and the community RPS kernel tables (DF1 / DRSLL:
/// Z_ab = −0.8491, q_cut = 5). Matches CP2K `VDW_POTENTIAL NON_LOCAL TYPE
/// DRSLL` and libvdwxc `VDWXC_DF1` energy on the same density/gradient.
///
/// Product path is **post-SCF energy-only** on the final total density
/// (non-self-consistent); potential back-propagation is a follow-on.
#pragma once

#include <vector>

#include "atomiverse/cpu/gpw/internal/collocate.h"

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief Result of a vdW-DF (DRSLL) nonlocal evaluation.
struct VdwDfResult {
    double energy = 0.0; ///< E_c^{nl} [Hartree].
};

/// @brief DF1 / DRSLL Z_ab gradient-correction prefactor (Dion et al.).
constexpr double kVdwDf1Zab = -0.8491;
/// @brief Soft q-cutoff [1/Bohr] matching libvdwxc / QE.
constexpr double kVdwDfQCut = 5.0;

/// @brief Compute q₀(ρ, σ=|∇ρ|²) for DF1 (DRSLL) at every grid point.
///
/// @param rho   Total density [e/Bohr³], length N.
/// @param sigma |∇ρ|² [e²/Bohr⁸], length N.
/// @param q0    Output q₀ [1/Bohr], length N (resized).
void vdw_df_q0(const std::vector<double>& rho, const std::vector<double>& sigma,
               std::vector<double>& q0);

/// @brief Nonlocal vdW-DF (DRSLL) energy on a periodic orthorhombic grid.
///
/// @param grid  Uniform orthorhombic cell (Bohr).
/// @param rho   Total density on the grid, length grid.npts().
/// @param sigma |∇ρ|² on the grid (spectral or finite-difference), same length.
/// @return E_c^{nl} [Hartree].
VdwDfResult vdw_df_drsll(const UniformGrid& grid, const std::vector<double>& rho,
                         const std::vector<double>& sigma);

} // namespace gpw
} // namespace cpu
} // namespace ATV
