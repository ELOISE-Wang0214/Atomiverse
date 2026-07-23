/// @file cpu/scf/ccsd.h
/// @brief Canonical closed-shell CCSD correlation energy (public API).
///
/// User path: `RHF::ccsd(frame)` after (or driving) a closed-shell RHF
/// reference. Spin-orbital Gauss–Stanton intermediates
/// (J. Chem. Phys. **103**, 3561 (1995); Stanton et al., J. Chem. Phys.
/// **94**, 4334 (1991)). Matches PySCF `cc.CCSD` on the same RHF orbitals.
#pragma once

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Closed-shell CCSD energy on a converged RHF reference.
struct CCSDResult {
    /// @brief CCSD correlation energy E_CCSD − E_HF [Hartree] (≤ 0).
    double e_corr = 0.0;
    /// @brief Total energy E_HF + E_corr [Hartree].
    double e_tot = 0.0;
    /// @brief Jacobi micro-iterations used to converge the amplitudes.
    int iterations = 0;
    /// @brief True when max|Δt₁| and max|Δt₂| fell below the internal tol.
    bool converged = false;
};

// Pack-based entry is private — see internal/ccsd_pack.h.
// Public: MeanFieldMethod::ccsd(frame) (closed-shell RHF).

} // namespace scf
} // namespace cpu
} // namespace ATV
