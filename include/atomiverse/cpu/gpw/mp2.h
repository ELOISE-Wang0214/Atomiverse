/// @file cpu/gpw/mp2.h
/// @brief Γ-point GPW MP2 correlation energy result (public value type).
///
/// User path: `GpwDft::mp2(atoms)` after a converged closed-shell `run()`
/// — the GPW analogue of the molecular `RHF::mp2(frame)`
/// (`cpu/scf/mp2.h`). The occupied–virtual MO repulsion integrals
/// \((ia|jb)\) are evaluated the Quickstep way: the pair density
/// \(\psi_i\psi_a\) is collocated on the uniform density grid, the periodic
/// FFT Poisson equation is solved per pair, and the potential is contracted
/// with \(\psi_j\psi_b\) by grid quadrature (Del Ben, Hutter, VandeVondele,
/// J. Chem. Theory Comput. 8, 4177 (2012); CP2K `&WF_CORRELATION` MP2
/// `DIRECT_CANONICAL`). Pair densities are neutral (\(S_{ia}=0\)), so the
/// G=0-dropped periodic Coulomb is divergence-free at Γ.
///
/// This header exposes only the `Mp2Result` value type; the MO-on-grid /
/// MO-ERI / doubles-sum kernels are private workspace in
/// `internal/mp2_kernels.h`.
#pragma once

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief Γ-point canonical GPW MP2 correlation energy (closed-shell RKS).
///
/// Spin-component channels `e_corr_os` / `e_corr_ss` are the **unscaled**
/// opposite-spin and same-spin pieces of the doubles sum (same convention as
/// molecular `cpu::scf::MP2Result`). The reported `e_corr` applies the
/// optional Grimme scales from the method members `mp2_os_scale` /
/// `mp2_ss_scale` (defaults 1.0 recover plain MP2).
struct Mp2Result {
    double e_corr = 0.0; ///< Scaled MP2 correlation energy [Hartree] (< 0).
    /// @brief Unscaled opposite-spin correlation [Hartree].
    double e_corr_os = 0.0;
    /// @brief Unscaled same-spin correlation [Hartree].
    double e_corr_ss = 0.0;
    double e_scf = 0.0;   ///< Converged GPW SCF total energy [Hartree].
    double e_total = 0.0; ///< e_scf + e_corr [Hartree].
    int n_occ = 0;  ///< Correlated occupied spatial orbitals.
    int n_virt = 0; ///< Virtual spatial orbitals in the sum (after any cap).
};

} // namespace gpw
} // namespace cpu
} // namespace ATV
