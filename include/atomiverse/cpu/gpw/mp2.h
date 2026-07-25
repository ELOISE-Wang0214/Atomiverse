/// @file cpu/gpw/mp2.h
/// @brief Γ-point GPW MP2 correlation energy result (public value type).
///
/// User path: `GpwDft::mp2(atoms)` after a converged closed-shell `run()`
/// — the GPW analogue of the molecular `RHF::mp2(frame)`
/// (`cpu/scf/mp2.h`). Product default (`mp2_ri = true`) is the RI-GPW
/// plane-wave resolution-of-identity assembly of \((ia|jb)\): FFT each
/// pair density \(\rho_{ia}=\psi_i\psi_a\), store
/// \(B_{ia}(\mathbf G)=\sqrt{4\pi/|\mathbf G|^2}\,\tilde\rho_{ia}\), contract
/// \((ia|jb)=(w/N)\sum_{\mathbf G} B_{ia}^* B_{jb}\) (Del Ben, Hutter,
/// VandeVondele, J. Chem. Theory Comput. 8, 4177 (2012); CP2K
/// `&WF_CORRELATION` `METHOD MP2_GPW`). The foothold Poisson real-space path
/// remains available with `mp2_ri = false` (`DIRECT_CANONICAL`-class).
/// Pair densities are neutral (\(S_{ia}=0\)), so the G=0-dropped periodic
/// Coulomb is divergence-free at Γ.
///
/// This header exposes only the `Mp2Result` value type; the MO-on-grid /
/// MO-ERI / doubles-sum kernels are private workspace in
/// `internal/mp2_kernels.h`.
#pragma once

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief Γ-point GPW MP2 correlation energy (closed-shell RKS; RI-GPW or
///        canonical Poisson assembly).
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
    int n_occ = 0;        ///< Correlated occupied spatial orbitals.
    int n_virt = 0; ///< Virtual spatial orbitals in the sum (after any cap).
    /// @brief true when (ia|jb) used the RI-GPW G-space product path.
    bool ri = false;
};

} // namespace gpw
} // namespace cpu
} // namespace ATV
