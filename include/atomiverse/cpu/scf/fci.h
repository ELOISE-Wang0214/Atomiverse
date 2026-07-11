/// @file cpu/scf/fci.h
/// @brief Tiny determinant-basis full configuration interaction (full CI)
///        reference solver for closed-shell RHF references.
///
/// The exact correlation energy in a finite basis, for MILLISECOND-scale tiny
/// systems (H2/STO-3G, He/STO-3G). Like `scf_gradient` and `mp2_energy`, the
/// entry points are free functions in `ATV::cpu::scf` that consume a converged
/// restricted `SCFResult` plus an explicitly materialized `IntegralPack::ERI`
/// -- pure functions of the SCF result. They are NOT Driver/Component
/// subclasses and register in no Plan; all determinant enumeration,
/// Slater-Condon assembly, and Jordan-Wigner phase logic is private to fci.cpp.
///
/// This is a reference/validation tool, not a production correlation method:
/// the dense N_det x N_det CI matrix and its `xt::linalg::eigh` diagonalization
/// grow factorially, so `fci_energy` throws for anything past a hard TINY cap
/// (n_mo > 14 or N_det > 2000). For real correlation at scale use MP2
/// (`mp2_energy`) or DFT.
///
/// References:
///   Slater, J. C. Phys. Rev. 34, 1293 (1929). DOI 10.1103/PhysRev.34.1293
///   Condon, E. U. Phys. Rev. 36, 1121 (1930). DOI 10.1103/PhysRev.36.1121
///     (Slater-Condon rules)
///   Knowles, P. J.; Handy, N. C. Chem. Phys. Lett. 111, 315 (1984).
///     DOI 10.1016/0009-2614(84)85513-X  (determinant-basis FCI)
///   Szabo, A.; Ostlund, N. S. Modern Quantum Chemistry, sec. 2.3 & 4.3.
#pragma once

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_solve.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Full-CI energy result.
struct FCIResult {
    double e_fci;  ///< Exact (full-CI) total energy [Hartree].
    double e_corr; ///< Correlation energy e_fci - scf.energy [Hartree] (<= 0).
};

/// @brief Build the dense determinant-basis full-CI Hamiltonian H_CI.
///
/// Enumerates all alpha/beta occupation bitstrings (canonical block order
/// 1a..n_mo a, 1b..n_mo b) and fills the symmetric N_det x N_det matrix via the
/// Slater-Condon rules (Case A diagonal, Case B single, Case C double, Case D
/// zero) in chemists'-notation MO integrals, with Jordan-Wigner parity phases.
/// Each diagonal element carries the nuclear repulsion `ints.V_nuc` once, so the
/// lowest eigenvalue is a total energy. Exposed alongside `fci_energy` so tests
/// can assert Hermiticity and the determinant count directly (mirrors how
/// `ao2mo_ovov` is exposed next to `mp2_energy`).
///
/// @param ints Molecular integrals; ERI must be materialized
///             (`build_integrals(..., compute_eri=true)`). Supplies Hcore, ERI,
///             V_nuc [Hartree].
/// @param scf  Converged restricted SCFResult (supplies C_alpha, n_alpha,
///             n_beta, energy). n_mo = C_alpha.shape(1).
/// @return Symmetric CI Hamiltonian [N_det x N_det], Hartree.
/// @throws std::invalid_argument if ERI is unmaterialized, the SCF result is not
///         converged, the reference is open-shell (n_alpha != n_beta), or the
///         system exceeds the TINY cap (n_mo > 14 or N_det > 2000).
xt::xtensor<double, 2> fci_hamiltonian(const IntegralPack& ints,
                                       const SCFResult& scf);

/// @brief Exact closed-shell full-CI energy from an RHF reference.
///
/// Diagonalizes the dense CI Hamiltonian from `fci_hamiltonian` and takes the
/// lowest eigenvalue: `auto [w, U] = xt::linalg::eigh(H_CI); e_fci = w(0)`. The
/// single-determinant boundary (N_det == 1, e.g. He/STO-3G with no virtuals)
/// short-circuits to exactly {scf.energy, 0.0}: full CI is identical to the RHF
/// reference and the correlation energy is exactly zero.
///
/// @param ints Molecular integrals; ERI must be materialized
///             (`build_integrals(..., compute_eri=true)`). Supplies Hcore, ERI,
///             V_nuc [Hartree].
/// @param scf  Converged restricted SCFResult (supplies C_alpha, n_alpha,
///             n_beta, energy [Hartree]).
/// @return Full-CI total energy and correlation energy [Hartree].
/// @throws std::invalid_argument if ERI is unmaterialized, the SCF result is not
///         converged, the reference is open-shell (n_alpha != n_beta), or the
///         system exceeds the TINY cap (n_mo > 14 or N_det > 2000).
FCIResult fci_energy(const IntegralPack& ints, const SCFResult& scf);

} // namespace scf
} // namespace cpu
} // namespace ATV
