/// @file cpu/scf/mp2.h
/// @brief AO->MO two-electron integral transformation (ovov block) and
///        closed-shell canonical MP2 correlation energy.
///
/// The first post-SCF (post-Hartree-Fock) correlation method in the Gaussian
/// SCF domain. Both entry points are free functions in `ATV::cpu::scf` that
/// consume a converged `SCFResult` and an explicitly materialized
/// `IntegralPack::ERI` -- pure functions of the SCF result, in the same family
/// as `scf_gradient`. They are NOT Driver/Component subclasses and register in
/// no Plan.
///
/// References:
///   Moller, C.; Plesset, M. S. Phys. Rev. 46, 618 (1934).
///     DOI 10.1103/PhysRev.46.618
///   Head-Gordon, M.; Pople, J. A.; Frisch, M. J. Chem. Phys. Lett. 153, 503
///     (1988). DOI 10.1016/0009-2614(88)85250-3  (quarter transformation)
#pragma once

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_types.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Transform AO two-electron integrals to the MO (ia|jb) ovov block.
///
/// Performs the four sequential O(N^5) quarter transforms
/// (mu nu|lambda sigma) -> (i nu|lambda sigma) -> (ia|lambda sigma)
/// -> (ia|j sigma) -> (ia|jb), occupied indices first to keep intermediates
/// small. Only the occupied x virtual x occupied x virtual block needed by MP2
/// is built; the full MO tensor is never formed. Each quarter transform is a
/// BLAS-backed matrix product (chemists' notation, real canonical orbitals).
///
/// @param ERI   AO integrals (mu nu|lambda sigma) [n_ao,n_ao,n_ao,n_ao], Hartree.
/// @param C     MO coefficients [n_ao, n_mo]; columns 0..n_occ-1 are occupied.
/// @param n_occ Number of doubly-occupied orbitals (dimensionless count).
/// @param n_mo  Number of molecular orbitals (dimensionless count);
///              n_vir = n_mo - n_occ.
/// @return (ia|jb) tensor [n_occ, n_vir, n_occ, n_vir], Hartree.
xt::xtensor<double, 4> ao2mo_ovov(const xt::xtensor<double, 4>& ERI,
                                  const xt::xtensor<double, 2>& C, int n_occ,
                                  int n_mo);

/// @brief Transform AO two-electron integrals to the full (pq|rs) MO tensor.
///
/// Performs the same four sequential BLAS-backed quarter transforms as
/// `ao2mo_ovov`, but with every index range widened to the full 0..n_mo span,
/// yielding the complete chemists'-notation MO tensor (pq|rs). Full CI needs all
/// MO orbital-pair categories (occupied and virtual on every index), so the full
/// tensor IS its narrowest sufficient interface. Each quarter transform is a
/// BLAS matrix product (real canonical orbitals); no scalar quadruple loop.
///
/// @param ERI  AO integrals (mu nu|lambda sigma) [n_ao,n_ao,n_ao,n_ao], Hartree.
/// @param C    MO coefficients [n_ao, n_mo]; the first n_mo columns are used.
/// @param n_mo Number of molecular orbitals to transform (dimensionless count).
/// @return (pq|rs) tensor [n_mo, n_mo, n_mo, n_mo], Hartree.
xt::xtensor<double, 4> ao2mo_full(const xt::xtensor<double, 4>& ERI,
                                  const xt::xtensor<double, 2>& C, int n_mo);

/// @brief Transform the AO core Hamiltonian to the MO basis: h_pq = C^T Hcore C.
///
/// Two BLAS gemms (`C^T (Hcore C)`); the one-electron companion to `ao2mo_full`
/// for post-SCF methods (e.g. full CI) that assemble a many-electron Hamiltonian
/// in the MO basis. Accumulation is in double precision (AccReal) throughout.
///
/// @param Hcore AO core Hamiltonian T + V [n_ao, n_ao], Hartree.
/// @param C     MO coefficients [n_ao, n_mo]; the first n_mo columns are used.
/// @param n_mo  Number of molecular orbitals to transform (dimensionless count).
/// @return MO core Hamiltonian h_pq [n_mo, n_mo], Hartree.
xt::xtensor<double, 2> mo_hcore(const xt::xtensor<double, 2>& Hcore,
                                const xt::xtensor<double, 2>& C, int n_mo);

/// @brief Closed-shell canonical MP2 energy.
struct MP2Result {
    double e_corr; ///< MP2 correlation energy [Hartree] (<= 0).
    double e_tot;  ///< Total MP2 energy = scf.energy + e_corr [Hartree].
};

/// @brief Closed-shell canonical MP2 correlation energy from an RHF reference.
///
/// Evaluates E_MP2 = sum_ij^occ sum_ab^vir
///   (ia|jb) [2(ia|jb) - (ib|ja)] / (eps_i + eps_j - eps_a - eps_b)
/// over the (ia|jb) ovov block from `ao2mo_ovov`, using the alpha channel of a
/// restricted reference (n_alpha == n_beta). All-electron (no frozen core),
/// matching the PySCF MP2 default. The n_vir == 0 boundary (e.g. He/STO-3G)
/// returns exactly {0.0, scf.energy}.
///
/// @param ints Molecular integrals; ERI must be materialized
///             (`build_integrals(..., compute_eri=true)`).
/// @param scf  Converged restricted SCFResult (supplies C_alpha,
///             eigenvalues_alpha, n_alpha, energy).
/// @return MP2 correlation energy and total energy [Hartree].
/// @throws std::invalid_argument if ERI is unmaterialized, the SCF result is
///         not converged, or the reference is unrestricted (n_alpha != n_beta).
MP2Result mp2_energy(const IntegralPack& ints, const SCFResult& scf);

} // namespace scf
} // namespace cpu
} // namespace ATV
