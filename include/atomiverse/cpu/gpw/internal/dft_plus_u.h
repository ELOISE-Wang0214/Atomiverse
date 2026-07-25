/// @file cpu/gpw/internal/dft_plus_u.h
/// @brief Dudarev DFT+U (PLUS_U) helpers for GPW — energy + potential.
///
/// Occupation matrix on the AO subspace of atom A with angular momentum L:
///   \( n_{mm'} = \tfrac12\big[(PS)_{mm'} + (SP)_{mm'}\big] \)  (Mulliken, symmetrized).
/// Energy (one spin channel, U in Hartree):
///   \( E_U = \tfrac{U}{2}\big[\mathrm{Tr}\,n - \mathrm{Tr}(n n)\big] \).
/// Potential added to the spin Fock matrix on that subspace:
///   \( V_{mm'} = U\big(\tfrac12\delta_{mm'} - n_{mm'}\big) \).
///
/// Multi-k (collinear): Brillouin-zone average
///   \( n=\sum_k w_k\,\tfrac12\mathrm{Herm}\big[(P_k S_k)+(S_k P_k)\big] \)
/// (real symmetric after Hermitian projection), one \(V\) for every \(H_k\).
///
/// Matches CP2K Dudarev / `U_MINUS_J` for collinear UKS (one channel at a time).
/// Not a public header — included only from the GPW solver.
#pragma once

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#ifdef ATV_HAS_MPI
#include <mpi.h>
#endif

#include "atomiverse/cpu/linalg/matrix.h"

namespace ATV {
namespace cpu {
namespace gpw {
namespace internal {

/// @brief One AO: owning atom index and angular momentum ℓ.
struct AoShellTag {
    int atom = -1;
    int l = -1;
};

/// @brief Collect global AO indices belonging to atom `iat` with ℓ = `L`.
inline std::vector<int>
ao_indices_for_shell(const std::vector<AoShellTag>& tags, int iat, int L) {
    std::vector<int> idx;
    for (int i = 0; i < static_cast<int>(tags.size()); ++i)
        if (tags[static_cast<std::size_t>(i)].atom == iat
            && tags[static_cast<std::size_t>(i)].l == L)
            idx.push_back(i);
    return idx;
}

/// @brief Symmetrized Mulliken occupation block on AO subspace `idx`.
/// @param P Spin density matrix [electrons] (nb × nb).
/// @param S Overlap (nb × nb).
/// @param idx AO indices of the correlated shell (size n).
/// @return n × n occupation matrix.
inline linalg::RealMatrix mulliken_occupation(const linalg::RealMatrix& P,
                                              const linalg::RealMatrix& S,
                                              const std::vector<int>& idx) {
    const int n = static_cast<int>(idx.size());
    linalg::RealMatrix occ = xt::zeros<double>({n, n});
    if (n == 0)
        return occ;
    // n_ij = ½ [(PS)_ij + (SP)_ij] with i,j ∈ idx
    for (int a = 0; a < n; ++a) {
        const int i = idx[static_cast<std::size_t>(a)];
        for (int b = 0; b < n; ++b) {
            const int j = idx[static_cast<std::size_t>(b)];
            double ps = 0.0;
            double sp = 0.0;
            const int nb = static_cast<int>(P.shape(0));
            for (int k = 0; k < nb; ++k) {
                ps += P(i, k) * S(k, j);
                sp += S(i, k) * P(k, j);
            }
            occ(a, b) = 0.5 * (ps + sp);
        }
    }
    return occ;
}

/// @brief Dudarev energy and Fock potential for one spin on one shell.
/// @param U_ha U_eff [Hartree].
/// @param occ Occupation matrix (n × n), modified only as const view.
/// @param[out] V_sub Potential on the subspace (n × n); empty if n==0.
/// @return Energy contribution [Hartree].
inline double dudarev_energy_potential(double U_ha,
                                       const linalg::RealMatrix& occ,
                                       linalg::RealMatrix& V_sub) {
    const int n = static_cast<int>(occ.shape(0));
    V_sub = xt::zeros<double>({n, n});
    if (n == 0 || !(U_ha > 0.0))
        return 0.0;
    double tr = 0.0;
    double tr2 = 0.0;
    for (int i = 0; i < n; ++i) {
        tr += occ(i, i);
        for (int j = 0; j < n; ++j)
            tr2 += occ(i, j) * occ(j, i);
    }
    const double e_u = 0.5 * U_ha * (tr - tr2);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            const double delta = (i == j) ? 0.5 : 0.0;
            V_sub(i, j) = U_ha * (delta - occ(i, j));
        }
    return e_u;
}

/// @brief Scatter subspace potential into full Fock matrix.
inline void scatter_potential(linalg::RealMatrix& H,
                              const std::vector<int>& idx,
                              const linalg::RealMatrix& V_sub) {
    const int n = static_cast<int>(idx.size());
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b)
            H(idx[static_cast<std::size_t>(a)],
              idx[static_cast<std::size_t>(b)]) += V_sub(a, b);
}

/// @brief Mulliken spin population m_A = Σ_{μ∈A,ν} (Pα−Pβ)_{μν} S_{νμ}.
inline std::vector<double>
mulliken_spin_populations(const linalg::RealMatrix& Pa,
                          const linalg::RealMatrix& Pb,
                          const linalg::RealMatrix& S,
                          const std::vector<AoShellTag>& tags, int n_atoms) {
    std::vector<double> m(static_cast<std::size_t>(n_atoms), 0.0);
    const int nb = static_cast<int>(Pa.shape(0));
    for (int i = 0; i < nb; ++i) {
        const int a = tags[static_cast<std::size_t>(i)].atom;
        if (a < 0 || a >= n_atoms)
            continue;
        double acc = 0.0;
        for (int j = 0; j < nb; ++j)
            acc += (Pa(i, j) - Pb(i, j)) * S(j, i);
        m[static_cast<std::size_t>(a)] += acc;
    }
    return m;
}

/// @brief Apply DFT+U for all configured (element,L,U) entries on one spin.
/// @param z_of_atom Atomic numbers per atom (length n_atoms).
/// @param U_ha_of_element Maps atomic number → (l, U_ha); U_ha==0 skips.
/// @return Total E_U for this spin [Hartree].
struct PlusUHa {
    int l = 2;
    double U_ha = 0.0;
};

inline double apply_dft_plus_u(linalg::RealMatrix& H,
                               const linalg::RealMatrix& P,
                               const linalg::RealMatrix& S,
                               const std::vector<AoShellTag>& tags,
                               const std::vector<int>& z_of_atom,
                               const std::vector<PlusUHa>& by_z_slot,
                               const std::vector<int>& z_keys) {
    double e_u = 0.0;
    const int n_atoms = static_cast<int>(z_of_atom.size());
    for (int iat = 0; iat < n_atoms; ++iat) {
        const int z = z_of_atom[static_cast<std::size_t>(iat)];
        for (std::size_t k = 0; k < z_keys.size(); ++k) {
            if (z_keys[k] != z)
                continue;
            const PlusUHa& pu = by_z_slot[k];
            if (!(pu.U_ha > 0.0))
                continue;
            const auto idx = ao_indices_for_shell(tags, iat, pu.l);
            if (idx.empty())
                continue;
            const linalg::RealMatrix occ = mulliken_occupation(P, S, idx);
            linalg::RealMatrix V;
            e_u += dudarev_energy_potential(pu.U_ha, occ, V);
            scatter_potential(H, idx, V);
        }
    }
    return e_u;
}

/// @brief Dudarev DFT+U energy only (no Fock scatter) — fixed-\(P\) stress.
///
/// Same occupation / energy formula as `apply_dft_plus_u` for one spin
/// channel. Used by the Hellmann–Feynman cell stress Lagrangian where \(P\)
/// is held fixed and \(S\) (geometry) varies.
inline double dft_plus_u_energy(const linalg::RealMatrix& P,
                                const linalg::RealMatrix& S,
                                const std::vector<AoShellTag>& tags,
                                const std::vector<int>& z_of_atom,
                                const std::vector<PlusUHa>& by_z_slot,
                                const std::vector<int>& z_keys) {
    double e_u = 0.0;
    const int n_atoms = static_cast<int>(z_of_atom.size());
    for (int iat = 0; iat < n_atoms; ++iat) {
        const int z = z_of_atom[static_cast<std::size_t>(iat)];
        for (std::size_t k = 0; k < z_keys.size(); ++k) {
            if (z_keys[k] != z)
                continue;
            const PlusUHa& pu = by_z_slot[k];
            if (!(pu.U_ha > 0.0))
                continue;
            const auto idx = ao_indices_for_shell(tags, iat, pu.l);
            if (idx.empty())
                continue;
            const linalg::RealMatrix occ = mulliken_occupation(P, S, idx);
            linalg::RealMatrix V;
            e_u += dudarev_energy_potential(pu.U_ha, occ, V);
        }
    }
    return e_u;
}

// ── Multi-k (complex Bloch) Dudarev helpers ──────────────────────────────

/// @brief Weighted Hermitian Mulliken occupation block from one \(P_k,S_k\).
///
/// Adds \(w\cdot\tfrac12\mathrm{Herm}[(PS)+(SP)]\) (real part after
/// Hermitian projection) into `occ` (n×n, pre-sized / accumulated).
inline void accumulate_mulliken_occupation_k(
    linalg::RealMatrix& occ, const linalg::CplxMatrix& P,
    const linalg::CplxMatrix& S, const std::vector<int>& idx, double w) {
    using Cplx = std::complex<double>;
    const int n = static_cast<int>(idx.size());
    if (n == 0 || !(w > 0.0))
        return;
    const int nb = static_cast<int>(P.shape(0));
    for (int a = 0; a < n; ++a) {
        const int i = idx[static_cast<std::size_t>(a)];
        for (int b = 0; b < n; ++b) {
            const int j = idx[static_cast<std::size_t>(b)];
            Cplx ps = 0.0;
            Cplx sp = 0.0;
            for (int k = 0; k < nb; ++k) {
                ps += P(i, k) * S(k, j);
                sp += S(i, k) * P(k, j);
            }
            // Hermitian projection of ½(PS+SP), keep real part (Dudarev real V).
            const Cplx m = 0.5 * (ps + sp);
            const Cplx mh = 0.5 * (m + std::conj(m)); // force Hermitian
            // For a,b we accumulate only Re; full Herm is restored by
            // symmetrizing after the k-sum (see build_dftu_potential_k).
            occ(a, b) += w * mh.real();
        }
    }
}

/// @brief Scatter a real subspace potential into a complex Fock matrix.
inline void scatter_potential_cplx(linalg::CplxMatrix& H,
                                   const std::vector<int>& idx,
                                   const linalg::RealMatrix& V_sub) {
    const int n = static_cast<int>(idx.size());
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b)
            H(idx[static_cast<std::size_t>(a)],
              idx[static_cast<std::size_t>(b)]) += V_sub(a, b);
}

/// @brief Add a full real potential matrix onto a complex \(H_k\).
inline void add_real_potential(linalg::CplxMatrix& H,
                               const linalg::RealMatrix& V) {
    const int nb = static_cast<int>(H.shape(0));
    for (int i = 0; i < nb; ++i)
        for (int j = 0; j < nb; ++j)
            if (V(i, j) != 0.0)
                H(i, j) += V(i, j);
}

/// @brief BZ-averaged Dudarev energy + full-AO real potential for one spin.
///
/// Accumulates Mulliken occupations over k indices
/// `ik = ik_begin, ik_begin+ik_stride, …`. When `allreduce_occ` is true the
/// occupation blocks are summed across ranks (MPI k-RR) **before** the
/// nonlinear Dudarev energy/potential — partial \(n\) must not enter \(E_U\).
///
/// Serial: `ik_begin=0`, `ik_stride=1`, `allreduce_occ=false`.
inline double build_dftu_potential_k(
    const std::vector<linalg::CplxMatrix>& P_k,
    const std::vector<linalg::CplxMatrix>& S_k,
    const std::vector<double>& weights, int ik_begin, int ik_stride,
    bool allreduce_occ, const std::vector<AoShellTag>& tags,
    const std::vector<int>& z_of_atom, const std::vector<PlusUHa>& by_z_slot,
    const std::vector<int>& z_keys, linalg::RealMatrix& V_full) {
    (void)allreduce_occ;
    const int nk = static_cast<int>(S_k.size());
    const int nbas = nk > 0 ? static_cast<int>(S_k[0].shape(0)) : 0;
    V_full = xt::zeros<double>({nbas, nbas});
    if (nbas == 0 || by_z_slot.empty() || nk == 0)
        return 0.0;

    double e_u = 0.0;
    const int n_atoms = static_cast<int>(z_of_atom.size());
    for (int iat = 0; iat < n_atoms; ++iat) {
        const int z = z_of_atom[static_cast<std::size_t>(iat)];
        for (std::size_t s = 0; s < z_keys.size(); ++s) {
            if (z_keys[s] != z)
                continue;
            const PlusUHa& pu = by_z_slot[s];
            if (!(pu.U_ha > 0.0))
                continue;
            const auto idx = ao_indices_for_shell(tags, iat, pu.l);
            if (idx.empty())
                continue;
            const int n = static_cast<int>(idx.size());
            linalg::RealMatrix occ = xt::zeros<double>({n, n});
            for (int ik = ik_begin; ik < nk; ik += ik_stride) {
                accumulate_mulliken_occupation_k(
                    occ, P_k[static_cast<std::size_t>(ik)],
                    S_k[static_cast<std::size_t>(ik)], idx,
                    weights[static_cast<std::size_t>(ik)]);
            }
#ifdef ATV_HAS_MPI
            if (allreduce_occ) {
                // E_U is nonlinear in n — reduce occupations, not energies.
                // Flatten through a contiguous buffer (xtensor layout-safe).
                std::vector<double> buf(static_cast<std::size_t>(n * n));
                for (int a = 0; a < n; ++a)
                    for (int b = 0; b < n; ++b)
                        buf[static_cast<std::size_t>(a * n + b)] = occ(a, b);
                MPI_Allreduce(MPI_IN_PLACE, buf.data(), n * n, MPI_DOUBLE,
                              MPI_SUM, MPI_COMM_WORLD);
                for (int a = 0; a < n; ++a)
                    for (int b = 0; b < n; ++b)
                        occ(a, b) = buf[static_cast<std::size_t>(a * n + b)];
            }
#endif
            // Symmetrize accumulated real occupation (numerical Hermiticity).
            for (int a = 0; a < n; ++a)
                for (int b = a + 1; b < n; ++b) {
                    const double avg = 0.5 * (occ(a, b) + occ(b, a));
                    occ(a, b) = avg;
                    occ(b, a) = avg;
                }
            linalg::RealMatrix V_sub;
            e_u += dudarev_energy_potential(pu.U_ha, occ, V_sub);
            scatter_potential(V_full, idx, V_sub);
        }
    }
    return e_u;
}

/// @brief Apply multi-k Dudarev to every local \(H_k\) for one spin.
/// @return \(E_U\) [Hartree] for this spin (global when `allreduce_occ`).
inline double apply_dft_plus_u_k(std::vector<linalg::CplxMatrix>& H_k,
                                 const std::vector<linalg::CplxMatrix>& P_k,
                                 const std::vector<linalg::CplxMatrix>& S_k,
                                 const std::vector<double>& weights,
                                 int ik_begin, int ik_stride, bool allreduce_occ,
                                 const std::vector<AoShellTag>& tags,
                                 const std::vector<int>& z_of_atom,
                                 const std::vector<PlusUHa>& by_z_slot,
                                 const std::vector<int>& z_keys,
                                 linalg::RealMatrix* V_out = nullptr) {
    linalg::RealMatrix V_full;
    const double e_u = build_dftu_potential_k(
        P_k, S_k, weights, ik_begin, ik_stride, allreduce_occ, tags, z_of_atom,
        by_z_slot, z_keys, V_full);
    const int nk = static_cast<int>(H_k.size());
    for (int ik = ik_begin; ik < nk; ik += ik_stride)
        add_real_potential(H_k[static_cast<std::size_t>(ik)], V_full);
    if (V_out)
        *V_out = std::move(V_full);
    return e_u;
}

} // namespace internal
} // namespace gpw
} // namespace cpu
} // namespace ATV
