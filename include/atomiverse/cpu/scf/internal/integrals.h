/// @file cpu/scf/integrals.h
/// @brief IntegralPack — molecular integrals via native C++.
///
/// All arrays are xtensor. Computes overlap (S), kinetic (T), nuclear attraction (V),
/// core Hamiltonian (Hcore = T + V), two-electron repulsion integrals (ERI),
/// and nuclear repulsion energy (V_nuc) for a given molecular geometry.
///
/// Supports elements H through O (Z=1..8) with STO-3G basis (s and sp shells).
///
/// @see cpu::scf::SCFSolver for the CPU RHF solver that consumes IntegralPack
#pragma once

#include <array>
#include <vector>
#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/dense_matrix.h"
#include "atomiverse/cpu/scf/internal/basis_workspace.h"
#include "atomiverse/cpu/scf/native_integrals.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Molecular integrals for a Hartree-Fock / Kohn-Sham calculation.
///
/// Construct from a Frame + basis set. There is no free-function factory —
/// `IntegralPack{frame, basis}` is the only entry point.
struct IntegralPack {
    xt::xtensor<double, 2> S; ///< Overlap matrix [n_basis x n_basis].
    xt::xtensor<double, 2> T; ///< Kinetic energy matrix [n_basis x n_basis].
    xt::xtensor<double, 2>
        V; ///< Nuclear attraction matrix [n_basis x n_basis].
    xt::xtensor<double, 2>
        Hcore;                  ///< Core Hamiltonian T + V [n_basis x n_basis].
    xt::xtensor<double, 4> ERI; ///< Two-electron integrals (ij|kl) [n^4].
    double V_nuc = 0;           ///< Nuclear repulsion energy [Hartree].
    int n_basis = 0;            ///< Number of basis functions.
    int n_electrons = 0;        ///< Total number of electrons.

    /// @brief Atom descriptors (Z, position in Bohr). Filled by the
    ///        Frame-based constructor.
    std::vector<native::Atom> atoms;

    /// @brief Shell descriptors with owned coefficient storage.
    /// Each shell has n_prim, alpha (pointer into coeff_storage), coeff (pointer into coeff_storage), center.
    std::vector<native::Shell> shells;

    /// @brief Backing storage for shell exponents and normalized coefficients.
    /// One inner vector per shell, length `2 * n_prim`, laid out as
    /// `[alpha_0 ... alpha_{n_prim-1}, coeff_0 ... coeff_{n_prim-1}]`.
    /// `pack.shells[i].alpha` / `.coeff` point into this storage; any code
    /// that copies `shell_data` must also rebind those pointers.
    std::vector<std::vector<double>> shell_data;

    /// @brief Empty pack (fields filled by the caller or left default).
    IntegralPack() = default;

    /// @brief Compute integrals for a molecular geometry via native C++.
    ///
    /// Reads atom positions from the `kAtom_x/y/z` columns (Angstrom) and nuclear
    /// charges from `kAtom_type`. Converts to Bohr internally.
    ///
    /// By default the two-electron ERI tensor is NOT computed (`compute_eri=false`).
    /// The CPU SCF solver uses direct SCF via ShellPairList::g_rhf,
    /// so storing the full O(n^4) ERI tensor is unnecessary for CPU calculations.
    /// Pass `compute_eri=true` when the caller explicitly needs `ERI`:
    ///   - CUDA RHFFock (uploads ERI to device)
    ///   - Tests that verify ERI values directly
    ///
    /// @param frame       CPU Frame with atom positions and types (atoms.type = Z).
    /// @param basis       Which basis set to resolve (STO3G, 6-31G, ...).
    /// @param charge      Molecular charge (0 = neutral, +1 = cation, etc.).
    /// @param compute_eri If true, compute and store the full ERI tensor [default false].
    IntegralPack(const Frame& frame, BasisSet basis, int charge = 0,
                 bool compute_eri = false);

    /// @brief Number of doubly-occupied orbitals (RHF closed-shell).
    /// @return `n_electrons / 2`.
    int n_occ() const { return n_electrons / 2; }

    /// @brief Fused Coulomb + exchange G = J - 0.5*K (RHF stored-ERI path).
    /// @details G_ij = sum_kl P_kl [(ij|kl) - 0.5*(ik|jl)].
    ///          P is the total RHF density (includes factor of 2).
    ///          Selected by `HamiltonianBuilder` when `ERI` holds a
    ///          materialized dense tensor; otherwise the direct-SCF
    ///          `ShellBlockPairList::g_rhf` path runs. Requires
    ///          `compute_eri=true` at construction.
    /// @param P_total Total RHF density matrix [n x n, unitless].
    /// @param[out] G Fused two-electron matrix [n x n, Hartree].
    void coulomb_exchange_fused(const RealMatrix& P_total, RealMatrix& G) const;

  private:
    /// @brief Overlap matrix S_ij = <φ_i|φ_j>.
    /// @param shells Resolved shell list.
    /// @return (n_basis × n_basis) symmetric matrix [unitless].
    [[nodiscard]] static xt::xtensor<double, 2>
    overlap(const ShellList& shells);

    /// @brief Kinetic-energy matrix T_ij = <φ_i|-½∇²|φ_j>.
    /// @param shells Resolved shell list.
    /// @return (n_basis × n_basis) symmetric matrix [Hartree].
    [[nodiscard]] static xt::xtensor<double, 2>
    kinetic(const ShellList& shells);

    /// @brief Nuclear-attraction matrix V_ij = -Σ_A Z_A <φ_i|1/|r - R_A||φ_j>.
    /// @param shells Resolved shell list.
    /// @param atoms  Atoms supplying nuclear charges and positions (Bohr).
    /// @return (n_basis × n_basis) symmetric matrix [Hartree].
    [[nodiscard]] static xt::xtensor<double, 2>
    nuclear_attraction(const ShellList& shells, const AtomList& atoms);

    /// @brief Dense two-electron repulsion tensor ERI(i,j,k,l) = (ij|kl).
    ///
    /// Dense ERI storage is O(N⁴) in memory and only appropriate for small
    /// systems, validation, or code paths that must hand the full tensor to
    /// another backend (e.g. CUDA RHFFock). The production CPU direct-SCF
    /// path uses `ShellBlockPairList` in internal/shell_pair.h instead.
    ///
    /// @param shells Resolved shell list.
    /// @return (n × n × n × n) dense tensor in Hartree.
    [[nodiscard]] static xt::xtensor<double, 4>
    eri_dense(const ShellList& shells);
};

} // namespace scf
} // namespace cpu
} // namespace ATV
