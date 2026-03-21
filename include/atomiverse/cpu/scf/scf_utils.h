/// @file cpu/scf/scf_utils.h
/// @brief Shared SCF utilities — extracted from RHF/UHF/UKS to eliminate duplication.
///
/// Inline header (no .cpp needed). All functions operate on column-major matrices
/// for LAPACK compatibility via xtensor-blas.
///
/// @see cpu::scf::rhf_solve, cpu::scf::uhf_solve, cpu::scf::uks_solve
#pragma once

#include <cmath>
#include <utility>

#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>
#include <xtensor-blas/xlinalg.hpp>

#include "atomiverse/core/real.h"
#include "atomiverse/core/int.h"

namespace ATV {
namespace cpu {
namespace scf {

using RealMatrix = xt::xtensor<Real, 2, xt::layout_type::column_major>;
using RealVector = xt::xtensor<Real, 1>;

/// @brief Build canonical orthogonalization matrix X = S^{-1/2}.
/// @param S Overlap matrix [n×n, column-major].
/// @return Orthogonalization matrix X such that X^T S X = I.
inline RealMatrix build_orthogonalization_matrix(const RealMatrix& S) {
    Int n = static_cast<Int>(S.shape(0));
    auto [s_eig, U] = xt::linalg::eigh(S);
    RealMatrix s_inv_sqrt = xt::zeros<Real>({n, n});
    for (Int i = 0; i < n; ++i)
        s_inv_sqrt(i, i) = Real(1) / std::sqrt(s_eig(i));
    return xt::linalg::dot(xt::linalg::dot(U, s_inv_sqrt), xt::transpose(U));
}

/// @brief Diagonalize Fock matrix in the orthogonal basis.
/// @param F Fock matrix in AO basis [n×n].
/// @param X Orthogonalization matrix S^{-1/2} [n×n].
/// @return (orbital energies [n], AO MO coefficient matrix [n×n]).
inline std::pair<RealVector, RealMatrix> diagonalize_fock(
    const RealMatrix& F, const RealMatrix& X) {
    RealMatrix F_prime = xt::linalg::dot(xt::transpose(X), xt::linalg::dot(F, X));
    auto [eig, C_prime] = xt::linalg::eigh(F_prime);
    RealMatrix C = xt::linalg::dot(X, C_prime);
    return {eig, C};
}

/// @brief Build density matrix P = factor * C_occ * C_occ^T.
/// @param C MO coefficient matrix [n×n].
/// @param n_occ Number of occupied orbitals.
/// @param n Basis dimension.
/// @param factor Occupancy factor: 2 for RHF (closed-shell), 1 for UHF/UKS (per-spin).
/// @return Density matrix [n×n].
inline RealMatrix build_density_matrix(const RealMatrix& C, int n_occ, Int n,
                                       Real factor = Real(1)) {
    if (n_occ == 0) return xt::zeros<Real>({n, n});
    auto C_occ = xt::view(C, xt::all(), xt::range(0, n_occ));
    return factor * xt::linalg::dot(C_occ, xt::transpose(C_occ));
}

/// @brief Apply HOMO-LUMO mixing to break spin symmetry in the initial guess.
///
/// Rotates the HOMO and LUMO columns of C by angle theta. This creates an
/// asymmetric starting density that can locate broken-symmetry solutions
/// (e.g., stretched H2).
///
/// @param C MO coefficient matrix [n×n] — modified in-place.
/// @param n_occ Number of occupied orbitals (HOMO = n_occ-1, LUMO = n_occ).
/// @param n Basis dimension.
/// @param theta Rotation angle in radians (default 20° = 0.349 rad).
inline void apply_homo_lumo_mixing(RealMatrix& C, int n_occ, Int n,
                                   Real theta = Real(0.3490658503988659)) {
    if (n_occ <= 0 || n_occ >= n) return;
    const Real ca = std::cos(theta);
    const Real sa = std::sin(theta);
    Int homo = n_occ - 1;
    Int lumo = n_occ;
    for (Int r = 0; r < n; ++r) {
        Real ci = C(r, homo);
        Real cj = C(r, lumo);
        C(r, homo) =  ca * ci + sa * cj;
        C(r, lumo) = -sa * ci + ca * cj;
    }
}

/// @brief Compute spin contamination ⟨S²⟩ - S(S+1).
///
/// For pure spin states this is zero. Nonzero values indicate spin contamination
/// from the unrestricted wavefunction mixing in higher spin states.
///
/// @param C_alpha Alpha MO coefficients [n×n].
/// @param C_beta  Beta MO coefficients [n×n].
/// @param S       Overlap matrix [n×n].
/// @param n_alpha Number of alpha electrons.
/// @param n_beta  Number of beta electrons.
/// @return ⟨S²⟩ - S(S+1), where S = (n_alpha - n_beta) / 2.
inline Real compute_spin_contamination(
    const RealMatrix& C_alpha, const RealMatrix& C_beta,
    const RealMatrix& S, int n_alpha, int n_beta) {
    if (n_alpha == 0 || n_beta == 0) return Real(0);
    auto Ca = xt::view(C_alpha, xt::all(), xt::range(0, n_alpha));
    auto Cb = xt::view(C_beta,  xt::all(), xt::range(0, n_beta));
    RealMatrix SCb = xt::linalg::dot(S, Cb);
    RealMatrix O   = xt::linalg::dot(xt::transpose(Ca), SCb);
    return Real(n_beta) - xt::sum(O * O)();
}

} // namespace scf
} // namespace cpu
} // namespace ATV
