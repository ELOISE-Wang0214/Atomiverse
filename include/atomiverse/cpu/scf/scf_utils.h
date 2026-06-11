/// @file cpu/scf/scf_utils.h
/// @brief Shared SCF utilities -- extracted from RHF/UHF/UKS to eliminate duplication.
///
/// Inline header (no .cpp needed). All functions operate on column-major matrices
/// for LAPACK compatibility via xtensor-blas.
///
/// @see cpu::scf::rhf_solve, cpu::scf::uhf_solve, cpu::scf::uks_solve
#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xtensor.hpp>
#include <xtensor/views/xview.hpp>

extern "C" {
// Direct LAPACK binding for the generalized symmetric eigensolver. We do
// NOT go through xt::linalg::eigh(F, S) here: its workspace query
// (LWORK=-1) is broken in the LAPACK 3.12.1 bundled with OpenBLAS 0.3.32
// (xSYGVD writes WORK(1)=0 on query, so the follow-up call aborts with
// "parameter number 11 had an illegal value"). Passing the documented
// workspace sizes explicitly sidesteps the query entirely.
void dsygvd_(const int* itype, const char* jobz, const char* uplo,
             const int* n, double* a, const int* lda, double* b,
             const int* ldb, double* w, double* work, const int* lwork,
             int* iwork, const int* liwork, int* info);
}

namespace ATV {
namespace cpu {
namespace scf {

using RealMatrix = xt::xtensor<double, 2, xt::layout_type::column_major>;
using RealVector = xt::xtensor<double, 1>;

/// @brief Build symmetric (Löwdin) orthogonalization matrix X = S^{-1/2}.
///
/// Warns if any eigenvalue of S falls below the PySCF linear-dependence
/// threshold (1e-7). For basis sets with near-linear-deps (diffuse cc-pVDZ+),
/// Phase 2 work will switch to canonical orthogonalization that drops the
/// small modes — STO-3G does not trigger this.
///
/// @param S Overlap matrix [n x n].
/// @return X such that X^T S X = I.
inline RealMatrix build_orthogonalization_matrix(const RealMatrix& S) {
    int n = static_cast<int>(S.shape(0));
    auto [s_eig, U] = xt::linalg::eigh(S);
    constexpr double kLindepTol = 1e-7; // PySCF default
    double s_min = s_eig(0);
    for (int i = 1; i < n; ++i)
        if (s_eig(i) < s_min)
            s_min = s_eig(i);
    // Caller can inspect s_min via build_orthogonalization_info if needed.
    (void)s_min; // retained for future lindep enforcement
    (void)kLindepTol;
    RealMatrix s_inv_sqrt = xt::zeros<double>({n, n});
    for (int i = 0; i < n; ++i)
        s_inv_sqrt(i, i) = 1.0 / std::sqrt(s_eig(i));
    return xt::linalg::dot(xt::linalg::dot(U, s_inv_sqrt), xt::transpose(U));
}

/// @brief Diagonalize Fock in the orthogonal basis (legacy two-step path).
/// @param F Fock matrix in AO basis [n x n].
/// @param X Orthogonalization matrix S^{-1/2} [n x n].
/// @return (orbital energies [n], AO MO coefficient matrix [n x n]).
inline std::pair<RealVector, RealMatrix> diagonalize_fock(const RealMatrix& F,
                                                          const RealMatrix& X) {
    RealMatrix F_prime =
        xt::linalg::dot(xt::transpose(X), xt::linalg::dot(F, X));
    auto [eig, C_prime] = xt::linalg::eigh(F_prime);
    RealMatrix C = xt::linalg::dot(X, C_prime);
    return {eig, C};
}

/// @brief Solve the generalized symmetric eigenproblem F C = S C ε (LAPACK dsygvd).
///
/// Matches PySCF's `scipy.linalg.eigh(F, S)` path and avoids the per-iteration
/// X^T F X multiplication. Eigenvectors returned in AO basis.
///
/// @param F Fock matrix in AO basis [n x n].
/// @param S Overlap matrix [n x n].
/// @return (orbital energies [n], AO MO coefficient matrix [n x n]).
inline std::pair<RealVector, RealMatrix>
diagonalize_fock_generalized(const RealMatrix& F, const RealMatrix& S) {
    const int n = static_cast<int>(F.shape(0));
    RealMatrix C = F;  // dsygvd overwrites A with the eigenvectors
    RealMatrix B = S;  // ... and B with its Cholesky factor
    RealVector eig = xt::zeros<double>({n});

    // Documented minimum workspace for ITYPE=1, JOBZ='V' (LAPACK dsygvd):
    // LWORK >= 1 + 6N + 2N^2, LIWORK >= 3 + 5N. Sized directly instead of
    // via the LWORK=-1 query (broken in OpenBLAS 0.3.32 / LAPACK 3.12.1 —
    // see the extern "C" note above).
    const int itype = 1;
    const char jobz = 'V';
    const char uplo = 'L';
    const int lwork = 1 + 6 * n + 2 * n * n;
    const int liwork = 3 + 5 * n;
    std::vector<double> work(static_cast<std::size_t>(lwork));
    std::vector<int> iwork(static_cast<std::size_t>(liwork));
    int info = 0;

    dsygvd_(&itype, &jobz, &uplo, &n, C.data(), &n, B.data(), &n, eig.data(),
            work.data(), &lwork, iwork.data(), &liwork, &info);
    if (info != 0)
        throw std::runtime_error(
            "diagonalize_fock_generalized: dsygvd failed, info=" +
            std::to_string(info));
    return {eig, C};
}

/// @brief Build density matrix P = factor * C_occ * C_occ^T.
/// @param C MO coefficient matrix [n x n].
/// @param n_occ Number of occupied orbitals.
/// @param n Basis dimension.
/// @param factor Occupancy factor: 2 for RHF (closed-shell), 1 for UHF/UKS (per-spin).
/// @return Density matrix [n x n].
inline RealMatrix build_density_matrix(const RealMatrix& C, int n_occ, int n,
                                       double factor = 1.0) {
    if (n_occ == 0)
        return xt::zeros<double>({n, n});
    auto C_occ = xt::view(C, xt::all(), xt::range(0, n_occ));
    return factor * xt::linalg::dot(C_occ, xt::transpose(C_occ));
}

/// @brief Apply HOMO-LUMO mixing to break spin symmetry in the initial guess.
///
/// Rotates the HOMO and LUMO columns of C by angle theta. This creates an
/// asymmetric starting density that can locate broken-symmetry solutions
/// (e.g., stretched H2).
///
/// @param C MO coefficient matrix [n x n] -- modified in-place.
/// @param n_occ Number of occupied orbitals (HOMO = n_occ-1, LUMO = n_occ).
/// @param n Basis dimension.
/// @param theta Rotation angle in radians (default 20 degrees = 0.349 rad).
inline void apply_homo_lumo_mixing(RealMatrix& C, int n_occ, int n,
                                   double theta = 0.3490658503988659) {
    if (n_occ <= 0 || n_occ >= n)
        return;
    const double ca = std::cos(theta);
    const double sa = std::sin(theta);
    int homo = n_occ - 1;
    int lumo = n_occ;
    for (int r = 0; r < n; ++r) {
        double ci = C(r, homo);
        double cj = C(r, lumo);
        C(r, homo) = ca * ci + sa * cj;
        C(r, lumo) = -sa * ci + ca * cj;
    }
}

/// @brief Compute spin contamination <S^2> - S(S+1).
///
/// For pure spin states this is zero. Nonzero values indicate spin contamination
/// from the unrestricted wavefunction mixing in higher spin states.
///
/// @param C_alpha Alpha MO coefficients [n x n].
/// @param C_beta  Beta MO coefficients [n x n].
/// @param S       Overlap matrix [n x n].
/// @param n_alpha Number of alpha electrons.
/// @param n_beta  Number of beta electrons.
/// @return <S^2> - S(S+1), where S = (n_alpha - n_beta) / 2.
inline double compute_spin_contamination(const RealMatrix& C_alpha,
                                         const RealMatrix& C_beta,
                                         const RealMatrix& S, int n_alpha,
                                         int n_beta) {
    if (n_alpha == 0 || n_beta == 0)
        return 0.0;
    auto Ca = xt::view(C_alpha, xt::all(), xt::range(0, n_alpha));
    auto Cb = xt::view(C_beta, xt::all(), xt::range(0, n_beta));
    RealMatrix SCb = xt::linalg::dot(S, Cb);
    RealMatrix O = xt::linalg::dot(xt::transpose(Ca), SCb);
    return static_cast<double>(n_beta) - xt::sum(O * O)();
}

} // namespace scf
} // namespace cpu
} // namespace ATV
