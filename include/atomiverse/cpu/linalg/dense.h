/// @file cpu/linalg/dense.h
/// @brief Method-agnostic dense linear algebra for AO-basis self-consistent
///        field: Löwdin orthogonalization, the generalized symmetric eigenproblem
///        F C = S C ε, and the density matrix from occupied orbitals.
///
/// These primitives belong to no single method — they are shared by the
/// Gaussian SCF (`cpu::scf`) and the GPW path (`cpu::gpw`), so they live in this
/// top-level `cpu::linalg` module rather than inside any one solver. (Lifted out
/// of the former `cpu/scf/scf_utils.h`.)
///
/// Inline header (no .cpp). All matrices are column-major for LAPACK via
/// xtensor-blas.
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
// Direct LAPACK binding for the generalized symmetric eigensolver. We do NOT go
// through xt::linalg::eigh(F, S): its workspace query (LWORK=-1) is broken in the
// LAPACK 3.12.1 bundled with OpenBLAS 0.3.32 (xSYGVD writes WORK(1)=0 on query,
// so the follow-up call aborts with "parameter number 11 had an illegal value").
// Passing the documented workspace sizes explicitly sidesteps the query entirely.
void dsygvd_(const int* itype, const char* jobz, const char* uplo, const int* n,
             double* a, const int* lda, double* b, const int* ldb, double* w,
             double* work, const int* lwork, int* iwork, const int* liwork,
             int* info);
}

namespace ATV {
namespace cpu {
namespace linalg {

using RealMatrix = xt::xtensor<double, 2, xt::layout_type::column_major>;
using RealVector = xt::xtensor<double, 1>;

/// @brief Build the symmetric (Löwdin) orthogonalization matrix X = S^{-1/2}.
/// @param S Overlap matrix [n x n].
/// @return X such that X^T S X = I.
inline RealMatrix build_orthogonalization_matrix(const RealMatrix& S) {
    int n = static_cast<int>(S.shape(0));
    auto [s_eig, U] = xt::linalg::eigh(S);
    RealMatrix s_inv_sqrt = xt::zeros<double>({n, n});
    for (int i = 0; i < n; ++i) s_inv_sqrt(i, i) = 1.0 / std::sqrt(s_eig(i));
    return xt::linalg::dot(xt::linalg::dot(U, s_inv_sqrt), xt::transpose(U));
}

/// @brief Diagonalize Fock in the orthogonal basis (two-step path).
/// @param F Fock matrix in AO basis [n x n].
/// @param X Orthogonalization matrix S^{-1/2} [n x n].
/// @return (orbital energies [n], AO MO coefficient matrix [n x n]).
inline std::pair<RealVector, RealMatrix> diagonalize_fock(const RealMatrix& F,
                                                          const RealMatrix& X) {
    RealMatrix F_prime = xt::linalg::dot(xt::transpose(X), xt::linalg::dot(F, X));
    auto [eig, C_prime] = xt::linalg::eigh(F_prime);
    RealMatrix C = xt::linalg::dot(X, C_prime);
    return {eig, C};
}

/// @brief Solve the generalized symmetric eigenproblem F C = S C ε (LAPACK dsygvd).
/// @param F Fock matrix in AO basis [n x n].
/// @param S Overlap matrix [n x n].
/// @return (orbital energies [n], AO MO coefficient matrix [n x n]).
inline std::pair<RealVector, RealMatrix>
diagonalize_fock_generalized(const RealMatrix& F, const RealMatrix& S) {
    const int n = static_cast<int>(F.shape(0));
    RealMatrix C = F; // dsygvd overwrites A with the eigenvectors
    RealMatrix B = S; // ... and B with its Cholesky factor
    RealVector eig = xt::zeros<double>({n});

    // Documented minimum workspace for ITYPE=1, JOBZ='V' (LAPACK dsygvd):
    // LWORK >= 1 + 6N + 2N^2, LIWORK >= 3 + 5N. Sized directly (the LWORK=-1
    // query is broken — see the extern "C" note above).
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
        throw std::runtime_error("diagonalize_fock_generalized: dsygvd failed, info=" +
                                 std::to_string(info));
    return {eig, C};
}

/// @brief Build density matrix P = factor * C_occ * C_occ^T.
/// @param C MO coefficient matrix [n x n].
/// @param n_occ Number of occupied orbitals.
/// @param n Basis dimension.
/// @param factor Occupancy factor: 2 for RHF (closed-shell), 1 per spin (UHF/UKS).
inline RealMatrix build_density_matrix(const RealMatrix& C, int n_occ, int n,
                                       double factor = 1.0) {
    if (n_occ == 0) return xt::zeros<double>({n, n});
    auto C_occ = xt::view(C, xt::all(), xt::range(0, n_occ));
    return factor * xt::linalg::dot(C_occ, xt::transpose(C_occ));
}

} // namespace linalg
} // namespace cpu
} // namespace ATV
