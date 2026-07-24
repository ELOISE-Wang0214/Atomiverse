/// @file cpu/linalg/eigensolve.h
/// @brief Symmetric / generalized eigenproblem for the SCF Fock matrix:
///        Löwdin orthogonalization X = S^{-1/2} and the solutions of F C = S C ε
///        (both the two-step X^T F X path and the direct LAPACK dsygvd path).
///
/// Method-agnostic — shared by `cpu::scf` and `cpu::gpw`.
#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <xtensor-blas/xlinalg.hpp>

#include "atomiverse/cpu/linalg/matrix.h"

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

// Direct LAPACK binding for the real nonsymmetric eigensolver, used by the
// non-Hermitian response spectra of time-dependent Hartree-Fock (TDHF),
// the algebraic diagrammatic construction (ADC), and equation-of-motion
// coupled-cluster singles-and-doubles (EOM-CCSD). Same rationale as dsygvd_ above:
// the bundled LAPACK's LWORK=-1 workspace query is unreliable, so we size the
// workspace from the documented minimum explicitly and never issue a query.
void dgeev_(const char* jobvl, const char* jobvr, const int* n, double* a,
            const int* lda, double* wr, double* wi, double* vl, const int* ldvl,
            double* vr, const int* ldvr, double* work, const int* lwork,
            int* info);
}

namespace ATV {
namespace cpu {
namespace linalg {

/// @brief Symmetric matrix power M^p for a real symmetric matrix.
///
/// Eigendecompose M = U diag(λ) Uᵀ (LAPACK dsyevd via `xt::linalg::eigh`,
/// ascending eigenvalues), clamp each eigenvalue from below
/// (λ_eff = max(λ, floor)) so a p < 0 exponent cannot amplify a near-null
/// direction, then reassemble U diag(λ_eff^p) Uᵀ. Backs the Löwdin S^{-1/2}
/// orthogonalizer and the (A±B)^{±1/2} transforms of Casida linear-response /
/// random-phase-approximation (RPA) theory (the `floor` variant is what
/// time-dependent Hartree-Fock uses to tame near-null modes).
/// Header-inline free math with no owning type — the factory-ban exception for
/// ownerless pure math.
///
/// @param M     Real symmetric matrix, n × n (dimensionless or Hartree,
///              matrix-dependent). Read-only; the eigensolver copies internally.
/// @param p     Exponent applied to each clamped eigenvalue (dimensionless).
/// @param floor Lower clamp applied to the eigenvalues before exponentiation
///              (same units as M's eigenvalues); default 0.0 pins negative
///              eigenvalues to 0 but leaves positive ones untouched.
/// @param label Diagnostic tag folded into any error message.
/// @return M^p as a symmetric n × n matrix (units = M-eigenvalue-units^p).
/// @throws std::runtime_error when the clamp still yields a non-finite entry
///         (e.g. 0^p with p < 0 under the default floor) — no fallback value.
inline RealMatrix sym_power(const RealMatrix& M, double p, double floor = 0.0,
                            std::string_view label = "sym_power") {
    const int n = static_cast<int>(M.shape(0));
    auto [lambda, U] = xt::linalg::eigh(M); // ascending eigenvalues
    RealMatrix scaled = xt::zeros<double>({n, n});
    for (int i = 0; i < n; ++i) {
        const double lam_eff = std::max(lambda(i), floor);
        const double value = std::pow(lam_eff, p);
        if (!std::isfinite(value))
            throw std::runtime_error(
                std::string(label) + ": non-finite eigenvalue power ("
                + std::to_string(lam_eff) + "^" + std::to_string(p)
                + ") — raise floor to bound the near-null direction");
        scaled(i, i) = value;
    }
    return xt::linalg::dot(xt::linalg::dot(U, scaled), xt::transpose(U));
}

/// @brief Canonical orthogonalization of an overlap matrix with
///        eigenvalue thresholding. Direct instantiation — no factory
///        function: construction performs the eigendecomposition.
///
/// Eigendecompose S = U s Uᵀ (dimensionless overlap eigenvalues,
/// ascending), keep the directions with s_i >= tau, and form the
/// rectangular transform X = Ũ s̃^{-1/2} (n × m, m <= n) satisfying
/// XᵀSX = I_m. Near-linear-dependent directions (s_i → 0) are dropped
/// instead of amplified, which is what keeps Roothaan–Hall solvable when
/// the basis approaches linear dependence.
struct CanonicalOrth {
    RealMatrix X;      ///< n × m transform, XᵀSX = I_m (dimensionless).
    int n_dropped = 0; ///< Number of eigen-directions dropped (s_i < tau).
    double s_min = 0;  ///< Raw smallest eigenvalue of S, before filtering.
    /// @brief Löwdin S^{-1/2} (n × n, dimensionless), computed in the
    ///        constructor through the shared `sym_power(S, -0.5)` primitive and
    ///        cached here so callers that want the symmetric orthogonalizer on
    ///        a healthy basis do not recompute it. Populated only when
    ///        n_dropped == 0 (empty otherwise — a filtered basis has no
    ///        meaningful full-rank S^{-1/2}). Bit-identical to
    ///        `build_orthogonalization_matrix(S)`: both route through
    ///        `sym_power`.
    RealMatrix lowdin;

    /// @brief Eigendecompose and threshold.
    /// @param S   Symmetric overlap matrix (dimensionless), n × n. Read-only.
    /// @param tau Keep-threshold on the raw S eigenvalues (dimensionless);
    ///            s_i >= tau is kept. tau <= 0 disables filtering (all n
    ///            columns kept). PySCF >= 2.13 uses 1e-6 on S eigenvalues;
    ///            Atomiverse's SCF default is 1e-8 (mean-field lindep_threshold).
    /// @throws std::runtime_error when every column would be dropped
    ///         (no fallback values).
    CanonicalOrth(const RealMatrix& S, double tau) {
        const int n = static_cast<int>(S.shape(0));
        auto [s_eig, U] = xt::linalg::eigh(S); // ascending eigenvalues
        s_min = s_eig(0);

        int first_kept = 0;
        if (tau > 0) {
            while (first_kept < n && s_eig(first_kept) < tau)
                ++first_kept;
        }
        n_dropped = first_kept;
        const int m = n - n_dropped;
        if (m == 0)
            throw std::runtime_error(
                "CanonicalOrth: every overlap eigenvalue below tau="
                + std::to_string(tau) + " (s_max="
                + std::to_string(s_eig(n - 1)) + ") — no basis left");

        X = xt::zeros<double>({n, m});
        for (int j = 0; j < m; ++j) {
            const double inv_sqrt = 1.0 / std::sqrt(s_eig(first_kept + j));
            for (int i = 0; i < n; ++i)
                X(i, j) = U(i, first_kept + j) * inv_sqrt;
        }
        // Healthy basis: the full-rank symmetric orthogonalizer is exactly
        // S^{-1/2}, so route it through the shared `sym_power(S, -0.5)`
        // primitive (bit-identical to build_orthogonalization_matrix, which
        // uses the same call). Filtered bases keep an empty lowdin.
        if (n_dropped == 0)
            lowdin = sym_power(S, -0.5);
    }
};

/// @brief Build the symmetric (Löwdin) orthogonalization matrix X = S^{-1/2}.
/// @param S Symmetric positive-definite overlap matrix (dimensionless), n × n.
/// @return X = S^{-1/2} (n × n, dimensionless).
inline RealMatrix build_orthogonalization_matrix(const RealMatrix& S) {
    return sym_power(S, -0.5);
}

/// @brief Diagonalize Fock in the orthogonal basis (two-step path).
/// @param F Fock matrix [Hartree], n × n.
/// @param X Orthogonalization matrix (dimensionless), n × n.
/// @return (orbital energies [Hartree, n], AO MO coefficients [n x n]).
inline std::pair<RealVector, RealMatrix> diagonalize_fock(const RealMatrix& F,
                                                          const RealMatrix& X) {
    RealMatrix F_prime =
        xt::linalg::dot(xt::transpose(X), xt::linalg::dot(F, X));
    auto [eig, C_prime] = xt::linalg::eigh(F_prime);
    RealMatrix C = xt::linalg::dot(X, C_prime);
    return {eig, C};
}

/// @brief Solve the generalized symmetric eigenproblem F C = S C ε (LAPACK dsygvd).
/// @param F Fock matrix [Hartree], n × n.
/// @param S Overlap matrix (dimensionless), n × n.
/// @return (orbital energies [Hartree, n], AO MO coefficients [n x n]).
inline std::pair<RealVector, RealMatrix>
diagonalize_fock_generalized(const RealMatrix& F, const RealMatrix& S) {
    const int n = static_cast<int>(F.shape(0));
    RealMatrix C = F; // dsygvd overwrites A with the eigenvectors
    RealMatrix B = S; // ... and B with its Cholesky factor
    RealVector eig = xt::zeros<double>({n});

    // Documented minimum workspace for ITYPE=1, JOBZ='V' (LAPACK dsygvd):
    // LWORK >= 1 + 6N + 2N^2, LIWORK >= 3 + 5N (the LWORK=-1 query is broken).
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
            "diagonalize_fock_generalized: dsygvd failed, info="
            + std::to_string(info));
    return {eig, C};
}

/// @brief Real nonsymmetric eigenproblem A v = λ v (LAPACK dgeev).
///
/// Follows the explicit-workspace / throw-on-info convention of
/// `diagonalize_fock_generalized` above — no LWORK=-1 query (the bundled
/// LAPACK's workspace query is unreliable, same failure mode documented on
/// `dsygvd_`). The eigenvalues of a real matrix are real or come in conjugate
/// pairs, so they are returned split into real (`wr`) and imaginary (`wi`)
/// parts; this backs the non-Hermitian response spectra of TDHF, ADC, and
/// EOM-CCSD (acronyms expanded at the `dgeev_` binding above). Left
/// eigenvectors are never requested (JOBVL='N').
///
/// @param A     Real nonsymmetric matrix, n × n (matrix-dependent units).
///              Read-only — dgeev overwrites its input, so a copy is taken.
/// @param wr    Out: real parts of the n eigenvalues (resized to n).
/// @param wi    Out: imaginary parts of the n eigenvalues (resized to n); a
///              real eigenvalue has wi(j) == 0, a conjugate pair has
///              wi(j) == -wi(j+1).
/// @param vr    Optional out: right eigenvectors, column-major n × n (resized
///              to n × n). nullptr requests eigenvalues only (JOBVR='N'). For a
///              real eigenvalue λ_j column j is its eigenvector; for a
///              conjugate pair (j, j+1) columns j and j+1 hold the LAPACK real
///              and imaginary parts (v(j) = vr[:,j] + i·vr[:,j+1],
///              v(j+1) = vr[:,j] - i·vr[:,j+1]).
/// @param label Diagnostic tag folded into the throw message on dgeev failure.
/// @throws std::runtime_error when dgeev returns info != 0 — no fallback value.
inline void eig_nonsymmetric(const RealMatrix& A, RealVector& wr, RealVector& wi,
                             RealMatrix* vr = nullptr,
                             std::string_view label = "eig_nonsymmetric") {
    const int n = static_cast<int>(A.shape(0));
    RealMatrix a_work = A; // dgeev overwrites its input; keep the caller's A
    wr = xt::zeros<double>({n});
    wi = xt::zeros<double>({n});

    const char jobvl = 'N';
    const char jobvr = (vr != nullptr) ? 'V' : 'N';
    const int lda = n;
    const int ldvl = 1; // VL unreferenced (JOBVL='N'), LDVL >= 1
    const int ldvr = (vr != nullptr) ? n : 1;
    double vl_unused = 0.0;
    double* vr_ptr = &vl_unused; // dummy target when JOBVR='N'
    if (vr != nullptr) {
        *vr = xt::zeros<double>({n, n}); // column-major, matches dgeev's VR
        vr_ptr = vr->data();
    }

    // Documented dgeev minimum: LWORK >= 4N when eigenvectors are requested,
    // >= 3N otherwise; a fixed margin is added and no LWORK=-1 query is made.
    const int lwork = (vr != nullptr ? 4 : 3) * n + 32;
    std::vector<double> work(static_cast<std::size_t>(lwork));
    int info = 0;

    dgeev_(&jobvl, &jobvr, &n, a_work.data(), &lda, wr.data(), wi.data(),
           &vl_unused, &ldvl, vr_ptr, &ldvr, work.data(), &lwork, &info);
    if (info != 0)
        throw std::runtime_error(std::string(label) + ": dgeev failed, info="
                                 + std::to_string(info));
}

} // namespace linalg
} // namespace cpu
} // namespace ATV
