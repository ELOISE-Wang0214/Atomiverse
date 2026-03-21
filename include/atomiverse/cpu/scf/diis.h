/// @file cpu/scf/diis.h
/// @brief Combined DIIS accelerator for spin-resolved SCF.
///
/// Pulay's Direct Inversion of the Iterative Subspace (DIIS) with combined
/// error metric across alpha/beta spin channels. Shared extrapolation
/// coefficients ensure consistent convergence for both spins.
///
/// References:
///   - Pulay, Chem. Phys. Lett. 73, 393 (1980)
///   - Pulay, J. Comp. Chem. 3, 556 (1982)
#pragma once

#include <deque>
#include <utility>

#include <xtensor/xtensor.hpp>
#include <xtensor-blas/xlinalg.hpp>

#include "atomiverse/core/real.h"
#include "atomiverse/core/int.h"
#include "atomiverse/cpu/scf/scf_utils.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Compute DIIS error vector: e = F P S - S P F (commutator in AO basis).
/// @param F Fock matrix [n×n].
/// @param P Density matrix [n×n].
/// @param S Overlap matrix [n×n].
/// @return DIIS error matrix [n×n]. At self-consistency, e → 0.
inline RealMatrix diis_error(const RealMatrix& F, const RealMatrix& P,
                             const RealMatrix& S) {
    return xt::linalg::dot(F, xt::linalg::dot(P, S))
         - xt::linalg::dot(S, xt::linalg::dot(P, F));
}

/// @brief Combined DIIS accelerator for spin-resolved Fock matrices.
///
/// Stores history of Fock matrices and DIIS error vectors for both spin channels.
/// Solves the Pulay B-matrix equations with combined error metric:
///   B_ij = Tr(e^α_i · e^α_j) + Tr(e^β_i · e^β_j)
///
/// For restricted (single-spin), pass the same F/e for both channels or use
/// the single-channel overload. The combined B matrix naturally reduces to
/// standard single-channel DIIS when e_alpha == e_beta.
class DIIS {
public:
    /// @param max_vectors Maximum number of stored (F, e) pairs.
    explicit DIIS(int max_vectors = 6) : max_vectors_(max_vectors) {}

    /// @brief Store a new Fock/error pair and return extrapolated Fock matrices.
    ///
    /// @param F_alpha Alpha Fock matrix.
    /// @param e_alpha Alpha DIIS error (from diis_error()).
    /// @param F_beta  Beta Fock matrix.
    /// @param e_beta  Beta DIIS error.
    /// @return (F_alpha_extrapolated, F_beta_extrapolated).
    std::pair<RealMatrix, RealMatrix> extrapolate(
        const RealMatrix& F_alpha, const RealMatrix& e_alpha,
        const RealMatrix& F_beta,  const RealMatrix& e_beta);

    /// @brief Single-channel convenience (for restricted SCF).
    RealMatrix extrapolate(const RealMatrix& F, const RealMatrix& e);

    /// @brief Number of stored Fock/error pairs.
    int size() const { return static_cast<int>(Fa_hist_.size()); }

    /// @brief Clear all history.
    void reset() {
        Fa_hist_.clear(); Fb_hist_.clear();
        ea_hist_.clear(); eb_hist_.clear();
    }

private:
    int max_vectors_;
    std::deque<RealMatrix> Fa_hist_, Fb_hist_;
    std::deque<RealMatrix> ea_hist_, eb_hist_;
};

} // namespace scf
} // namespace cpu
} // namespace ATV
