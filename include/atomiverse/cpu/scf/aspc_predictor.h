/// @file cpu/scf/aspc_predictor.h
/// @brief Always-Stable Predictor-Corrector (ASPC) for density matrix extrapolation.
///
/// Implements the Kolafa predictor-corrector scheme for second-generation
/// Car-Parrinello AIMD. Extrapolates density matrices from previous MD steps
/// to provide a good initial guess, reducing the number of SCF iterations
/// needed per step.
///
/// Real \(\Gamma\) path: `ASPCPredictor` on real AO densities.
/// Multi-k path: `ASPCPredictorCplx` on Hermitian Bloch \(P_k\) (same real
/// \(\omega_k\) coefficients, element-wise on complex entries).
///
/// Reference: Kolafa, J. Comput. Chem. 25, 335 (2004)
///            Kuhne et al., PRL 98, 066401 (2007)
#pragma once

#include <complex>
#include <deque>
#include <vector>

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/linalg/matrix.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief ASPC density matrix predictor-corrector (real AO density).
///
/// Usage:
/// @code{.cpp}
/// ASPCPredictor pred(3);  // order K=3
/// for (int step = 0; step < n_steps; ++step) {
///     xt::xtensor<double, 2> P_guess;
///     if (pred.predict(P_guess)) {
///         // Use P_guess as SCF initial guess
///         auto result = rhf_solve_with_guess(integrals, P_guess, tol, max_iter);
///         pred.correct(result.P, P_guess, P_corrected);
///         pred.record(P_corrected);
///     } else {
///         // Not enough history, do full SCF
///         auto result = SCFSolver(integrals, tol, max_iter);
///         pred.record(result.P);
///     }
/// }
/// @endcode
class ASPCPredictor {
  public:
    /// @brief Construct with extrapolation order.
    /// @param order Lagrange extrapolation order K (default 3). Uses K+1 history points.
    explicit ASPCPredictor(int order = 3);

    /// @brief Record a converged density matrix into history.
    /// @param P_scf Converged density matrix from SCF.
    void record(const xt::xtensor<double, 2>& P_scf);

    /// @brief Predict density for the next step via Lagrange extrapolation.
    /// @param[out] P_predicted Extrapolated density matrix.
    /// @return True if prediction was made (sufficient history), false otherwise.
    bool predict(xt::xtensor<double, 2>& P_predicted) const;

    /// @brief Blend corrector result with predictor for time-reversibility.
    ///
    /// P_corrected = (K+2)/(2K+3) * P_scf + (K+1)/(2K+3) * P_predicted
    ///
    /// @param P_scf       SCF-converged density at current step.
    /// @param P_predicted Predictor output (from predict()).
    /// @param[out] P_corrected Blended density to record.
    void correct(const xt::xtensor<double, 2>& P_scf,
                 const xt::xtensor<double, 2>& P_predicted,
                 xt::xtensor<double, 2>& P_corrected) const;

    /// @brief Reset history (e.g., after geometry discontinuity or restart).
    void reset();

    /// @brief Number of history points currently stored.
    /// @return Count of stored density matrices.
    int history_size() const { return static_cast<int>(history_.size()); }

    /// @brief Extrapolation order K.
    /// @return Lagrange extrapolation order.
    int order() const { return K_; }

  private:
    int K_; ///< Extrapolation order.
    std::deque<xt::xtensor<double, 2>>
        history_;               ///< Stored density matrices (newest last).
    std::vector<double> omega_; ///< Lagrange extrapolation coefficients.

    /// @brief Compute Lagrange extrapolation coefficients for order K.
    void compute_coefficients();
};

/// @brief ASPC on complex Hermitian densities \(P_k\) (multi-k SGCP).
///
/// Same Kolafa Lagrange \(\omega_k\) and corrector weights as
/// `ASPCPredictor`; coefficients act element-wise on
/// `linalg::CplxMatrix`. Used by `GPWEvaluator` with `second_gen_cp` and a
/// non-Γ k-mesh (one instance per sampled k, dual for UKS).
class ASPCPredictorCplx {
  public:
    /// @brief Construct with extrapolation order.
    /// @param order Lagrange order K (≥ 0); needs K+1 history points.
    /// @throws std::invalid_argument when order < 0.
    explicit ASPCPredictorCplx(int order = 3);

    /// @brief Record a (post-corrector or cold-SCF) density into history.
    /// @param P_scf Complex Hermitian density (n_ao × n_ao).
    void record(const linalg::CplxMatrix& P_scf);

    /// @brief Predict the next-step density via Lagrange extrapolation.
    /// @param[out] P_predicted Extrapolated complex density.
    /// @return true when history depth ≥ K+1.
    bool predict(linalg::CplxMatrix& P_predicted) const;

    /// @brief Kolafa blend \(P_\mathrm{corr} = w_\mathrm{scf} P_\mathrm{scf}
    ///        + w_\mathrm{pred} P_\mathrm{pred}\).
    /// @param P_scf       SCF density at the current step.
    /// @param P_predicted Predictor output from `predict()`.
    /// @param[out] P_corrected Blended density to `record()`.
    void correct(const linalg::CplxMatrix& P_scf,
                 const linalg::CplxMatrix& P_predicted,
                 linalg::CplxMatrix& P_corrected) const;

    /// @brief Drop all history (geometry discontinuity / restart).
    void reset();

    /// @brief Number of stored densities. @return history length.
    int history_size() const { return static_cast<int>(history_.size()); }

    /// @brief Extrapolation order K. @return K ≥ 0.
    int order() const { return K_; }

  private:
    int K_;
    std::deque<linalg::CplxMatrix> history_;
    std::vector<double> omega_;

    void compute_coefficients();
};

} // namespace scf
} // namespace cpu
} // namespace ATV
