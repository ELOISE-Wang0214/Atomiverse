/// @file cpu/gpw/xl_bomd.h
/// @brief XlBomd / XlBomdCplx — dissipative extended-Lagrangian BOMD
///        auxiliary-density propagators (Niklasson XL-BOMD).
///
/// Propagates an auxiliary AO density matrix P alongside the nuclei and
/// serves it as the SCF initial guess each MD step:
///
///   P_{n+1} = 2 P_n − P_{n−1} + κ (D_n − P_n) + α Σ_{k=0}^{K} c_k P_{n−k}
///
/// where D_n is the SCF-converged density at step n and (κ, α, c_k) are the
/// dissipation coefficients of Table I in Niklasson et al., J. Chem. Phys.
/// 130, 214109 (2009) (time-reversible scheme: Niklasson, Phys. Rev. Lett.
/// 100, 123004 (2008)). Unlike the Kolafa ASPC *extrapolation*
/// (`cpu::scf::ASPCPredictor`), the auxiliary density is a dynamical
/// variable: the recursion is time-reversible up to the O(α) dissipation,
/// which is what removes the systematic energy drift that plain guess reuse
/// causes in long BOMD trajectories.
///
/// Unit-agnostic dense linear algebra on AO density matrices; consumed by
/// `cpu::gpw::GPWEvaluator` (`xl_bomd = true`):
/// - **Γ RKS:** total real \(P\) via `XlBomd`
/// - **Γ UKS:** one `XlBomd` per spin channel
/// - **multi-k RKS:** one `XlBomdCplx` per sampled \(k\) (complex \(P_k\))
#pragma once

#include <complex>
#include <deque>
#include <vector>

#include "atomiverse/cpu/linalg/matrix.h"

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief Dissipative XL-BOMD auxiliary-density propagator (real AO density).
///
/// Lifecycle: construct with the dissipation order K, then once per MD step
/// (after the SCF solve) call `advance(D_scf)`. The first `advance` seeds
/// the history (stationary start: P ≡ D); every later call propagates. The
/// SCF guess for the *next* step is `aux()`.
class XlBomd {
  public:
    /// @brief Construct with the Niklasson dissipation order.
    /// @param k Dissipation order K (history depth K+1); Table I of
    ///          JCP 130, 214109 (2009) covers K ∈ [3, 9].
    /// @throws std::invalid_argument when k is outside [3, 9].
    explicit XlBomd(int k = 5);

    /// @brief Dissipation order K. @return K in [3, 9].
    int order() const { return k_; }
    /// @brief Verlet coupling κ (dimensionless). @return κ of Table I.
    double kappa() const { return kappa_; }
    /// @brief Dissipation strength α (dimensionless). @return α of Table I.
    double alpha() const { return alpha_; }
    /// @brief Dissipation weights c_0..c_K (dimensionless).
    /// @return K+1 coefficients of Table I (integer-valued).
    const std::vector<double>& c() const { return c_; }

    /// @brief Whether an auxiliary density is available (history seeded).
    /// @return true after the first `advance`.
    bool ready() const { return ready_; }

    /// @brief Current auxiliary density — the SCF initial guess for the
    ///        next MD step (same scale/convention as the SCF density fed to
    ///        `advance`: RKS total \(P\), or one UKS spin channel).
    /// @return P_{n+1} of the recursion above.
    /// @throws std::logic_error before the first `advance`.
    const linalg::RealMatrix& aux() const;

    /// @brief Record this step's SCF density and propagate the auxiliary
    ///        density one step.
    ///
    /// First call: seeds the whole history with `d_scf` (stationary start —
    /// the recursion's fixed point for a constant density). Later calls:
    /// applies the dissipative Verlet recursion.
    /// @param d_scf SCF-converged density matrix D_n (n_ao × n_ao).
    /// @throws std::invalid_argument when `d_scf` is empty or its shape
    ///         differs from the seeded history.
    void advance(const linalg::RealMatrix& d_scf);

    /// @brief Drop all history (geometry discontinuity / restart policy).
    void reset();

  private:
    int k_;                 ///< dissipation order K.
    double kappa_ = 0.0;    ///< Verlet coupling κ.
    double alpha_ = 0.0;    ///< dissipation strength α.
    std::vector<double> c_; ///< dissipation weights c_0..c_K.
    /// Auxiliary-density history, newest first: P_n, P_{n−1}, …, P_{n−K}.
    std::deque<linalg::RealMatrix> hist_;
    linalg::RealMatrix aux_; ///< P_{n+1} (next-step SCF guess).
    bool ready_ = false;     ///< history seeded.
};

/// @brief Dissipative XL-BOMD propagator on complex Bloch densities \(P_k\).
///
/// Same Table-I recursion as `XlBomd`; real \((\kappa,\alpha,c_k)\) act
/// element-wise on `linalg::CplxMatrix`. A Hermitian history remains
/// Hermitian (up to floating-point roundoff). Used for multi-k RKS XL-BOMD
/// (one instance per sampled k).
class XlBomdCplx {
  public:
    /// @brief Construct with the Niklasson dissipation order.
    /// @param k Dissipation order K ∈ [3, 9] (Table I).
    /// @throws std::invalid_argument when k is outside [3, 9].
    explicit XlBomdCplx(int k = 5);

    int order() const { return k_; }
    double kappa() const { return kappa_; }
    double alpha() const { return alpha_; }
    const std::vector<double>& c() const { return c_; }
    bool ready() const { return ready_; }

    /// @brief Current auxiliary complex density \(P_{n+1}\).
    /// @throws std::logic_error before the first `advance`.
    const linalg::CplxMatrix& aux() const;

    /// @brief Record SCF \(D_k\) and propagate one step.
    /// @param d_scf Complex Hermitian density (n_ao × n_ao).
    /// @throws std::invalid_argument on empty or shape-mismatched input.
    void advance(const linalg::CplxMatrix& d_scf);

    /// @brief Drop all history.
    void reset();

  private:
    int k_;
    double kappa_ = 0.0;
    double alpha_ = 0.0;
    std::vector<double> c_;
    std::deque<linalg::CplxMatrix> hist_;
    linalg::CplxMatrix aux_;
    bool ready_ = false;
};

} // namespace gpw
} // namespace cpu
} // namespace ATV
