/// @file cpu/scf/stages/mix_stage.h
/// @brief cpu::scf::MixStage — the `scf::Phase::mix` stage component: DIIS
///        error + Fock-space extrapolation + level shift, bit-identical to the
///        shipping SCF loop, with the DIIS history owned inside the component.
///
/// Wraps the shipping internal helpers `internal::diis_error`,
/// `internal::DIISMixer` and `internal::apply_level_shift` behind a
/// `cpu::Component` bound at `scf::Phase::mix`. State flows through named
/// constructor references (const& reads, & writes) — no config struct, no
/// workspace bag, no Context slice. `launch()` reproduces the shipping loop
/// body (`scf_solver.cpp` steps 2 + 2b) byte-for-byte: build the DIIS error
/// from the UNSHIFTED Fock, extrapolate F in Fock space (dual channel when
/// unrestricted, single channel when restricted), then apply the
/// Saunders–Hillier level shift.
///
/// The DIIS subspace history is a private accelerator state of this component
/// (the `DIISMixer` deque); `reset()` clears it at the start of each run. Its
/// per-iteration deque growth in `launch()` is the documented allocation
/// exception — no other fixed-buffer allocation happens in `launch()`.
///
/// DIIS: Pulay, Chem. Phys. Lett. 73, 393 (1980); J. Comp. Chem. 3, 556 (1982).
/// Level shift: Saunders & Hillier, Int. J. Quantum Chem. 7, 699 (1973).
#pragma once

#include <string>
#include <vector>

#include "atomiverse/core/binding.h"
#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/scf/dense_matrix.h"
#include "atomiverse/cpu/scf/internal/diis_mixer.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief SCF `mix` stage: DIIS error, Fock extrapolation and level shift.
class MixStage : public cpu::Component {
  public:
    /// @brief Wire the mix stage to the method driver's iteration state.
    /// @param restricted       True drives a single (alpha) spin channel; false
    ///                         drives both alpha and beta.
    /// @param use_diis         Enable DIIS extrapolation.
    /// @param diis_start_iter  First iteration index at which DIIS activates.
    /// @param diis_max_vectors Maximum stored (F, e) pairs in the DIIS space.
    /// @param level_shift      Virtual-orbital level shift lambda [Hartree];
    ///                         |lambda| <= 1e-4 is a strict no-op.
    /// @param S                Overlap matrix (dimensionless), n x n [read].
    /// @param X                Canonical orthogonalization transform, n x n
    ///                         [read].
    /// @param P_alpha          Alpha density matrix [read].
    /// @param P_beta           Beta density matrix [read].
    /// @param iter             Current SCF iteration index [read].
    /// @param F_alpha          Alpha Fock matrix [Hartree], in/out [write].
    /// @param F_beta           Beta Fock matrix [Hartree], in/out [write].
    /// @param err_alpha        Alpha DIIS error matrix, output [write].
    /// @param err_beta         Beta DIIS error matrix, output [write].
    MixStage(bool restricted, bool use_diis, int diis_start_iter,
             int diis_max_vectors, double level_shift, const RealMatrix& S,
             const RealMatrix& X, const RealMatrix& P_alpha,
             const RealMatrix& P_beta, const int& iter, RealMatrix& F_alpha,
             RealMatrix& F_beta, RealMatrix& err_alpha, RealMatrix& err_beta)
        : restricted_(restricted), use_diis_(use_diis),
          diis_start_iter_(diis_start_iter),
          diis_max_vectors_(diis_max_vectors), level_shift_(level_shift), S_(S),
          X_(X), P_alpha_(P_alpha), P_beta_(P_beta), iter_(iter),
          F_alpha_(F_alpha), F_beta_(F_beta), err_alpha_(err_alpha),
          err_beta_(err_beta), diis_(diis_max_vectors) {}

    /// @brief Component identity.
    /// @return The literal "scf_mix_stage".
    std::string name() const override { return "scf_mix_stage"; }

    /// @brief Default binding: `scf::Phase::mix` @ 0.
    /// @return Single binding placing this at the SCF mix phase.
    std::vector<DefaultBinding> default_bindings() const override {
        return {ATV::make_binding(ATV::scf::Phase::mix, 0)};
    }

    /// @brief DIIS error + Fock-space extrapolation + level shift, in place.
    /// @param ctx Launch context; work runs only at `scf::Phase::mix`.
    void launch(const cpu::LaunchContext& ctx) override;

    /// @brief Clear the DIIS subspace history (called once per run start).
    void reset() { diis_ = internal::DIISMixer(diis_max_vectors_); }

  private:
    bool restricted_;
    bool use_diis_;
    int diis_start_iter_;
    int diis_max_vectors_;
    double level_shift_;

    const RealMatrix& S_;
    const RealMatrix& X_;
    const RealMatrix& P_alpha_;
    const RealMatrix& P_beta_;
    const int& iter_;

    RealMatrix& F_alpha_;
    RealMatrix& F_beta_;
    RealMatrix& err_alpha_;
    RealMatrix& err_beta_;

    internal::DIISMixer diis_; ///< Private DIIS subspace history.
};

} // namespace scf
} // namespace cpu
} // namespace ATV
