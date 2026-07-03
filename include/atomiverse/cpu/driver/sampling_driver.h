/// @file cpu/driver/sampling_driver.h
/// @brief cpu::SamplingDriver — Monte Carlo (accept/reject) driver category.
///
/// `SamplingDriver` is the fifth CPU outer-driver category (siblings:
/// Evaluate / Dynamics / Optimize / Response), reserved by CLAUDE.md's
/// driver taxonomy for "Sampling (MC/HMC/REMD)" on a genuinely distinct
/// outer-loop semantic: a fixed sweep budget of stochastic trial moves, each
/// kept or discarded by a per-move acceptance test. It runs the requested
/// number of sweeps; like `DynamicsDriver` it does **not** interpret
/// `outer_status.converged` (sampling has no convergence success path) but
/// honours `outer_status.stop_requested` for hard aborts.
///
/// Per-sweep semantics (phases `propose` → `energy` → `end_of_step`):
///   1. `propose` — a `Mover` component displaces atoms in place and records
///      the changeset into `Context::proposal`.
///   2. `energy`  — `ZeroEnergy` then the registered `Capability::Energy`
///      evaluator recompute the **full-frame** trial energy `frame.scalars.energy`.
///      Full-frame (rather than a local sub-frame) is the only form valid for
///      *all* potentials, including long-range/kspace where a single-atom
///      energy is undefined.
///   3. accept/reject — the driver applies the inline **Metropolis** test
///      `min(1, exp(-ΔE/(k_B T)))` against the cached accepted energy `e_cur_`
///      (`ΔE = E_trial - e_cur_`) and either commits (`e_cur_ = E_trial`, frame
///      already holds the move) or rolls back the moved atoms from
///      `Context::proposal->old_*`.
///
/// **Bootstrap invariant:** `e_cur_` is seeded from one full-frame energy
/// evaluation of the *initial* configuration before the first trial move
/// (step index `< 0` is the energy-only bootstrap pass — a Mover proposes
/// nothing on it). Without this the first ΔE is meaningless.
///
/// Validation (build-time): ≥ 1 component with `role()==ComponentRole::Mover`
/// AND ≥ 1 component advertising `Capability::Energy`.
///
/// Concrete user-facing class: `cpu::mc::MC` (added in mc-sampling-02).
///
/// *Alternative considered — rejected for v1:* folding accept/reject inside
/// the Mover component (driver a dumb loop, nothing published on Context).
/// Rejected in favour of the inline-in-driver form here; revisit when a
/// swappable acceptance policy (REMD/GCMC) makes the seam worthwhile.
#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <stdexcept>

#include "atomiverse/core/precision.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/driver/cpu_driver.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/mc/transaction.h"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {

/// @brief CPU Monte Carlo (accept/reject) driver category base.
class SamplingDriver : public CPUDriver<ATV::sampling::Phase> {
  public:
    /// @brief Advance the system by `n_sweeps` Monte Carlo trial moves.
    ///
    /// Seeds `e_cur_` from the initial configuration (bootstrap invariant),
    /// then per sweep proposes a move, recomputes the full-frame trial energy,
    /// and applies the inline Metropolis test — committing or rolling back.
    /// Honours `outer_status.stop_requested` (early break); ignores
    /// `outer_status.converged` (Sampling has no convergence success path).
    /// @param frame Host-resident system state, updated in place.
    /// @param n_sweeps Number of trial moves to attempt.
    void run(cpu::Frame& frame, int n_sweeps) {
        setup_pass(frame);

        // Publish the per-move changeset + RNG so a Mover can reach them
        // (mirrors NeighborList publishing pair_list). Addresses are stable
        // for the whole run.
        this->backend_->context().proposal = &proposal_;
        this->backend_->context().rng = &rng_;

        // Bootstrap invariant: seed e_cur_ from the initial configuration's
        // energy BEFORE the first trial move. Step index < 0 is the energy-only
        // bootstrap pass; a Mover proposes nothing on it (see TranslateMove).
        step_once(kBootstrapStep);
        e_cur_ = frame.scalars.energy;

        for (int s = 0; s < n_sweeps; ++s) {
            mutable_outer_status().reset();
            step_once(s); // propose (mutate + record) → energy (full recompute)
            const AccReal e_trial = frame.scalars.energy;
            const AccReal dE = e_trial - e_cur_;
            ++n_attempts_;
            if (metropolis_accept_(dE)) {
                e_cur_ = e_trial; // commit — frame already holds the move
                ++n_accepts_;
            } else {
                rollback_(frame);              // restore moved atoms
                frame.scalars.energy = e_cur_; // frame reflects the restored state
            }
            if (outer_status().stop_requested)
                break;
        }

        teardown_pass();
    }

    /// @brief Energy [eV] of the current (last-accepted) configuration —
    ///        equals a fresh full-frame recompute of the current state.
    AccReal current_energy() const { return e_cur_; }

    /// @brief Total trial moves attempted across all `run()` calls.
    long long attempts() const { return n_attempts_; }

    /// @brief Total trial moves accepted across all `run()` calls.
    long long accepts() const { return n_accepts_; }

    /// @brief Accepted fraction of attempted moves (0 if none attempted yet).
    double acceptance_ratio() const {
        return n_attempts_ > 0
                   ? static_cast<double>(n_accepts_) / static_cast<double>(n_attempts_)
                   : 0.0;
    }

  protected:
    /// @brief Construct the sampler with its ensemble temperature and RNG seed.
    /// @param temperature_K Ensemble temperature [K] (must be > 0).
    /// @param seed Seed for the single reproducible RNG stream.
    /// @throws std::invalid_argument if `temperature_K <= 0`.
    explicit SamplingDriver(AccReal temperature_K, std::uint64_t seed)
        : temperature_(temperature_K), rng_(seed) {
        if (!(temperature_K > AccReal(0))) {
            throw std::invalid_argument(
                "SamplingDriver: temperature must be > 0 K");
        }
    }

    /// @brief Require ≥ 1 Mover and ≥ 1 Capability::Energy evaluator.
    /// @throws std::runtime_error naming the first missing requirement.
    void validate_components() const override {
        bool has_mover = false;
        bool has_energy = false;
        for_each_pending([&](const cpu::ComponentBase& c) {
            if (c.role() == ComponentRole::Mover)
                has_mover = true;
            if (has_capability(c.capabilities(), Capability::Energy))
                has_energy = true;
        });
        if (!has_mover)
            throw_missing("Mover");
        if (!has_energy)
            throw_missing("Component with Capability::Energy");
    }

  private:
    /// @brief Metropolis accept test for a trial energy change.
    /// @param dE Trial energy change `E_trial - e_cur_` [eV].
    /// @return true iff the move is accepted (always true for `dE <= 0`).
    bool metropolis_accept_(AccReal dE) {
        // ponytail: canonical Metropolis, hardcoded for v1. A swappable
        // acceptance-policy seam (REMD replica swap, GCMC μVT bias) plugs in
        // here when those ensembles land — one branch until then.
        if (dE <= AccReal(0))
            return true;
        const double beta =
            1.0 / (ATV::units::kB_ev_per_k * static_cast<double>(temperature_));
        const double p = std::exp(-static_cast<double>(dE) * beta);
        std::uniform_real_distribution<double> u(0.0, 1.0);
        return u(rng_) < p;
    }

    /// @brief Undo a rejected move in place from the recorded changeset.
    /// @param frame System state to restore.
    void rollback_(cpu::Frame& frame) {
        for (std::size_t k = 0; k < proposal_.moved.size(); ++k) {
            const LocalInt i = proposal_.moved[k];
            frame.atoms.x[i] = proposal_.old_x[k];
            frame.atoms.y[i] = proposal_.old_y[k];
            frame.atoms.z[i] = proposal_.old_z[k];
        }
    }

    /// @brief Reserved step index for the energy-only bootstrap pass.
    static constexpr int kBootstrapStep = -1;

    AccReal temperature_;      ///< Ensemble temperature [K].
    std::mt19937_64 rng_;      ///< Single reproducible RNG stream.
    mc::Transaction proposal_; ///< Reused per-move changeset (published on Context).
    AccReal e_cur_ = 0;        ///< Cached energy of the last accepted state [eV].
    long long n_attempts_ = 0; ///< Trial moves attempted (lifetime).
    long long n_accepts_ = 0;  ///< Trial moves accepted (lifetime).
};

} // namespace cpu
} // namespace ATV
