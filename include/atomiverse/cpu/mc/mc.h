/// @file cpu/mc/mc.h
/// @brief cpu::mc::MC — concrete canonical (NVT) Metropolis Monte Carlo driver.
#pragma once

#include <cstdint>
#include <optional>
#include <random>
#include <string>
#include <utility>

#include <cstddef>

#include "atomiverse/core/precision.h"
#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/component/zero.h"
#include "atomiverse/cpu/driver/cpu_driver.h"
#include "atomiverse/cpu/driver/outcome_policy.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/mc/transaction.h"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace mc {

/// @brief CPU canonical (NVT) Metropolis Monte Carlo — a **preset**, not a category.
///
/// It runs a budget of trial moves, each proposed by a `Mover` (v1: `TranslateMove`),
/// scored by any registered `Capability::Energy` evaluator (e.g. `PairLJCut`), and kept
/// or discarded by an `OutcomePolicy`. It reuses the MD component set unchanged — energy
/// evaluator, neighbor list, frame — and mirrors the shape of `cpu::md::MD`.
///
/// ## Accept/reject is not a termination
///
/// `SamplingDriver` used to exist as a class whose entire contribution was a hardcoded
/// Metropolis `if` sitting between the stages and the loop's termination test. That `if`
/// is now an `OutcomePolicy` — a second, orthogonal policy slot that every loop has and
/// only this one uses:
///
///   - `LoopHandler`   — should the loop stop?      (MC attaches none; the budget is it)
///   - `OutcomePolicy` — does this iteration count? (MC attaches `Metropolis`, in its
///                                                   constructor, like any other preset)
///
/// Swapping in a GCMC or REMD acceptance rule is now a `set_outcome()` call rather than
/// a new driver category. @see cpu/driver/outcome_policy.h
///
/// The constructor binds `cpu::ZeroEnergy` at the head of the `energy` phase (order -999)
/// so the evaluator accumulates onto a freshly zeroed `frame.scalars.energy` each trial.
/// Categorical setters `set_pair` / `set_neighbor` / `set_observer` / `set_record` /
/// `set_thermo` are inherited from `CPUDriver`; `set_mover` is added here (mirroring MD's
/// `set_integrator`).
///
/// @note For MC the neighbor list should rebuild whenever an accepted move
///       displaces an atom past the skin — construct it with `delay = 0`
///       (`NeighborList<>{cutoff, skin, 1, 0, true}`) so the displacement
///       check is never deferred.
///
/// Usage:
/// @code{.cpp}
///   cpu::mc::MC mc{/*T [K]=*/300.0, /*seed=*/12345u};
///   mc.set_pair(cpu::PairLJCut{eps, sigma, cutoff});
///   mc.set_neighbor(cpu::md::NeighborList<>{cutoff, skin, 1, 0, true});
///   mc.set_mover(cpu::mc::TranslateMove{/*delta [Å]=*/0.2});
///   mc.set_thermo(100);
///   mc.build();
///   mc.run(frame, 100000);            // 100k trial moves
///   double e = mc.current_energy();   // [eV]; mc.acceptance_ratio() in [0,1]
/// @endcode
class MC : public CPUDriver<ATV::sampling::Phase> {
  public:
    /// @brief Construct the NVT sampler.
    /// @param temperature_K Ensemble temperature [K] (must be > 0).
    /// @param seed Seed for the single reproducible RNG stream.
    /// @throws std::invalid_argument if `temperature_K <= 0`.
    MC(AccReal temperature_K, std::uint64_t seed)
        : rng_(seed),
          metropolis_(this->set_outcome(cpu::Metropolis{temperature_K, rng_})) {
        // Zero the energy accumulator at the head of the trial-energy phase so
        // the evaluator accumulates onto a clean slate each move (mirrors MD's
        // ZeroEnergy at pre_force -999).
        this->on(ATV::sampling::Phase::energy).add(cpu::ZeroEnergy{}, -999);
    }

    /// @brief Flush the deferred component `destroy()` of a live continuable chain.
    ~MC() { this->flush_live_(); }

    /// @brief Advance the system by `n_sweeps` Monte Carlo trial moves.
    ///
    /// **Continuable.** Consecutive calls on the same frame extend one Markov chain —
    /// the standard `run(equilibrate); run(sample); …` pattern.
    ///
    /// Per sweep (phases `propose` → `energy` → `end_of_step`):
    ///   1. `propose` — a Mover displaces atoms in place and records the changeset into
    ///      `ctx.mc.proposal`.
    ///   2. `energy` — ZeroEnergy, then the registered evaluator recomputes the
    ///      **full-frame** trial energy. Full-frame (rather than a local sub-frame) is
    ///      the only form valid for *all* potentials, including long-range/kspace where
    ///      a single-atom energy is undefined.
    ///   3. the `OutcomePolicy` (Metropolis) keeps or discards it. On `Rollback` the
    ///      moved atoms are restored from the changeset — the policy decides, the driver
    ///      that owns the changeset undoes.
    ///
    /// **Bootstrap invariant:** the policy is seeded from one full-frame energy
    /// evaluation of the *initial* configuration before the first trial move (step index
    /// `kBootstrapStep` is the energy-only pass; a Mover proposes nothing on it). Without
    /// it the first ΔE is measured against zero and is meaningless.
    ///
    /// @param frame Host-resident system state, updated in place.
    /// @param n_sweeps Number of trial moves to attempt.
    void run(cpu::Frame& frame, int n_sweeps) {
        if (this->begin_run_(frame)) {
            // Publish the per-move changeset + the RNG so a Mover can reach them
            // (mirrors NeighborList publishing pair_list). Addresses are stable.
            this->backend_->context().mc.proposal = &proposal_;
            this->backend_->context().mc.rng = &rng_;
            this->step_once(kBootstrapStep);
            this->start_outcome_(frame); // seed Metropolis's reference energy
        }

        for (int s = 0; s < n_sweeps; ++s) {
            this->mutable_outer_status().reset();
            const int sweep = this->next_step_();
            this->step_once(sweep); // propose (mutate + record) → energy
            if (this->decide_outcome_(frame, sweep) == Outcome::Rollback) {
                rollback_(frame);
                // The frame must reflect the restored state, not the rejected trial.
                frame.scalars.energy = metropolis_.current_energy();
            }
            if (this->poll_handlers_(s) != Verdict::Continue)
                break;
        }
    }

    /// @brief Energy [eV] of the current (last-accepted) configuration.
    AccReal current_energy() const { return metropolis_.current_energy(); }

    /// @brief Total trial moves attempted across all `run()` calls.
    long long attempts() const { return metropolis_.attempts(); }

    /// @brief Total trial moves accepted across all `run()` calls.
    long long accepts() const { return metropolis_.accepts(); }

    /// @brief Accepted fraction of attempted moves, in [0,1] (0 if none attempted).
    double acceptance_ratio() const { return metropolis_.acceptance_ratio(); }

    /// @brief Configure the trial-move generator (mirrors MD `set_integrator`).
    /// Singular: re-call with the same Mover is idempotent, a different Mover
    /// throws. Binds the Mover at its `default_bindings()` (`propose` phase).
    /// @param m Mover component (e.g. `TranslateMove`).
    /// @return `*this` for chaining.
    template <typename M>
    MC& set_mover(M&& m) {
        this->set_singular_default_(std::forward<M>(m), mover_set_, "set_mover");
        return *this;
    }

  protected:
    /// @brief Enforce that a pair evaluator that needs a neighbor list has one.
    void install_defaults() override {
        const auto s = this->scan_pending_infra_();
        this->throw_if_missing_neighbor_list_(
            s, "MC",
            "Add `mc.set_neighbor(cpu::md::NeighborList<>{cutoff, skin, 1, 0, "
            "true})`.");
    }

    /// @brief A Mover to propose with, and something that can score the proposal.
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
    /// @brief Undo a rejected move in place from the recorded changeset.
    ///
    /// The policy decided; this undoes. Only the component that recorded the change
    /// knows how to reverse it, which is why this is not on `OutcomePolicy`.
    /// @param frame System state to restore.
    void rollback_(cpu::Frame& frame) {
        for (std::size_t k = 0; k < proposal_.moved.size(); ++k) {
            const LocalInt i = proposal_.moved[k];
            frame.atoms.x[i] = proposal_.old_x[k];
            frame.atoms.y[i] = proposal_.old_y[k];
            frame.atoms.z[i] = proposal_.old_z[k];
        }
    }

    /// @brief Reserved step index for the energy-only bootstrap pass (a Mover
    ///        proposes nothing on a negative step).
    static constexpr int kBootstrapStep = -1;

    // Declaration order is load-bearing: `rng_` is constructed before `metropolis_`,
    // which binds a reference to it in the mem-init list.

    /// The run's single reproducible RNG stream. Owned here because it is a run-level
    /// service: the Mover draws displacements from it (via `ctx.mc.rng`) and the
    /// acceptance policy draws from it too, interleaved in one sequence per seed.
    std::mt19937_64 rng_;
    /// The acceptance policy, owned by the base's outcome slot. This reference is how
    /// the preset reads back what its own policy decided (energy, acceptance ratio).
    cpu::Metropolis& metropolis_;
    mc::Transaction proposal_;             ///< Reused per-move changeset (on the Context).
    std::optional<std::string> mover_set_; ///< set_mover singular slot.
};

} // namespace mc
} // namespace cpu
} // namespace ATV
