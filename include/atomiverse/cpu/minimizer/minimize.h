/// @file cpu/minimizer/minimize.h
/// @brief cpu::minimizer::Minimize — concrete energy-minimizer driver.
///
/// `Minimize` is the user-facing CPU outer-driver under the
/// `OptimizeDriver` category. It iterates an optimizer component
/// (LBFGS / CG / future FIRE) until the optimizer reports convergence
/// through `outer_status`, or until `max_iter` is reached.
///
/// Exposes the same LAMMPS-shaped categorical setter API as
/// `cpu::md::MD` (minus the time-evolution setters): set_pair,
/// set_bond, set_angle, set_dihedral, set_neighbor, set_optimizer,
/// set_observer, set_record, set_thermo — all except set_optimizer
/// are inherited from CPUDriver.
///
/// `add()` and `on(phase).add()` from the underlying `Driver<>`
/// remain public for user-defined hooks.
///
/// **MPI / domain decomposition.** Minimize mirrors `cpu::md::MD`'s
/// infrastructure auto-injection, folded into the single `force` phase
/// (the minimizer skeleton has no dedicated exchange/neighbor phases):
///
/// @code
///   force @ -1000  ZeroForces        (ctor)
///   force @  -999  ZeroEnergy        (ctor)
///   force @   -50  comm::Migrate     (install_defaults; re-homes owners,
///                                     single-rank no-op)
///   force @   -30  comm::HaloSend    (install_defaults, ATV_HAS_MPI +
///                                     pair evaluator; deduped)
///   force @   -10  NeighborList      (set_neighbor)
///   force @     0  pair / bonded evaluators
///   force @   +10  comm::HaloRecv    (install_defaults, gated as above)
///   force @   +20  md::mpi::Reduce   (install_defaults; makes
///                                     frame.scalars global *before* the
///                                     convergence_check phase reads them)
/// @endcode
///
/// The optimizer components are themselves DD-aware (owned-row loops,
/// globally-reduced inner products and convergence metrics — see
/// `lbfgs.h`). The recorders are DD-correct: `TrajectoryWriter` /
/// `FinalGeometry` gather the owner rows to rank 0
/// (`md::mpi::gather_to_root`, tag-ordered) before writing, and
/// `OptLog` reduces its metrics across ranks; all three write on
/// rank 0 only.
#pragma once

#include <string>
#include <type_traits>

#include "atomiverse/cpu/async/snapshot_publisher.h"
#include "atomiverse/cpu/comm/atom_migration.h"
#include "atomiverse/cpu/comm/halo_exchange.h"
#include "atomiverse/cpu/comm/reduce.h"
#include "atomiverse/cpu/component/pair_evaluator.h"
#include "atomiverse/cpu/component/zero.h"
#include "atomiverse/cpu/driver/cpu_driver.h"
#include "atomiverse/cpu/frame/materialize.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace minimizer {

/// @brief Energy-minimization driver — convergence-aware Optimize loop.
///
/// Auto-adds ZeroForces and ZeroEnergy at force phase. Wires the shared
/// `cpu::Publish` into `end_of_step @ +10` so async observers
/// see snapshots. `install_defaults()` injects the MPI infrastructure
/// documented in the file header (single-rank: structural no-ops).
///
/// Validation rules inherited from `OptimizeDriver`: requires ≥ 1
/// `Optimizer` component (e.g. `LBFGS`, `CG`) and ≥ 1 `Evaluator`
/// with `Capability::Energy` or `Capability::Force` (e.g. `SCFEvaluator`).
///
/// Inheritance: Minimize → CPUDriver<minimizer::Phase> →
/// Driver<minimizer::Phase, cpu::device_tag>. There is no `OptimizeDriver` category in
/// between any more: "stop when converged" is a `ForceConverged` handler this preset
/// attaches in `install_defaults()`, and the user can replace it without touching a
/// driver class or an optimizer.
/// CPUDriver provides the shared categorical setters and slot variables.
class Minimize : public CPUDriver<ATV::minimizer::Phase> {
  public:
    /// @brief Iterate to convergence, or up to `max_iter` iterations.
    ///
    /// Unlike MD (continuable: consecutive calls extend one trajectory), each `run()`
    /// here is a complete, self-contained optimization — `setup_pass` re-runs component
    /// `build()`, wiping the optimizer's curvature history, and `teardown_pass` runs
    /// before returning. That is deliberate: a converged geometry has no useful
    /// "continuation state", so re-running means re-optimizing (possibly after the caller
    /// perturbed the system) and stale L-BFGS history would be wrong to keep.
    ///
    /// @param frame Host-resident system state, updated in place.
    /// @param max_iter Iteration budget.
    /// @return Number of iterations actually executed (≤ `max_iter`).
    int run(cpu::Frame& frame, int max_iter) {
        setup_pass(frame);
        // Publish the re-entrant force re-eval handle so an Optimizer component can
        // evaluate trial geometries during the line_search phase. Replays the force
        // phase (zeroing → evaluator → reduce) on the frame's current geometry.
        // Captures `this`; valid for the run.
        this->backend_->context().opt.reeval_forces = [this](bool /*need_grad*/) {
            this->backend_->run_stage(Traits::eval_phase, /*step_idx=*/0);
        };
        int iters = 0;
        for (int s = 0; s < max_iter; ++s) {
            mutable_outer_status().reset();
            step_once(s);
            ++iters;
            // The optimizer PUBLISHED its metrics during the step; the handlers judge
            // them here. Neither LBFGS nor CG decides its own convergence any more —
            // which is why a nested optimizer can no longer stop this loop.
            if (poll_handlers_(s) != Verdict::Continue)
                break;
        }
        teardown_pass();
        return iters;
    }

    /// @brief User-facing entry — optimize the system stored in a molrs.Frame.
    /// @param frame molrs.Frame handle (input geometry + output sink).
    /// @param max_iter Iteration budget.
    /// @return Number of iterations actually executed (≤ `max_iter`).
    int run(molrs::FrameRef& frame, int max_iter) {
        cpu::Frame work = cpu::materialize(frame);
        const int iters = run(work, max_iter);
        cpu::write_back(work, frame);
        return iters;
    }

    Minimize() {
        // Phase-agnostic backend-root Zero*, bound to minimizer phases.
        this->on(ATV::minimizer::Phase::force)
            .add(cpu::ZeroForces{}, -1000)
            .add(cpu::ZeroEnergy{}, -999);
        this->on(ATV::minimizer::Phase::end_of_step).add(cpu::Publish{}, 10);
    }

    // ─── Minimizer-specific categorical setters ─────────────────────
    //
    // set_pair, set_bond, set_angle, set_dihedral, set_neighbor,
    // set_observer, set_record, set_thermo are inherited from CPUDriver.

    /// @brief Configure the optimization algorithm (LAMMPS `min_style`).
    template <typename O>
    Minimize& set_optimizer(O&& o) {
        this->set_singular_default_(std::forward<O>(o), optimizer_set_,
                                    "set_optimizer");
        return *this;
    }

    /// @brief Owner-migration cadence in optimizer iterations (default 10,
    /// matching `HaloSend`'s full-rebuild cadence). Each actual migration
    /// discards the optimizer's curvature history (row indices reshuffle),
    /// so raise this when minimizing systems that barely drift; the
    /// exchange itself is skipped whenever no rank has a strayed atom.
    /// Call before `build()`.
    Minimize& set_migrate_every(int every) {
        migrate_every_ = every < 1 ? 1 : every;
        return *this;
    }

  protected:
    /// @brief An Optimizer, and an Evaluator that provides Energy or Force.
    void validate_components() const override {
        bool has_optimizer = false;
        bool has_energy_or_force_evaluator = false;
        for_each_pending([&](const cpu::ComponentBase& c) {
            if (c.role() == ComponentRole::Optimizer)
                has_optimizer = true;
            if (c.role() == ComponentRole::Evaluator
                && (has_capability(c.capabilities(), Capability::Energy)
                    || has_capability(c.capabilities(), Capability::Force)))
                has_energy_or_force_evaluator = true;
        });
        if (!has_optimizer)
            throw_missing("Optimizer");
        if (!has_energy_or_force_evaluator)
            throw_missing(
                "Evaluator with Capability::Energy or Capability::Force");
    }

    void install_defaults() override {
        // A Minimize with no handler would never stop early — the optimizer publishes
        // metrics and nobody judges them. So the PRESET supplies the policy it is a
        // preset for. This is what "MD-style pipelines become a default-constructed
        // driver + handler" means in practice: the class contributes a default, not a
        // category.
        //
        // Textbook LAMMPS/ASE-ish defaults; override wholesale with
        // `opt.add_handler(cpu::ForceConverged{ftol, ...})` before build().
        if (!this->has_handler())
            this->add_handler(cpu::ForceConverged{});

        const auto s = this->scan_pending_infra_();

        // Owner migration + global scalar reduction, explicitly bound to
        // minimizer phases (their default_bindings() are md::Phase values
        // and must not be used here). Both are single-rank no-ops. Reduce
        // sits after HaloRecv (+10) and before convergence_check so the
        // optimizer reads globally-summed energy.
        this->on(ATV::minimizer::Phase::force)
            .add(comm::Migrate{migrate_every_}, -50);
        this->on(ATV::minimizer::Phase::force).add(md::mpi::Reduce{}, +20);

#ifdef ATV_HAS_MPI
        if (s.neighbor_consumer_cutoff) {
            // Same 10% skin rule as MD::install_defaults — ghosts must stay
            // registered as atoms drift between the every-10-step full halo
            // rebuilds. The minimizer's step cap (0.1 Å) keeps per-step
            // drift well inside this band. Explicit halo registrations
            // suppress the auto-inject: the reverse exchange must run
            // exactly once (see halo_exchange.h).
            const double halo_cut =
                1.1 * static_cast<double>(*s.neighbor_consumer_cutoff);
            if (!s.has_halo_send)
                this->on(ATV::minimizer::Phase::force)
                    .add(comm::HaloSend{halo_cut}, -30);
            if (!s.has_halo_recv)
                this->on(ATV::minimizer::Phase::force)
                    .add(comm::HaloRecv{}, +10);
        }
#endif

        this->throw_if_missing_neighbor_list_(
            s, "Minimize",
            "Add `opt.set_neighbor(cpu::md::NeighborList<>{cutoff, skin})`.");
    }

  private:
    std::optional<std::string> optimizer_set_;
    int migrate_every_ = 10;
};

} // namespace minimizer
} // namespace cpu
} // namespace ATV
