/// @file cpu/driver/path_driver.h
/// @brief cpu::PathDriver — the 6th category: N coupled replicas, all must converge.
///
/// The category CLAUDE.md predicted ("Path for NEB — terminate when every image of a
/// coupled chain has converged"). The prediction was half right: the termination rule is
/// indeed new, but "N coupled replicas" turned out to be a **second axis** of the
/// taxonomy rather than a sixth point on the first one — see `cpu/replica/replica_set.h`,
/// which is where the replica machinery lives so that a future budget-terminating
/// replica category (REMD) reuses it instead of cloning it.
///
/// ## Termination
///
/// `ReplicaSet::all_converged()` — an AND over every replica. One unconverged image
/// keeps the whole chain running.
///
/// ## The loop
///
/// @code
///   expand one Frame -> N replicas (ReplicaExpansion)
///   repeat:
///     for each replica i:            <- force is PER-replica
///        point the backend at replica i's Frame; run the `force` phase
///     run `couple`                   <- band-wide: sees all replicas via ctx.path
///     run `search_direction` .. `end_of_step`
///     stop when every replica reports converged
///   collapse N replicas -> the caller's one Frame
/// @endcode
///
/// ## Why ONE backend and not N
///
/// The spec for this probe called for a Backend (and therefore a Context, and therefore
/// an `outer_status`) **per replica**, which would have given collective termination for
/// free. It is not implementable, and finding that out is one of this probe's results:
/// **components cannot be cloned.** `Driver::add()` *moves* a component into
/// `owned_components_`, and there is no `clone()` anywhere in the tree — so N backends
/// would need N independently-constructed component sets, which the user would have to
/// supply through a factory callback (against both the OOP style rule and any hope of a
/// decent API), or which would need a `clone()` virtual bolted onto `ComponentBase` —
/// a dead interface for every component that is not in a band.
///
/// So: one backend, one plan, and the band reached through `ctx.path`. The `force` phase
/// is replayed once per replica by re-pointing the backend's frame; the band-wide phases
/// run once and see everything. Per-replica convergence lives on `ReplicaSet`, not on N
/// copies of `IterationStatus`, and `PathDriver` ANDs it into the one `outer_status`.
///
/// Known ceiling of the one-backend route: force-phase components that cache per-frame
/// state (a neighbor list's bins) get invalidated on every replica switch and rebuild N
/// times per iteration. Correct, but O(N) rebuilds. A band of Müller-Brown particles has
/// no neighbor list at all; a band of a real solid would want a per-replica neighbor
/// cache, which is the point at which a component *would* need to be replica-aware.
#pragma once

#include <stdexcept>

#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/driver/cpu_driver.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/replica/replica_set.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {

class PathDriver : public CPUDriver<ATV::path::Phase> {
  public:
    /// @brief Advance the band until every replica converges, or `max_iter` runs out.
    ///
    /// @param frame Initial state. Expanded into the band on the way in; **overwritten**
    ///        with the collapsed result on the way out (for CI-NEB: the saddle-point
    ///        geometry). Frame in, Frame out — the facade currency is unchanged.
    /// @param expansion The strategy that generates the band. Borrowed for the call.
    /// @param max_iter Iteration budget for the whole chain.
    /// @return Iterations actually executed.
    int run(cpu::Frame& frame, const ReplicaExpansion& expansion, int max_iter) {
        ReplicaSet band(frame, expansion);

        // setup_pass binds the backend to replica 0's frame. Every replica has the same
        // atom count and topology (they are expansions of ONE frame), so components size
        // their buffers once and stay valid as the active frame is swapped underneath.
        setup_pass(band.frame(0));

        auto& ctx = this->backend_->context();
        ctx.path.replicas = &band;

        int iters = 0;
        for (int s = 0; s < max_iter; ++s) {
            mutable_outer_status().reset();

            // ── Per-replica: force ────────────────────────────────────────────
            for (int i = 0; i < band.size(); ++i) {
                ctx.path.active = i;
                this->backend_->set_frame(band.frame(i));
                this->backend_->run_stage(ATV::path::Phase::force, s);
            }
            ctx.path.active = -1;

            // ── Band-wide: couple → direction → line search → convergence → eos ──
            // These see every replica through ctx.path.replicas. The active frame is
            // left pointing at the last replica and is deliberately not meaningful
            // here; a band-wide component that reads ctx.frame is reading one arbitrary
            // image, which is a bug in that component.
            for (const auto p : {ATV::path::Phase::couple,
                                 ATV::path::Phase::search_direction,
                                 ATV::path::Phase::line_search,
                                 ATV::path::Phase::convergence_check,
                                 ATV::path::Phase::end_of_step}) {
                this->backend_->run_stage(p, s);
            }

            ++iters;
            if (band.all_converged()) {
                mutable_outer_status().converged = true;
                break;
            }
            if (outer_status().stop_requested)
                break;
        }

        band.collapse_into(frame, expansion);
        ctx.path.replicas = nullptr;
        teardown_pass();
        return iters;
    }

  protected:
    /// A Path run needs something to compute forces on each replica and something to
    /// move them. The coupling component is NOT required here: an uncoupled band is
    /// just N independent minimizations, which is a legitimate (if pointless) thing to
    /// run, and refusing it would mean this driver knew what NEB is.
    void validate_components() const override {
        bool has_optimizer = false;
        bool has_force_evaluator = false;
        for_each_pending([&](const cpu::ComponentBase& c) {
            if (c.role() == ComponentRole::Optimizer)
                has_optimizer = true;
            if (c.role() == ComponentRole::Evaluator
                && has_capability(c.capabilities(), Capability::Force))
                has_force_evaluator = true;
        });
        if (!has_force_evaluator)
            throw_missing("Evaluator with Capability::Force");
        if (!has_optimizer)
            throw_missing("Optimizer");
    }
};

} // namespace cpu
} // namespace ATV
