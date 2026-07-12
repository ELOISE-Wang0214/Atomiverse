/// @file cpu/driver/cpu_driver.h
/// @brief cpu::CPUDriver — the CPU-backend layer of the driver hierarchy.
///
/// Hierarchy (two orthogonal axes — backend & category):
///
/// @code
///   Driver<PhaseEnum, DeviceTag>          ← plan / dispatch primitive
///     ├─ cpu::CPUDriver<PhaseEnum>         ← CPU backend layer  (this file)
///     │    ├─ cpu::EvaluateDriver
///     │    ├─ cpu::DynamicsDriver
///     │    ├─ cpu::OptimizeDriver
///     │    ├─ cpu::ResponseDriver
///     │    └─ cpu::SamplingDriver
///     └─ cuda::CUDADriver<PhaseEnum>       ← CUDA backend layer
///          ├─ cuda::EvaluateDriver           (cuda_driver.cuh)
///          ├─ cuda::DynamicsDriver
///          └─ cuda::OptimizeDriver
/// @endcode
///
/// `CPUDriver<P>` inherits the underlying `Driver<P, cpu::device_tag>`
/// (which holds the Plan and Backend ownership) and adds the
/// CPU-specific machinery the five category bases share:
///
///   1. A `build()` override that runs `validate_components()` before
///      the Plan is compiled. Subclasses override
///      `validate_components()` to enforce their required role +
///      capability set.
///   2. Protected helpers `setup_pass`, `step_once`, `teardown_pass`
///      that subclass `run(...)` implementations use to drive the
///      per-iteration loop with category-correct termination.
///   3. Read/write access to the outer-driver `IterationStatus`
///      (`cpu::Context::outer_status`), which `OptimizeDriver` polls
///      for `converged` and `DynamicsDriver` polls for
///      `stop_requested`.
///   4. Unified categorical setter API (set_pair, set_bond, set_angle,
///      set_dihedral, set_neighbor, set_observer, set_record, set_thermo)
///      shared across MD / Minimize / Evaluate via `CategoryTraits`.
///
/// User code does not subclass this directly — the five concrete CPU
/// category bases (`cpu::EvaluateDriver` / `DynamicsDriver` /
/// `OptimizeDriver` / `ResponseDriver` / `SamplingDriver`) do.
#pragma once

#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

#include "atomiverse/core/device.h"
#include "atomiverse/core/real.h"
#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/comm/domain.h"
#include "atomiverse/cpu/component/constraint_protocol.h"
#include "atomiverse/cpu/component/halo_participant.h"
#include "atomiverse/cpu/component/neighbor_protocol.h"
#include "atomiverse/cpu/component/pair_evaluator.h"
#include "atomiverse/cpu/driver/backend_traits.h"
#include "atomiverse/cpu/driver/dof.h"
#include "atomiverse/cpu/driver/loop_handler.h"
#include "atomiverse/cpu/recorder/thermo.h"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/driver.h"
#include "atomiverse/facade/iteration_status.h"
#include "atomiverse/facade/loop.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {

/// @brief Phase traits for each driver category.
///
/// Maps PhaseEnum → the eval phase (where pair/bonded forces run) and
/// the end-of-step phase (where thermo output runs). Specializations
/// below for md, minimizer, and evaluate drivers.
template <typename PhaseEnum>
struct CategoryTraits;

/// `neighbor_phase` is where the neighbor-list rebuild binds. MD has a
/// dedicated `neighbor` phase (LAMMPS-granular skeleton); the other drivers
/// have no separate phase, so it folds into `eval_phase` just before the
/// force evaluators (order -10).
template <>
struct CategoryTraits<ATV::md::Phase> {
    static constexpr auto eval_phase = ATV::md::Phase::force;
    static constexpr auto neighbor_phase = ATV::md::Phase::neighbor;
    static constexpr auto end_phase = ATV::md::Phase::end_of_step;
};

template <>
struct CategoryTraits<ATV::minimizer::Phase> {
    static constexpr auto eval_phase = ATV::minimizer::Phase::force;
    static constexpr auto neighbor_phase = ATV::minimizer::Phase::force;
    static constexpr auto end_phase = ATV::minimizer::Phase::end_of_step;
};

template <>
struct CategoryTraits<ATV::evaluate::Phase> {
    static constexpr auto eval_phase = ATV::evaluate::Phase::evaluate;
    static constexpr auto neighbor_phase = ATV::evaluate::Phase::evaluate;
    static constexpr auto end_phase = ATV::evaluate::Phase::finalize;
};

template <>
struct CategoryTraits<ATV::vib::Phase> {
    static constexpr auto eval_phase = ATV::vib::Phase::force;
    static constexpr auto neighbor_phase = ATV::vib::Phase::force;
    static constexpr auto end_phase = ATV::vib::Phase::end_of_step;
};

/// Monte Carlo (Sampling category): the energy evaluator + neighbor list bind
/// in the `energy` phase (the trial-state recompute), and thermo/observers in
/// `end_of_step`. There is no separate neighbor phase, so it folds into
/// `energy` just ahead of the evaluator (order -10, before ZeroEnergy at -999
/// runs first and the evaluator at 0 last).
template <>
struct CategoryTraits<ATV::sampling::Phase> {
    static constexpr auto eval_phase = ATV::sampling::Phase::energy;
    static constexpr auto neighbor_phase = ATV::sampling::Phase::energy;
    static constexpr auto end_phase = ATV::sampling::Phase::end_of_step;
};

/// Path category: force evaluators bind at `force` (as in Minimize), the
/// neighbor list folds in just ahead of them, and thermo/observers at
/// `end_of_step`. Without this specialization the categorical setters
/// inherited from CPUDriver do not instantiate — the whole reason the
/// specialization exists as a separate step from adding the enum.
template <>
struct CategoryTraits<ATV::path::Phase> {
    static constexpr auto eval_phase = ATV::path::Phase::force;
    static constexpr auto neighbor_phase = ATV::path::Phase::force;
    static constexpr auto end_phase = ATV::path::Phase::end_of_step;
};

/// @brief CPU-backend driver layer — shared by all three CPU category bases.
template <typename PhaseEnum>
class CPUDriver : public Driver<PhaseEnum, cpu::device_tag> {
  public:
    using Base = Driver<PhaseEnum, cpu::device_tag>;
    using Traits = CategoryTraits<PhaseEnum>;

    /// @brief Compile the Plan after installing infrastructure and
    /// validating the required role + capability set.
    ///
    /// Order: `install_defaults()` first (so auto-injected components are
    /// visible to validation and the compiled Plan), `validate_components()`
    /// second (so the auto-injected set counts toward role/capability
    /// requirements), `Base::build()` last.
    void build() {
        install_defaults();
        bind_halo_participants_();
        validate_components();
        // Capture the Constraint components NOW, while `pending_` still holds them:
        // `Base::build()` consumes the pending list into the compiled Plan, so a scan
        // in setup_pass() would find an empty list and publish Ndof = 3N-3 with zero
        // constraints. That is not hypothetical — it shipped, and the rigid-water
        // regression caught it: the geometry was rigid but the thermostat drove the
        // system to the UNCONSTRAINED kinetic target, landing at 447 K against a
        // 300 K setpoint. The count itself still has to wait for setup_pass, because
        // it depends on the frame, which is not bound yet.
        // Dedup by pointer: `for_each_pending` walks BINDINGS, not unique components,
        // and FixRattle declares three (snapshot / positions / velocities). Counting
        // it once per binding tripled the constraint count and made the driver throw
        // "no degrees of freedom left" on a perfectly ordinary water box.
        // `Backend::collect_component_memory` dedups the same way, for the same reason.
        constraint_providers_.clear();
        std::unordered_set<const cpu::ConstraintProvider*> seen;
        for_each_pending([&](const cpu::ComponentBase& c) {
            if (const cpu::ConstraintProvider* cp = c.as_constraint()) {
                if (seen.insert(cp).second) constraint_providers_.push_back(cp);
            }
        });
        Base::build();
    }

    /// @brief Enable/disable opt-in per-component timing for this run.
    ///
    /// Eagerly creates the backend so the flag survives into `run()`; call before
    /// `run()`. Also enabled by the `ATV_TIMER` env var.
    ///
    /// Lives on the category base, not on one concrete driver: the timing registry
    /// is Backend infrastructure and every category has one. (It sat on
    /// `cpu::md::MD` alone until an out-of-library regression case needed to quiet
    /// a chain of `Evaluate` single points and found it missing — the second caller
    /// that justifies hoisting it.)
    /// @param on True to record per-component timings.
    void enable_timing(bool on = true) {
        if (!this->backend_)
            this->backend_ = std::make_unique<typename Base::BackendT>();
        this->backend_->enable_timing(on);
    }

    /// @brief Reset the backend's accumulated step timing (e.g. per stage).
    void reset_timing() {
        if (this->backend_) this->backend_->reset_timing();
    }

    /// @brief Cap the OpenMP thread count for this driver.
    ///
    /// The engine otherwise takes whatever `omp_get_max_threads()` reports, which on a
    /// fat node is every core — and for a small system the fork/join overhead then
    /// dwarfs the physics by orders of magnitude. Measured on a 24-atom rigid-water
    /// test on a 48-core node: 3974 ms at 48 threads, 5 ms at 4, 2 ms at 1. The
    /// physics costs microseconds; the threads cost seconds.
    ///
    /// Until now the only way out was the `OMP_NUM_THREADS` environment variable, which
    /// an in-process consumer (a library user, a test, a regression case) cannot set for
    /// itself after startup. Call before `build()`.
    ///
    /// @param n Thread count; values < 1 are ignored.
    void set_threads(int n) {
        if (n < 1) return;
        if (!this->backend_)
            this->backend_ = std::make_unique<typename Base::BackendT>();
        this->backend_->set_threads(n);
    }

  protected:
    /// @brief Subclass hook — install opt-in infrastructure components
    /// derived from declared user state (e.g. `NoseHooverNPT::reads()
    /// == Stress` → inject `Virial`; pair evaluator + multi-rank
    /// build → inject MPI halo exchange). Default no-op. Subclasses
    /// override to apply category-specific rules. Called from `build()`
    /// before `validate_components()` so the auto-injected set counts
    /// toward validation.
    /// @see .claude/specs/driver-infra-auto-injection.md
    virtual void install_defaults() {}

    /// @brief Subclass hook — validate the component set and throw on
    /// missing roles / capabilities. Called *before* the Plan is built,
    /// so `this->pending_` still carries the registered components
    /// (including any auto-injected by `install_defaults()`).
    virtual void validate_components() const = 0;

    /// @brief Walk every (component, phase, order) pending binding.
    /// Components with multiple `default_bindings()` entries appear once
    /// per binding — duplicate-tolerant predicates expected.
    template <typename Fn>
    void for_each_pending(Fn&& fn) const {
        for (const auto& pb : this->pending_) {
            fn(*pb.component);
        }
    }

    /// @brief Mutable variant (named per the `cuda_driver.cuh`
    /// `for_each_pending_mut` precedent) — the halo-protocol probes
    /// (`halo_participant()`) are non-const because the bind pass hands
    /// endpoints a participant list to store.
    template <typename Fn>
    void for_each_pending_mut(Fn&& fn) {
        for (auto& pb : this->pending_) {
            fn(*pb.component);
        }
    }

    /// @brief Build-time halo wiring: collect every halo participant among
    /// the pending components (deduplicated — multi-binding components
    /// appear once) and hand the list to each exchange endpoint via
    /// `bind_participants()`. Runs after `install_defaults()` (so
    /// auto-injected endpoints are covered) and before component `build()`
    /// (endpoints flatten field descriptors lazily at first launch).
    void bind_halo_participants_() {
        std::vector<const cpu::HaloParticipant*> participants;
        std::vector<cpu::HaloParticipant*> endpoints;
        std::unordered_set<const cpu::ComponentBase*> seen;
        for (auto& pb : this->pending_) {
            cpu::ComponentBase* c = pb.component;
            if (!seen.insert(c).second)
                continue;
            cpu::HaloParticipant* hp = c->halo_participant();
            if (!hp)
                continue;
            participants.push_back(hp);
            if (hp->halo_endpoint() != cpu::HaloParticipant::Endpoint::none)
                endpoints.push_back(hp);
        }
        for (auto* e : endpoints)
            e->bind_participants(participants);
    }

    // ── Driver loop helpers (protected — subclasses only) ──────────

    /// @brief One-time setup — backend init, plan load, set frame, build components.
    /// Resets outer_status before the run so a stale flag from a prior
    /// invocation cannot leak into this one.
    void setup_pass(cpu::Frame& frame) {
        if (!this->backend_)
            this->backend_ = std::make_unique<typename Base::BackendT>();
        this->backend_->init();
        this->backend_->load(this->plan_);
        this->backend_->set_frame(frame);
        publish_dof_(frame);
        // Point the Context at THIS loop's status + metrics. The loop owns them; the
        // Context only borrows. A nested driver will re-point them at its own and
        // restore ours on exit, which is what makes a nested loop physically unable to
        // stop its parent.
        this->backend_->context().loop.status = &status_;
        this->backend_->context().loop.metrics = &metrics_;
        this->backend_->build_components();
        status_.reset();
        metrics_ = LoopMetrics{};
    }

    /// @brief Single iteration through the compiled Plan.
    void step_once(int iter_index) { this->backend_->step(iter_index); }

    /// @brief Emit the built-in end-of-run summary — the LAMMPS-style time +
    /// memory breakdown every driver reports. Called by the category driver's
    /// `run()` after the loop (all ranks, in lockstep — the reduction is
    /// collective). Default no-op; concrete drivers override to render their
    /// own category taxonomy. @see cpu/md/timing_report.h
    virtual void emit_run_summary(std::string_view /*stage*/) {}

    /// @brief Tear down — destroy components + shut down backend.
    void teardown_pass() {
        this->backend_->destroy_components();
        this->backend_->shutdown();
    }

    /// @brief This loop's iteration status. **Owned here, not on the Context.**
    const IterationStatus& outer_status() const { return status_; }

    /// @brief Mutable handle to this loop's status (for the per-iteration `reset()`).
    IterationStatus& mutable_outer_status() { return status_; }

    /// @brief This loop's metrics — what its components measured, for handlers to judge.
    LoopMetrics& metrics() { return metrics_; }

    /// @brief Ask every handler about the iteration just completed, and reduce.
    ///
    /// Composition: **Stop > Converged > Continue** — the most severe verdict wins, so a
    /// run that has gone numerically insane cannot report success by also happening to
    /// satisfy a force tolerance.
    ///
    /// The verdict is all-reduced across ranks HERE, once, rather than in each handler.
    /// Ranks that disagree about when to stop deadlock at the next collective, and that
    /// is not a hazard a handler author should have to know about. `Verdict` is ordered
    /// by severity precisely so the reduction is a `MAX`.
    ///
    /// @param iteration Index of the iteration just completed.
    /// @return The reduced verdict; also folds `stop_requested` (which a *component* may
    ///         set — components request stops, they never declare convergence).
    Verdict poll_handlers_(int iteration) {
        const cpu::LoopState st{iteration, this->backend_->frame(),
                                this->backend_->context(), metrics_};
        Verdict v = status_.stop_requested ? Verdict::Stop : Verdict::Continue;
        for (auto& h : handlers_) {
            const Verdict hv = h->after_iteration(st);
            if (hv != Verdict::Continue && worse_of(v, hv) == hv && v != hv) {
                status_.reason = h->name();
            }
            v = worse_of(v, hv);
        }
        v = reduce_verdict_(v);
        if (v == Verdict::Converged) status_.converged = true;
        if (v == Verdict::Stop) status_.stop_requested = true;
        return v;
    }

    /// @brief Convenience: throw a clear exception when a required role or
    /// capability is absent at `validate_components()` time.
    [[noreturn]] static void throw_missing(const char* what) {
        throw std::runtime_error(std::string("driver build: missing required ")
                                 + what);
    }

  public:
    // ─── Categorical setter API ─────────────────────────────────────
    //
    // These are the LAMMPS-shaped setters shared across MD, Minimize,
    // and Evaluate. Each concrete driver inherits them; driver-specific
    // additions (set_integrator, set_optimizer, etc.) are declared in
    // the concrete class.

    /// @brief Configure the pair-cutoff force evaluator (LAMMPS `pair_style`).
    /// Singular: re-call with same class idempotent, different class throws.
    template <typename P>
    void set_pair(P&& p) {
        using D = std::decay_t<P>;
        static_assert(std::is_base_of_v<ATV::cpu::PairEvaluator, D>,
                      "set_pair requires a cpu::PairEvaluator subclass");
        set_singular_at_(std::forward<P>(p), pair_set_, "set_pair",
                         Traits::eval_phase, 0);
    }

    /// @brief Configure the bond force evaluator (LAMMPS `bond_style`).
    template <typename B>
    void set_bond(B&& b) {
        set_singular_at_(std::forward<B>(b), bond_set_, "set_bond",
                         Traits::eval_phase, 0);
    }

    /// @brief Configure the angle force evaluator (LAMMPS `angle_style`).
    template <typename A>
    void set_angle(A&& a) {
        set_singular_at_(std::forward<A>(a), angle_set_, "set_angle",
                         Traits::eval_phase, 0);
    }

    /// @brief Configure the dihedral force evaluator (LAMMPS `dihedral_style`).
    template <typename D>
    void set_dihedral(D&& d) {
        set_singular_at_(std::forward<D>(d), dihedral_set_, "set_dihedral",
                         Traits::eval_phase, 0);
    }

    /// @brief Configure the neighbor list (LAMMPS `neighbor` + `neigh_modify`).
    /// Bound at the driver's `neighbor_phase` @ -10: for MD that is the
    /// dedicated `neighbor` phase (runs after ghost build, before `force`); for
    /// the other drivers it folds into `eval_phase` just ahead of the force
    /// evaluators.
    template <typename N>
    void set_neighbor(N&& n_comp) {
        set_singular_at_(std::forward<N>(n_comp), neighbor_set_, "set_neighbor",
                         Traits::neighbor_phase, -10);
    }

    /// @brief Add a derived-quantity observer (LAMMPS `compute`). Plural.
    template <typename Obs>
    void set_observer(Obs&& obs) {
        this->add(std::forward<Obs>(obs));
    }

    /// @brief Add a recorder / I/O sink (LAMMPS `dump`). Plural.
    template <typename R>
    void set_record(R&& r) {
        this->add(std::forward<R>(r));
    }

    /// @brief Attach a termination policy. Plural — handlers compose.
    ///
    /// `opt.add_handler(ForceConverged{0.01});`
    /// `sim.add_handler(StopOnNaN{}).add_handler(StopWhenTooHot{2000});`
    ///
    /// This is the seam that used to be a class hierarchy. @see cpu/driver/loop_handler.h
    /// @param h A LoopHandler subclass, moved in.
    /// @return Reference to self, for chaining.
    template <typename H>
    CPUDriver& add_handler(H&& h) {
        using D = std::decay_t<H>;
        static_assert(std::is_base_of_v<cpu::LoopHandler, D>,
                      "add_handler requires a cpu::LoopHandler subclass");
        handlers_.push_back(std::make_unique<D>(std::forward<H>(h)));
        return *this;
    }

    /// @brief Whether any handler is attached (a preset injects a default if not).
    bool has_handler() const { return !handlers_.empty(); }

    /// @brief Did this loop terminate because a handler declared success?
    ///
    /// **The loop knows; the algorithm does not.** `LBFGS::converged()` and
    /// `CG::converged()` are deleted: an optimizer that publishes metrics and does not
    /// judge them has no business claiming to know the verdict, and asking it was how
    /// the policy got welded to the algorithm in the first place.
    bool converged() const { return status_.converged; }

    /// @brief Which handler ended the loop (empty if it simply ran out of budget).
    const std::string& stop_reason() const { return status_.reason; }

    /// @brief Configure the thermo print cadence (LAMMPS `thermo N`).
    /// Singular: re-call with a different cadence throws.
    /// Binds `cpu::Thermo{every}` at the driver's end-of-step phase @ +20.
    void set_thermo(int every) {
        if (thermo_every_) {
            if (*thermo_every_ == every)
                return;
            throw std::invalid_argument("set_thermo already configured with "
                                        + std::to_string(*thermo_every_)
                                        + ", rejected re-call with "
                                        + std::to_string(every));
        }
        thermo_every_ = every;
        this->on(Traits::end_phase).add(ATV::cpu::Thermo{every}, +20);
    }

    // ─── Singular-setter helpers (used by the categorical setters above) ──

    /// @brief Singular setter using the component's `default_bindings()`.
    /// Idempotent on same-name re-call; throws on different-name re-call.
    template <typename C>
    void set_singular_default_(C&& comp, std::optional<std::string>& slot,
                               const char* setter_name) {
        const std::string n = comp.name();
        if (slot) {
            if (*slot == n)
                return;
            throw std::invalid_argument(
                std::string(setter_name) + " already configured with '" + *slot
                + "', rejected re-call with '" + n + "'");
        }
        slot = n;
        this->add(std::forward<C>(comp));
    }

    /// @brief Singular setter at an explicit phase / order. Used for
    /// components without `default_bindings()` (NeighborList) or where
    /// the categorical API pins a specific slot.
    template <typename C>
    void set_singular_at_(C&& comp, std::optional<std::string>& slot,
                          const char* setter_name, PhaseEnum phase, int order) {
        const std::string n = comp.name();
        if (slot) {
            if (*slot == n)
                return;
            throw std::invalid_argument(
                std::string(setter_name) + " already configured with '" + *slot
                + "', rejected re-call with '" + n + "'");
        }
        slot = n;
        this->on(phase).add(std::forward<C>(comp), order);
    }

  protected:
    // ── Infrastructure scan (subclass install_defaults() only) ──────────

    /// @brief Aggregate of pending-component facts that `install_defaults()`
    /// uses to decide what to auto-inject and what to throw about. Single
    /// pass over `pending_` populates everything.
    struct InfraScan {
        bool needs_stress = false; ///< Some pending component reads Stress.
        bool has_virial_stress =
            false; ///< Canonical Stress provider already there.
        bool has_neighbor_list = false; ///< NeighborList already there.
        bool has_constraint = false;    ///< A Constraint component is registered.
        bool has_halo_send = false;     ///< A HaloSend is already registered.
        bool has_halo_recv = false;     ///< A HaloRecv is already registered.
        std::string neighbor_consumer_name; ///< First `NeighborConsumer`.
        std::optional<Real>
            neighbor_consumer_cutoff; ///< Set when it declares a cutoff > 0.
    };

    /// @brief Walk the pending list once, collect everything install_defaults
    /// rules might need.
    ///
    /// Keyed entirely on capability-typed probes (`reads()`,
    /// `capabilities()`, `role()`, the `as_neighbor_consumer()` /
    /// `as_neighbor_provider()` protocol accessors, and the
    /// `halo_participant()`/`halo_endpoint()` protocol accessor) — no
    /// `name()` string matching and no RTTI cast, so user-supplied
    /// neighbor structures / stress providers / halo endpoints participate
    /// in auto-injection whatever their `name()`. Same pattern as
    /// `cuda::md::MD::install_infra_`. Non-const because the halo probe is
    /// non-const (its other caller, the bind pass, mutates endpoints).
    /// @brief Compute and publish the global DOF count — the single owner.
    ///
    /// `Ndof = 3N - 3 - n_constraints`. Runs between `set_frame()` and
    /// `build_components()`, and that ordering is load-bearing: thermostats cache
    /// Ndof in their own `build()`, so a count published any later would be ignored
    /// by exactly the components that most need it.
    ///
    /// Until constraints arrived, five separate sites each hardcoded `3N - 3`
    /// (Thermo, InitVelocity, NoseHoover, NoseHooverNPT, FixCSVR). That formula is
    /// simply wrong once anything is constrained, and the failure is silent: rigid
    /// water reports T ~33 % low and every thermostat drives to the wrong setpoint.
    void publish_dof_(cpu::Frame& frame) {
        int n_con = 0;
        for (const cpu::ConstraintProvider* cp : constraint_providers_)
            n_con += cp->n_constraints(this->backend_->context(), frame);
        cpu::publish_global_dof(this->backend_->context(), frame, n_con);
    }

    /// Constraint components, captured at `build()` while the pending list still
    /// exists. Non-owning: the Plan owns the components.
    std::vector<const cpu::ConstraintProvider*> constraint_providers_;

    InfraScan scan_pending_infra_() {
        InfraScan s;
        for_each_pending_mut([&](cpu::ComponentBase& c) {
            if (has_capability(c.reads(), Capability::Stress))
                s.needs_stress = true;
            if (c.role() == ComponentRole::Observer
                && has_capability(c.capabilities(), Capability::Stress))
                s.has_virial_stress = true;
            if (c.as_neighbor_provider())
                s.has_neighbor_list = true;
            if (c.role() == ComponentRole::Constraint)
                s.has_constraint = true;
            if (cpu::HaloParticipant* hp = c.halo_participant()) {
                if (hp->halo_endpoint() == cpu::HaloParticipant::Endpoint::send)
                    s.has_halo_send = true;
                if (hp->halo_endpoint() == cpu::HaloParticipant::Endpoint::recv)
                    s.has_halo_recv = true;
            }
            if (const cpu::NeighborConsumer* nc = c.as_neighbor_consumer()) {
                if (s.neighbor_consumer_name.empty())
                    s.neighbor_consumer_name = c.name();
                // MAX over every consumer, not the first one that happens to be
                // scanned. The halo has to satisfy all of them at once: a plan with
                // PairLJCut (8 A) and CoulLong (12 A) used to size its ghost shell
                // from whichever landed first in the pending list. And a many-body
                // style asks for MORE than its interaction cutoff — EAM needs a
                // ghost-of-ghost shell of 2 r_c to complete rho on ghosts (see
                // NeighborConsumer::halo_cutoff).
                const Real hc = nc->halo_cutoff();
                if (hc > Real(0)
                    && (!s.neighbor_consumer_cutoff || hc > *s.neighbor_consumer_cutoff))
                    s.neighbor_consumer_cutoff = hc;
            }
        });
        return s;
    }

    /// @brief Shared missing-NL check. `setter_help` is the
    /// driver-specific instruction printed in the throw message
    /// (`set_pair` is the same, but the variable name differs:
    /// `sim` / `opt` / `eval`).
    void throw_if_missing_neighbor_list_(const InfraScan& s,
                                         const char* driver_name,
                                         const char* setter_help) const {
        if (!s.neighbor_consumer_name.empty() && !s.has_neighbor_list) {
            throw std::runtime_error(std::string(driver_name)
                                     + "::build(): component '"
                                     + s.neighbor_consumer_name
                                     + "' is a NeighborConsumer but no "
                                       "neighbor list was registered. "
                                     + setter_help);
        }
    }

    /// @brief Reduce a verdict across ranks (MAX by severity). @see poll_handlers_
    ///
    /// Severity is the enum's integer value, so "the most severe verdict any rank
    /// reached" is a plain MPI_MAX. Every rank therefore leaves the loop on the same
    /// iteration — the alternative is a rank that keeps looping into a collective nobody
    /// else will ever enter.
    Verdict reduce_verdict_(Verdict v) const {
#ifdef ATV_HAS_MPI
        const auto* dd = this->backend_->context().md.dd;
        if (this->backend_->context().mpi_size > 1 && dd
            && dd->cart_comm != MPI_COMM_NULL) {
            int sev = static_cast<int>(v);
            MPI_Allreduce(MPI_IN_PLACE, &sev, 1, MPI_INT, MPI_MAX, dd->cart_comm);
            return static_cast<Verdict>(sev);
        }
#endif
        return v;
    }

  protected:
    /// This loop's own status and metrics. **Not on the Context** — see LoopServices.
    IterationStatus status_{};
    LoopMetrics metrics_{};
    std::vector<std::unique_ptr<cpu::LoopHandler>> handlers_;

  private:
    // ── Categorical setter slots ────────────────────────────────────

    std::optional<std::string> pair_set_;
    std::optional<std::string> bond_set_;
    std::optional<std::string> angle_set_;
    std::optional<std::string> dihedral_set_;
    std::optional<std::string> neighbor_set_;
    std::optional<int> thermo_every_;
};

} // namespace cpu
} // namespace ATV
