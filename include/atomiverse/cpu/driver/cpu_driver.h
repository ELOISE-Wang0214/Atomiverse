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
///     │    └─ cpu::OptimizeDriver
///     └─ cuda::CUDADriver<PhaseEnum>       ← CUDA backend layer (planned —
///          ├─ cuda::EvaluateDriver           cuda-driver-taxonomy.md)
///          ├─ cuda::DynamicsDriver
///          └─ cuda::OptimizeDriver
/// @endcode
///
/// `CPUDriver<P>` inherits the underlying `Driver<P, cpu::device_tag>`
/// (which holds the Plan and Backend ownership) and adds the
/// CPU-specific machinery the three category bases share:
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
/// User code does not subclass this directly — the three concrete CPU
/// category bases (`cpu::EvaluateDriver` / `DynamicsDriver` /
/// `OptimizeDriver`) do.
#pragma once

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
#include "atomiverse/cpu/component/halo_participant.h"
#include "atomiverse/cpu/component/neighbor_protocol.h"
#include "atomiverse/cpu/component/pair_evaluator.h"
#include "atomiverse/cpu/driver/backend_traits.h"
#include "atomiverse/cpu/recorder/thermo.h"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/driver.h"
#include "atomiverse/facade/iteration_status.h"
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
        Base::build();
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
        this->backend_->build_components();
        mutable_outer_status().reset();
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

    /// @brief Read the outer-driver iteration status (set by ConvergenceCheck
    /// components during the most recent `step_once`).
    const IterationStatus& outer_status() const {
        return this->backend_->context().outer_status;
    }

    /// @brief Mutable handle to outer status (for per-iteration `reset()`).
    IterationStatus& mutable_outer_status() {
        return this->backend_->context().outer_status;
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
            if (cpu::HaloParticipant* hp = c.halo_participant()) {
                if (hp->halo_endpoint() == cpu::HaloParticipant::Endpoint::send)
                    s.has_halo_send = true;
                if (hp->halo_endpoint() == cpu::HaloParticipant::Endpoint::recv)
                    s.has_halo_recv = true;
            }
            if (const cpu::NeighborConsumer* nc = c.as_neighbor_consumer()) {
                if (s.neighbor_consumer_name.empty()) {
                    s.neighbor_consumer_name = c.name();
                    if (nc->neighbor_cutoff() > Real(0))
                        s.neighbor_consumer_cutoff = nc->neighbor_cutoff();
                }
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
