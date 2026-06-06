/// @file cpu/md/md.h
/// @brief cpu::md::MD — concrete time-advance driver under DynamicsDriver.
#pragma once

#include <cstdlib>
#include <string>
#include <string_view>
#include <type_traits>

#include "atomiverse/cpu/backend.h"
#include "atomiverse/cpu/md/timing_report.h"
#include "atomiverse/cpu/dynamics_driver.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/comm/atom_migration.h"
#include "atomiverse/cpu/md/mpi/domain.h"
#include "atomiverse/cpu/comm/halo_exchange.h"
#include "atomiverse/cpu/md/mpi/reduce.h"
#include "atomiverse/cpu/md/virial_stress.h"
#include "atomiverse/cpu/pair_evaluator.h"
#include "atomiverse/cpu/snapshot_publisher.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief CPU MD simulation driver.
///
/// Convenience class that auto-adds always-on infrastructure components
/// in the constructor, exposes a LAMMPS-shaped categorical setter API
/// (set_pair, set_bond, set_neighbor, set_integrator, set_thermostat,
/// set_barostat, set_observer, set_record, set_thermo — the first six
/// are inherited from CPUDriver), and runs `install_defaults()` at `build()`
/// time to derive opt-in infrastructure from declared component requirements:
///
///   - Stress capability dependency → auto-injects `Virial`
///     (deduped against an explicitly-registered `Virial`).
///   - `#ifdef ATV_HAS_MPI` and a pair evaluator was set →
///     auto-injects `HaloSend{cutoff}` and
///     `HaloRecv{}` at the documented orders.
///   - A pair evaluator with `needs_neighbor_list()==true` was
///     registered without a corresponding `set_neighbor()` →
///     `build()` throws naming the offending evaluator.
///
/// Inheritance: MD → DynamicsDriver → CPUDriver<md::Phase> →
/// Driver<md::Phase, cpu::device_tag>.
/// CPUDriver provides the shared categorical setter API (set_pair,
/// set_bond, set_angle, set_dihedral, set_neighbor, set_observer,
/// set_record, set_thermo) and the slot variables. MD adds
/// driver-specific setters (set_integrator, set_thermostat, set_barostat).
///
/// `add()` and `on(phase).add()` from the underlying `Driver<>`
/// remain public — they are the documented escape hatch for
/// user-defined hooks (custom observers, recorders, etc.) that do
/// not yet have a categorical setter.
///
/// Phase skeleton (LAMMPS-granular; see
/// `.claude/notes/cpu-md-execution-model.md`): initial_integrate →
/// post_integrate → pre_exchange (Migrate) → pre_neighbor (HaloSend) →
/// neighbor (NeighborList) → pre_force (Zero{Forces,Energy,Stress}) → force
/// (pair + bonded) → post_force (HaloRecv, Virial) → final_integrate →
/// end_of_step. End-of-step ordering: ZeroKinetic(-1000) → KineticEnergy(-100)
/// → Reduce(-50) → Publish(+10) → user observers. Categorical thermo at
/// `set_thermo(N)` binds `cpu::Thermo{N}` at `end_of_step @ +20`.
///
/// Usage (LAMMPS-shaped):
/// @code{.cpp}
/// cpu::md::MD sim;
/// sim.set_pair(cpu::md::LJ{0.0104, 3.4, 6.8});
/// sim.set_neighbor(cpu::md::NeighborList<>{6.8, 0.3});
/// sim.set_integrator(cpu::md::VelocityVerlet{0.001});
/// sim.set_thermo(100);
/// sim.build();
/// sim.run(frame, 1000);
/// @endcode
class MD : public DynamicsDriver {
  public:
    MD() {
        // Always-on zero / kinetic infrastructure. Phase-agnostic
        // Zero* classes are bound explicitly to MD phases.
        this->on(ATV::md::Phase::pre_force)
            .add(ZeroForces{}, -1000)
            .add(ZeroEnergy{}, -999)
            .add(ZeroStress{}, -998);
        this->on(ATV::md::Phase::end_of_step)
            .add(ZeroKinetic{}, -1000)
            .add(KineticEnergy{}, -100);

        // MPI infrastructure (single-rank no-ops in Phase 1). `Reduce` is
        // injected in install_defaults() so set_scalar_reduce_every() (which
        // runs after construction) can configure its cadence.
        add(comm::Migrate{10});

        // Async-analysis publisher — shared infrastructure, explicitly
        // bound to *this driver's* end_of_step @ +10. Subscribers
        // inheriting `cpu::md::AsyncObserver` default to +20 on the same
        // phase and see this publisher's slot.
        this->on(ATV::md::Phase::end_of_step).add(cpu::Publish{}, 10);
    }

    // ─── MD-specific categorical setters ────────────────────────────
    //
    // set_pair, set_bond, set_angle, set_dihedral, set_neighbor,
    // set_observer, set_record, set_thermo are inherited from CPUDriver.

    /// @brief Configure the time integrator (LAMMPS `fix … nve|nvt|npt`).
    template <typename I>
    MD& set_integrator(I&& integ) {
        this->set_singular_default_(std::forward<I>(integ), integrator_set_,
                                    "set_integrator");
        return *this;
    }

    /// @brief Add a thermostat (LAMMPS `fix … langevin / temp/berendsen`).
    /// Plural: each call appends; insertion order preserved.
    template <typename T>
    MD& set_thermostat(T&& t) {
        this->add(std::forward<T>(t));
        return *this;
    }

    /// @brief Add a barostat (LAMMPS `fix … press/berendsen`). Plural.
    template <typename Bs>
    MD& set_barostat(Bs&& b) {
        this->add(std::forward<Bs>(b));
        return *this;
    }

    /// @brief Set the cadence of the global scalar `MPI_Allreduce` (energy /
    /// kinetic / stress). Default 1 = every step. A value > 1 fires the reduce
    /// only on step 0 and every Nth step thereafter — eliding the per-step
    /// latency-bound collective. ONLY valid when no per-step consumer reads the
    /// reduced scalars (NVE / per-atom Langevin, where the global energy/stress
    /// are pure diagnostics sampled every N steps). MUST stay 1 whenever a
    /// NoseHoover thermostat or barostat is active — those integrate the global
    /// temperature/pressure every step. Set N to align with the thermo /
    /// recorder sampling stride so sampled steps see globally-reduced scalars.
    MD& set_scalar_reduce_every(int every) {
        scalar_reduce_every_ = every < 1 ? 1 : every;
        return *this;
    }

    /// @brief Hand off the spatial `DomainDecomposition` from a previously-
    ///        torn-down MD instance.
    ///
    /// Call BEFORE `build()` to make the new MD reuse the partition state
    /// (`n_local`, `tag_to_local`, `sublo/subhi`, `cuts`) of the prior
    /// stage's DD. The companion `release_dd()` is invoked on the prior
    /// MD just before its destructor; the orchestrator (e.g.
    /// polyethylene_tension's per-stage loop) shuttles the unique_ptr
    /// through. Without this, every new MD ⇒ new Backend ⇒ new DD ⇒
    /// `first_partition=true` re-scan that misinterprets the previous
    /// stage's stale frame tail as owned atoms.
    void adopt_dd(std::unique_ptr<cpu::md::mpi::DomainDecomposition> dd) {
        using Base = Driver<ATV::md::Phase, cpu::device_tag>;
        if (!this->backend_)
            this->backend_ = std::make_unique<typename Base::BackendT>();
        this->backend_->adopt_dd(std::move(dd));
    }

    /// @brief Release the DD for the next stage to adopt.
    ///
    /// Call AFTER the last `run()` of this stage and BEFORE letting this
    /// MD destruct. Returns nullptr if no DD is owned.
    std::unique_ptr<cpu::md::mpi::DomainDecomposition> release_dd() {
        if (!this->backend_)
            return nullptr;
        return this->backend_->release_dd();
    }

    /// @brief Enable/disable opt-in per-component timing for this run. Eagerly
    /// creates the backend (like `adopt_dd`) so the flag survives into `run()`.
    /// Call before `run()`. (Also enabled by the `ATV_TIMER` env var.)
    void enable_timing(bool on = true) {
        if (!this->backend_)
            this->backend_ = std::make_unique<typename Base::BackendT>();
        this->backend_->enable_timing(on);
    }

    /// @brief Reset the backend's accumulated step timing (e.g. per stage).
    void reset_timing() {
        if (this->backend_)
            this->backend_->reset_timing();
    }

    /// @brief Raw access to the timing registry (null before the backend
    /// exists). @see core/timer.h
    const ATV::TimerRegistry* timing() const {
        return this->backend_ ? &this->backend_->timing() : nullptr;
    }

    /// @brief Cross-rank-reduced timing report for this run/stage.
    cpu::md::ReducedReport timing_report() const {
        return this->backend_ ? cpu::md::reduce_timing(this->backend_->timing())
                              : cpu::md::ReducedReport{};
    }

    /// @brief Emit the timing breakdown: human table (console, rank 0) and,
    /// if a Csv sink is given, a machine row (metrics file). Routed through
    /// the log module. @see cpu/md/timing_report.h
    void emit_timing(std::string_view stage,
                     ATV::report::Csv* csv = nullptr) const {
        if (!this->backend_)
            return;
        const auto rep = cpu::md::reduce_timing(this->backend_->timing());
        cpu::md::emit_timing_table(rep, stage);
        if (std::getenv("ATV_TIMER_VERBOSE"))
            cpu::md::emit_per_component(this->backend_->timing(), stage);
        if (csv)
            cpu::md::emit_timing_csv(rep, *csv, stage);
    }

  protected:
    void install_defaults() override {
        const auto s = this->scan_pending_infra_();

        // Global scalar reduction (single-rank no-op). Injected here rather
        // than in the ctor so set_scalar_reduce_every() can set the cadence.
        this->add(mpi::Reduce{scalar_reduce_every_});

        if (s.needs_stress && !s.has_virial_stress)
            this->add(Virial{});

#ifdef ATV_HAS_MPI
        if (s.pair_evaluator_cutoff) {
            // Halo cutoff must include a SKIN beyond the pair cutoff so
            // ghost atoms stay registered as they drift in/out of the
            // pair-cutoff zone between full halo rebuilds. Without the
            // skin band, atoms slip out of `tag_to_local` at each
            // `every_=10` rebuild and any cross-rank bond/angle/dihedral
            // involving them is silently dropped — observed under MSU
            // PE stage-2 NPT as a +1-3 eV per-rebuild jump in the MTK
            // extended Hamiltonian (E_bond drops by ~2.5 eV at rebuild
            // boundaries because ~16 chain bonds lose their ghost
            // partner each cycle). The 10% multiplier matches the
            // ~0.1 Å typical 10-step drift at PE-thermalized 500 K with
            // good margin while staying under `2× pair_cutoff` for the
            // 144-atom-in-48Å unit tests (those have slabs of 24 Å, so
            // a 13.65 Å halo would duplicate ghosts and overflow
            // `frame.atoms.size()`). Callers running with a custom NL
            // skin > 10% of cutoff should register
            // `HaloSend{halo_cut}` explicitly before the
            // pair evaluator — `install_defaults` will see an existing
            // halo and skip the auto-inject.
            // Skip the auto-inject when the caller already registered halo
            // components (custom cutoff / skin / reneighbor cadence). Auto-
            // injecting on top of an explicit halo would run the REVERSE
            // exchange (HaloRecv) twice, which double-applies the already-
            // reduced ghost forces to their owners — a per-step momentum leak
            // (Σ F ≠ 0) that pumps energy and blows the NVE trajectory up.
            // `reverse exchange must run exactly once` (see halo_exchange.h).
            const double halo_cut =
                1.1 * static_cast<double>(*s.pair_evaluator_cutoff);
            if (!s.has_halo_send)
                this->add(comm::HaloSend{halo_cut});
            if (!s.has_halo_recv)
                this->add(comm::HaloRecv{});
        }
#endif

        this->throw_if_missing_neighbor_list_(
            s, "MD",
            "Add `sim.set_neighbor(cpu::md::NeighborList<>{cutoff, skin})`.");
    }

  private:
    std::optional<std::string> integrator_set_;
    int scalar_reduce_every_ = 1;
};

} // namespace md
} // namespace cpu
} // namespace ATV
