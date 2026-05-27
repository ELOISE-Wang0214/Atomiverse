/// @file cpu/md/md.h
/// @brief cpu::md::MD — concrete time-advance driver under DynamicsDriver.
#pragma once

#include <string>
#include <type_traits>

#include "atomiverse/cpu/backend.h"
#include "atomiverse/cpu/dynamics_driver.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/md/mpi/atom_migration.h"
#include "atomiverse/cpu/md/mpi/domain.h"
#include "atomiverse/cpu/md/mpi/halo_exchange.h"
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
/// are inherited from CPUDriver), and runs `install_infra()` at `build()`
/// time to derive opt-in infrastructure from declared component requirements:
///
///   - Stress capability dependency → auto-injects `VirialStress`
///     (deduped against an explicitly-registered `VirialStress`).
///   - `#ifdef ATV_HAS_MPI` and a pair evaluator was set →
///     auto-injects `HaloExchangeForward{cutoff}` and
///     `HaloExchangeReverse{}` at the documented orders.
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
/// End-of-step ordering: ZeroKinetic(-1000) → KineticEnergy(-100) →
/// ReduceEnergy(-50) → SnapshotPublisher(+10) → user observers.
/// Categorical thermo at `set_thermo(N)` binds `cpu::Thermo{N}` at
/// `end_of_step @ +20`. pre_force ordering includes AtomMigration(-30).
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

        // MPI infrastructure (single-rank no-ops in Phase 1).
        add(mpi::AtomMigration{10});
        add(mpi::ReduceEnergy{});

        // Async-analysis publisher — shared infrastructure, explicitly
        // bound to *this driver's* end_of_step @ +10. Subscribers
        // inheriting `cpu::md::AsyncObserver` default to +20 on the same
        // phase and see this publisher's slot.
        this->on(ATV::md::Phase::end_of_step).add(cpu::SnapshotPublisher{}, 10);
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

  protected:
    void install_infra() override {
        const auto s = this->scan_pending_infra_();

        if (s.needs_stress && !s.has_virial_stress)
            this->add(VirialStress{});

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
            // `HaloExchangeForward{halo_cut}` explicitly before the
            // pair evaluator — `install_infra` will see an existing
            // halo and skip the auto-inject.
            const double halo_cut =
                1.1 * static_cast<double>(*s.pair_evaluator_cutoff);
            this->add(mpi::HaloExchangeForward{halo_cut});
            this->add(mpi::HaloExchangeReverse{});
        }
#endif

        this->throw_if_missing_neighbor_list_(
            s, "MD",
            "Add `sim.set_neighbor(cpu::md::NeighborList<>{cutoff, skin})`.");
    }

  private:
    std::optional<std::string> integrator_set_;
};

} // namespace md
} // namespace cpu
} // namespace ATV
