/// @file cpu/md/md.h
/// @brief cpu::md::MD — concrete time-advance driver under DynamicsDriver.
#pragma once

#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "atomiverse/core/device.h"
#include "atomiverse/cpu/backend_traits.h"
#include "atomiverse/cpu/dynamics_driver.h"
#include "atomiverse/cpu/snapshot_publisher.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/md/mpi/atom_migration.h"
#include "atomiverse/cpu/md/mpi/halo_exchange.h"
#include "atomiverse/cpu/md/mpi/reduce.h"
#include "atomiverse/cpu/pair_evaluator.h"
#include "atomiverse/cpu/md/virial_stress.h"
#include "atomiverse/cpu/observer/thermo_component.h"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief CPU MD simulation driver.
///
/// Convenience class that auto-adds always-on infrastructure components
/// in the constructor, exposes a LAMMPS-shaped categorical setter API
/// (`set_pair`, `set_bond`, `set_neighbor`, `set_integrator`,
/// `set_thermostat`, `set_barostat`, `set_observer`, `set_record`,
/// `set_thermo`), and runs `install_infra()` at `build()` time to
/// derive opt-in infrastructure from declared component requirements:
///
///   - Stress capability dependency → auto-injects `VirialStress`
///     (deduped against an explicitly-registered `VirialStress`).
///   - `#ifdef ATV_BUILD_MPI` and a pair evaluator was set →
///     auto-injects `HaloExchangeForward{cutoff}` and
///     `HaloExchangeReverse{}` at the documented orders.
///   - A pair evaluator with `needs_neighbor_list()==true` was
///     registered without a corresponding `set_neighbor()` →
///     `build()` throws naming the offending evaluator.
///
/// `add()` and `on(phase).add()` from the underlying `Driver<>`
/// remain public — they are the documented escape hatch for
/// user-defined hooks (custom observers, recorders, etc.) that do
/// not yet have a categorical setter. Mixing categorical + raw calls
/// is supported and produces a single coherent Plan.
///
/// End-of-step ordering: `ZeroKinetic(-1000) → KineticEnergy(-100) →
/// ReduceEnergy(-50) → SnapshotPublisher(+10) → user observers`.
/// Categorical thermo at `set_thermo(N)` binds `cpu::Thermo{N}` at
/// `end_of_step @ +20`. pre_force ordering includes `AtomMigration(-30)`,
/// which follows the LAMMPS reneighbor cadence by attempting exchange
/// every 10 steps.
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
        // Always-on zero / kinetic infrastructure.
        add(ZeroForces{});
        add(ZeroEnergy{});
        add(ZeroStress{});   // pair-virial accumulator must start at 0
        add(ZeroKinetic{});
        add(KineticEnergy{});

        // MPI infrastructure (single-rank no-ops in Phase 1).
        add(mpi::AtomMigration{10});
        add(mpi::ReduceEnergy{});

        // Async-analysis publisher -- shared infrastructure, explicitly
        // bound to *this driver's* end_of_step @ +10. Subscribers
        // inheriting `cpu::md::AsyncObserver` default to +20 on the same
        // phase and see this publisher's slot.
        this->on(ATV::md::Phase::end_of_step)
            .add(cpu::SnapshotPublisher{}, 10);
    }

    // ─── Categorical setters (LAMMPS-shaped) ──────────────────────────

    /// @brief Configure the pair-cutoff force evaluator
    /// (LAMMPS analogue: `pair_style …`). Singular: re-call with the
    /// same component class is idempotent; re-call with a different
    /// class throws `std::invalid_argument`.
    /// @tparam P Subclass of `cpu::md::PairEvaluator`.
    template <typename P>
    MD& set_pair(P&& p) {
        using D = std::decay_t<P>;
        static_assert(std::is_base_of_v<ATV::cpu::PairEvaluator, D>,
                      "set_pair requires a cpu::PairEvaluator subclass");
        this->set_singular_at_(std::forward<P>(p), pair_set_, "set_pair",
                               ATV::md::Phase::force, 0);
        return *this;
    }

    /// @brief Configure the bond force evaluator (LAMMPS `bond_style`).
    template <typename B>
    MD& set_bond(B&& b) {
        this->set_singular_at_(std::forward<B>(b), bond_set_, "set_bond",
                               ATV::md::Phase::force, 0);
        return *this;
    }

    /// @brief Configure the angle force evaluator (LAMMPS `angle_style`).
    template <typename A>
    MD& set_angle(A&& a) {
        this->set_singular_at_(std::forward<A>(a), angle_set_, "set_angle",
                               ATV::md::Phase::force, 0);
        return *this;
    }

    /// @brief Configure the dihedral force evaluator (LAMMPS `dihedral_style`).
    template <typename D>
    MD& set_dihedral(D&& d) {
        this->set_singular_at_(std::forward<D>(d), dihedral_set_, "set_dihedral",
                               ATV::md::Phase::force, 0);
        return *this;
    }

    /// @brief Configure the neighbor list (LAMMPS `neighbor` + `neigh_modify`).
    /// Bound at `force @ -10` so it runs before pair / bonded force evaluators.
    template <typename N>
    MD& set_neighbor(N&& n_comp) {
        this->set_singular_at_(std::forward<N>(n_comp), neighbor_set_, "set_neighbor",
                               ATV::md::Phase::force, -10);
        return *this;
    }

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

    /// @brief Add a barostat (LAMMPS `fix … press/berendsen`).
    /// Plural.
    template <typename Bs>
    MD& set_barostat(Bs&& b) {
        this->add(std::forward<Bs>(b));
        return *this;
    }

    /// @brief Add a derived-quantity observer (LAMMPS `compute …`).
    /// Plural.
    template <typename O>
    MD& set_observer(O&& o) {
        this->add(std::forward<O>(o));
        return *this;
    }

    /// @brief Add a recorder / I/O sink (LAMMPS `dump …`). Plural.
    template <typename R>
    MD& set_record(R&& r) {
        this->add(std::forward<R>(r));
        return *this;
    }

    /// @brief Configure the thermo print cadence (LAMMPS `thermo N`).
    /// Singular: re-call with a different cadence throws.
    /// Binds `cpu::Thermo{every}` at `end_of_step @ +20`.
    MD& set_thermo(int every) {
        if (thermo_every_) {
            if (*thermo_every_ == every) return *this;
            throw std::invalid_argument(
                "set_thermo already configured with " +
                std::to_string(*thermo_every_) +
                ", rejected re-call with " + std::to_string(every));
        }
        thermo_every_ = every;
        this->on(ATV::md::Phase::end_of_step)
             .add(ATV::cpu::Thermo{every}, +20);
        return *this;
    }

protected:
    /// @brief Apply the three install_infra rules:
    ///   - Stress capability dependency → canonical `VirialStress`.
    ///   - `#ifdef ATV_BUILD_MPI` halo exchange when a pair was set.
    ///   - Missing-NL throw for any pair-cutoff evaluator without a
    ///     corresponding neighbor list.
    ///
    /// Stress dedup is by canonical-provider name, not by capability
    /// bit: LJ / bonded forces also write `Capability::Stress` (the
    /// pair / bond *virial* term), but the *kinetic* term
    /// `Σ m·v_α·v_β` that NPT pressure feedback needs is
    /// `VirialStress`'s sole job.
    void install_infra() override {
        const auto s = this->scan_pending_infra_();

        if (s.needs_stress && !s.has_virial_stress)
            this->add(VirialStress{});

#ifdef ATV_BUILD_MPI
        if (s.pair_evaluator_cutoff) {
            this->add(mpi::HaloExchangeForward{static_cast<double>(*s.pair_evaluator_cutoff)});
            this->add(mpi::HaloExchangeReverse{});
        }
#endif

        this->throw_if_missing_neighbor_list_(
            s, "MD",
            "Add `sim.set_neighbor(cpu::md::NeighborList<>{cutoff, skin})`.");
    }

private:
    std::optional<std::string> pair_set_;
    std::optional<std::string> bond_set_;
    std::optional<std::string> angle_set_;
    std::optional<std::string> dihedral_set_;
    std::optional<std::string> neighbor_set_;
    std::optional<std::string> integrator_set_;
    std::optional<int>         thermo_every_;
};

} // namespace md
} // namespace cpu
} // namespace ATV
