/// @file cuda/md/md.cuh
/// @brief cuda::md::MD -- molecular dynamics driver with LAMMPS-style phases.
///
/// `MD` is the concrete CUDA molecular-dynamics driver. It is a
/// `cuda::DynamicsDriver` — the time-advance outer category — and so
/// inherits that base's continuable `run(frame, n_steps)` (full step
/// budget, honours `outer_status.stop_requested`, ignores
/// `outer_status.converged`) and its build-time `validate_components()`
/// contract (≥ 1 `Integrator`-role component and ≥ 1 component carrying
/// `Capability::Force`).
///
/// The constructor adds only infrastructure components — it contributes
/// no runtime behaviour beyond the base. The Zero* initializers are
/// phase-agnostic (`cuda/md/init.cuh`); `MD` binds them explicitly to
/// `md::Phase` via `on(phase).add()`:
/// - ZeroForces at pre_force (order -1000) -- zeros fx, fy, fz.
/// - ZeroEnergy at pre_force (order -999) -- zeros potential energy.
/// - ZeroKinetic at end_of_step (order -1000) -- zeros kinetic energy.
/// - KineticEnergy at end_of_step (order -100) -- accumulates KE.
///
/// ## Execution model (shared with the CPU backend)
///
/// `MD` walks the same LAMMPS-granular `md::Phase` skeleton documented in
/// `docs/developer/md-execution-model.md` — every phase is always walked,
/// even when empty. This single-GPU driver binds work to `neighbor`
/// (NeighborList), `pre_force` (Zero*), `force` (pair/bonded),
/// `post_force` (Virial), `initial_integrate`/`final_integrate`
/// (integrator), and `end_of_step` (KE, thermo, observers). The
/// communication phases `pre_exchange` (atom migration) and
/// `pre_neighbor` (ghost-atom halo send), plus `post_integrate`, carry
/// no component on a single GPU — multi-GPU halo/migration is out of
/// scope per `msu-pe-validation.md` decision #1 — but the phases still
/// exist as documented "nothing happens here" slots, ready for a
/// future multi-GPU arm without a phase-skeleton change.
///
/// ## LAMMPS-shaped categorical setter API
///
/// `MD` exposes the same 11-method categorical surface as
/// `cpu::md::MD` (chain link `msu-pe-cuda-components-02`, Phase G.4):
/// `set_pair / set_bond / set_angle / set_dihedral / set_neighbor`
/// (singular, idempotent on same-class re-call, throw on
/// different-class re-call), `set_integrator` (singular via the
/// component's own `default_bindings()`), `set_thermostat /
/// set_barostat / set_observer / set_record` (plural — each call
/// appends, insertion order preserved), and `set_thermo(int every)`
/// (singular cadence; binds `cuda::Thermo{every}` at
/// `end_of_step @ +20`).
///
/// ## Auto-injection at build()
///
/// `MD::build()` runs an `install_infra_()` pass *before* the
/// `CUDADriver::build()` chain (so auto-injected components reach
/// `validate_components()`). Single-GPU subset of the CPU rules,
/// keyed on capability-typed virtuals (no string proxies, no
/// `dynamic_cast`):
///
/// - **Virial** — if any pending component declares
///   `(reads() & Capability::Stress) != None` and no pending component
///   already advertises `(capabilities() & Capability::Stress) != None`
///   with `role() == Observer`, a fresh `cuda::Virial{}` is
///   added at `post_force @ -50`. Covers both `cuda::md::Barostat`
///   (Berendsen, `role()==Barostat`) and `cuda::md::NoseHooverNPT`
///   (MTK NPT, `role()==Integrator`) uniformly via `reads()`.
/// - **NeighborList missing** — if any pending component reports
///   `needs_neighbor_list() == true` (pair-style evaluators do, via
///   `cuda::Pair`) without a registered `cuda::NeighborList`,
///   `build()` throws naming the offending evaluator.
///
/// Multi-GPU halo exchange auto-injection is **out of scope** for the
/// CUDA arm per `msu-pe-validation.md` fixed-scope decision #1.
///
/// The raw `add()` / `on(phase).add()` API remains public — it is the
/// documented escape hatch for components without a categorical slot,
/// and the additive setter surface must not break the existing
/// `tests/cuda/md/test_md.cu` paths that use it directly.
///
/// Usage:
/// @code{.cpp}
/// cuda::md::MD sim;
/// sim.set_pair(cuda::LJ{1.0, 1.0, 3.0});
/// sim.set_neighbor(cuda::NeighborList{3.0});
/// sim.set_integrator(cuda::md::VelocityVerlet{0.001});
/// sim.set_thermo(5);
/// sim.build();
/// sim.run(frame, 10000);
/// @endcode
///
/// @see cuda::DynamicsDriver for the inherited loop + validation
/// @see md::Phase for the LAMMPS-style phase sequence
/// @see cpu::md::MD for the CPU sibling driver with the same API surface
#pragma once

#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include "atomiverse/cuda/dynamics_driver.cuh"
#include "atomiverse/cuda/md/init.cuh"
#include "atomiverse/cuda/md/kinetic_energy.cuh"
#include "atomiverse/cuda/observer/stress.cuh"
#include "atomiverse/cuda/recorder/thermo.cuh"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cuda {
namespace md {

/// @brief CUDA MD driver — a `cuda::DynamicsDriver` with LAMMPS-style
///        phase sequence and the 11-method categorical setter API
///        that mirrors `cpu::md::MD`.
///
/// Inherits the time-advance loop and `validate_components()` from
/// `cuda::DynamicsDriver`; the constructor binds the zeroing and
/// kinetic-energy infrastructure components, then user code feeds the
/// rest of the simulation in via the categorical setters (or the raw
/// `add()` / `on(phase).add()` escape hatch).
///
/// @see cuda::minimizer::Minimize for the energy minimizer driver
class MD : public cuda::DynamicsDriver {
  public:
    /// @brief Construct MD driver, binding the phase-agnostic Zero* and
    ///        KineticEnergy infrastructure components to `md::Phase`.
    MD() {
        on(ATV::md::Phase::pre_force).add(ZeroForces{}, -1000);
        on(ATV::md::Phase::pre_force).add(ZeroEnergy{}, -999);
        on(ATV::md::Phase::end_of_step).add(ZeroKinetic{}, -1000);
        on(ATV::md::Phase::end_of_step).add(KineticEnergy{}, -100);
    }

    // ─── Categorical setter API (LAMMPS-shaped) ─────────────────────

    /// @brief Configure the pair-cutoff force evaluator (LAMMPS
    ///        `pair_style`). Bound at `md::Phase::force @ 0`.
    /// @tparam P Pair evaluator type (e.g. `cuda::LJ`).
    /// @param p Pair evaluator (passed by forwarding reference).
    template <typename P>
    MD& set_pair(P&& p) {
        set_singular_at_(std::forward<P>(p), pair_set_, "set_pair",
                         ATV::md::Phase::force, 0);
        return *this;
    }

    /// @brief Configure the bond force evaluator (LAMMPS `bond_style`).
    ///        Bound at `md::Phase::force @ 0`.
    template <typename B>
    MD& set_bond(B&& b) {
        set_singular_at_(std::forward<B>(b), bond_set_, "set_bond",
                         ATV::md::Phase::force, 0);
        return *this;
    }

    /// @brief Configure the angle force evaluator (LAMMPS
    ///        `angle_style`). Bound at `md::Phase::force @ 0`.
    template <typename A>
    MD& set_angle(A&& a) {
        set_singular_at_(std::forward<A>(a), angle_set_, "set_angle",
                         ATV::md::Phase::force, 0);
        return *this;
    }

    /// @brief Configure the dihedral force evaluator (LAMMPS
    ///        `dihedral_style`). Bound at `md::Phase::force @ 0`.
    template <typename D>
    MD& set_dihedral(D&& d) {
        set_singular_at_(std::forward<D>(d), dihedral_set_, "set_dihedral",
                         ATV::md::Phase::force, 0);
        return *this;
    }

    /// @brief Configure the neighbor list (LAMMPS `neighbor`).
    ///        Bound at the dedicated `md::Phase::neighbor @ -10`, which
    ///        runs just before `pre_force` (force/energy zeroing) and the
    ///        `force` phase — mirroring the CPU `cpu::md::MD` neighbor
    ///        placement (`CategoryTraits<md::Phase>::neighbor_phase`).
    template <typename N>
    MD& set_neighbor(N&& n_comp) {
        set_singular_at_(std::forward<N>(n_comp), neighbor_set_, "set_neighbor",
                         ATV::md::Phase::neighbor, -10);
        return *this;
    }

    /// @brief Configure the time integrator (LAMMPS `fix … nve|nvt|npt`).
    ///        Bound via the integrator's own `default_bindings()`.
    template <typename I>
    MD& set_integrator(I&& integ) {
        set_singular_default_(std::forward<I>(integ), integrator_set_,
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
    /// `install_infra_()` auto-injects a `cuda::Virial{}` if a
    /// Barostat-role component is present and no `virial_stress` is
    /// already registered.
    template <typename Bs>
    MD& set_barostat(Bs&& b) {
        this->add(std::forward<Bs>(b));
        return *this;
    }

    /// @brief Add a derived-quantity observer (LAMMPS `compute`). Plural.
    template <typename Obs>
    MD& set_observer(Obs&& obs) {
        this->add(std::forward<Obs>(obs));
        return *this;
    }

    /// @brief Add a recorder / I/O sink (LAMMPS `dump`). Plural.
    template <typename R>
    MD& set_record(R&& r) {
        this->add(std::forward<R>(r));
        return *this;
    }

    /// @brief Configure the thermo print cadence (LAMMPS `thermo N`).
    /// Singular: re-call with a different cadence throws.
    /// Binds `cuda::Thermo{every}` at `end_of_step @ +20`.
    /// @param every Print cadence [steps]. `0` is treated as "every step"
    ///        per the `Thermo` constructor contract.
    MD& set_thermo(int every) {
        if (thermo_every_) {
            if (*thermo_every_ == every)
                return *this;
            throw std::invalid_argument("set_thermo already configured with "
                                        + std::to_string(*thermo_every_)
                                        + ", rejected re-call with "
                                        + std::to_string(every));
        }
        thermo_every_ = every;
        this->on(ATV::md::Phase::end_of_step).add(cuda::Thermo{every}, +20);
        return *this;
    }

    /// @brief Override `CUDADriver::build()` so the auto-injection pass
    ///        runs *before* `validate_components()` sees pending_.
    ///
    /// `install_infra_()` may add a `cuda::Virial{}` if any pending
    /// component declares `(reads() & Capability::Stress) != None`
    /// (covers both `cuda::md::Barostat` and `cuda::md::NoseHooverNPT`).
    /// After the install pass we chain to `cuda::DynamicsDriver::build()`
    /// (which calls `validate_components()` + `Base::build()` in that
    /// order — see `cuda/cuda_driver.cuh`).
    void build() {
        install_infra_();
        cuda::DynamicsDriver::build();
    }

  protected:
    /// @brief Auto-injection pass — single-GPU subset of the CPU rules.
    ///
    /// Two checks, both keyed on `ComponentBase` capability-typed
    /// virtuals (`reads()`, `capabilities()`, `needs_neighbor_list()`,
    /// `provides_neighbor_list()`) — no `name()` string matching,
    /// no `dynamic_cast`:
    ///   1. If any pending component declares `(reads() & Stress) != None`
    ///      AND no pending component already declares
    ///      `(capabilities() & Stress) != None && role()==Observer`,
    ///      add a fresh `cuda::Virial{}` at `post_force @ -50`.
    ///   2. If any pending component reports
    ///      `needs_neighbor_list() == true` without any pending component
    ///      reporting `provides_neighbor_list() == true`, throw naming
    ///      the first offender.
    void install_infra_() {
        bool needs_stress_input = false;
        bool has_stress_producer = false;
        bool has_neighbor_list = false;
        std::string pair_evaluator_name; // first needs_neighbor_list consumer
        this->for_each_pending([&](const cuda::ComponentBase& c) {
            if ((c.reads() & Capability::Stress) != Capability::None)
                needs_stress_input = true;
            if (c.role() == ComponentRole::Observer
                && (c.capabilities() & Capability::Stress) != Capability::None)
                has_stress_producer = true;
            if (c.needs_neighbor_list()) {
                if (pair_evaluator_name.empty())
                    pair_evaluator_name = c.name();
            }
            if (c.provides_neighbor_list())
                has_neighbor_list = true;
        });

        if (needs_stress_input && !has_stress_producer) {
            // `cuda::Virial` has no default_bindings(); bind it
            // explicitly to `post_force @ -50` — the slot the CPU
            // sibling `cpu::md::Virial::default_bindings()` uses,
            // where pair / bonded virials are already accumulated and
            // the velocities are post-integrator-kick so Σ m·v·v
            // reflects the appropriate time-slice.
            this->on(ATV::md::Phase::post_force).add(cuda::Virial{}, -50);
        }

        if (!pair_evaluator_name.empty() && !has_neighbor_list) {
            throw std::runtime_error(
                std::string("cuda::md::MD::build(): component '")
                + pair_evaluator_name
                + "' is a pair evaluator but no neighbor list was registered. "
                  "Add `sim.set_neighbor(cuda::NeighborList{cutoff})`.");
        }
    }

    // ── Singular-setter helpers (CUDA side, mirroring cpu::CPUDriver) ──

    /// @brief Singular setter at an explicit phase / order. Used for
    /// components without `default_bindings()` (NeighborList) or where
    /// the categorical API pins a specific slot.
    template <typename C>
    void set_singular_at_(C&& comp, std::optional<std::string>& slot,
                          const char* setter_name, ATV::md::Phase phase,
                          int order) {
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

    /// @brief Singular setter using the component's `default_bindings()`.
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

  private:
    // ── Categorical setter slots ────────────────────────────────────
    std::optional<std::string> pair_set_;
    std::optional<std::string> bond_set_;
    std::optional<std::string> angle_set_;
    std::optional<std::string> dihedral_set_;
    std::optional<std::string> neighbor_set_;
    std::optional<std::string> integrator_set_;
    std::optional<int> thermo_every_;
};

} // namespace md
} // namespace cuda
} // namespace ATV
