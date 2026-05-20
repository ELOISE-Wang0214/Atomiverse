/// @file cpu/md/nose_hoover_npt.h
/// @brief CPU isotropic NPT integrator (single-chain MTK Nose-Hoover).
///        LAMMPS-equivalent `fix npt iso`.
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/component.h"
#include "atomiverse/core/binding.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief Martyna-Tobias-Klein isotropic NPT integrator with M-chain NHC.
///
/// Mirrors LAMMPS `fix npt temp T_target T_target tau_T iso P_target
/// P_target tau_P mtk`. Each of the two NHC objects (one for the atom
/// KE, one for the barostat KE) is a **Nosé-Hoover chain of length M**
/// (default `M = 3`, matching LAMMPS `mtchain`/`mpchain`). The chain
/// provides natural damping of the innermost thermostat-momentum
/// oscillation that a single-chain integrator cannot absorb, which is
/// what keeps an over-compressed start-state (large `P_target − P_int`)
/// from producing the 8000 K-then-freeze ringing the debug_amorphous
/// harness exhibited under length-1 chains.
///
/// Algorithm (Tuckerman 2006, *J. Phys. A* 39, 5629; Shinoda 2004,
/// *Phys. Rev. B* 69, 134103):
///
/// Trotter splitting per timestep `dt`. The outer-NHC pair runs
/// baro-then-atoms on the front half and atoms-then-baro on the back
/// half — required for symmetric Trotter and matches LAMMPS
/// `fix_nh.cpp::initial_integrate` / `final_integrate` ordering:
///
///     [P1] NHC-baro step       (dt/2)   updates η_b, η̇_b; scales ε̇
///     [P2] NHC-atom step       (dt/2)   updates η_t, η̇_t; scales v
///     [P3] barostat-momentum   (dt/2)   ε̇ += G_ε · dt/2
///     [P4] barostat-velocity   (dt/2)   v *= exp(-(1+3/N) · ε̇ · dt/2)
///     [P5] velocity half-kick  (dt/2)   v += (dt/2) · F/m
///     [P6] position+box update (dt)     x = x · e^{ε̇dt} + v · dt · (e^{ε̇dt}-1)/(ε̇dt);
///                                       box *= e^{ε̇dt}
///                                  ↓ (force evaluation; pair virial)
///     [P7] velocity half-kick  (dt/2)   v += (dt/2) · F/m
///     [P8] barostat-velocity   (dt/2)
///     [P9] barostat-momentum   (dt/2)   ε̇ += G_ε(P_new, KE_new) · dt/2
///     [P10] NHC-atom step      (dt/2)
///     [P11] NHC-baro step      (dt/2)
///
/// where the driving force is
///
///     G_ε = [3·V·(P_internal - P_target) + (6/N_dof)·KE] / W
///
/// with `P_internal = (σ_xx + σ_yy + σ_zz) / (3V)` — the LAMMPS-standard
/// pressure-tensor convention (positive trace ⇔ compressive). `P_target`
/// uses the same sign convention as LAMMPS `fix npt iso`, so passing
/// `P_target = 0` matches LAMMPS at the same setpoint. Note: this is
/// the OPPOSITE sign convention from `Barostat` (Berendsen), which still
/// uses `P = -trace/(3V)`; the Berendsen formula is self-consistent
/// under that sign so it does not need to change.
///
/// Phase distribution:
///
///   * `initial_integrate @ 0`  -- steps P1 through P6
///   * `final_integrate   @ 0`  -- step  P7
///   * `end_of_step      @ -40` -- steps P8 through P11
///                                 (after `ReduceEnergy @ -50` so
///                                  pressure & KE are global)
///
/// ── LAMMPS correspondence ───────────────────────────────────────────
///
/// Reference: LAMMPS `src/fix_nh.cpp` (github.com/lammps/lammps,
/// release 22 Jul 2025). Atomiverse implements the same MTK scheme;
/// each helper maps to a `FixNH` method:
///
///   | Atomiverse (this class)        | LAMMPS `FixNH` (`src/fix_nh.cpp`) |
///   |--------------------------------|-----------------------------------|
///   | `build()` (Q/W mass init)      | `FixNH::setup()`                  |
///   | `launch()` initial_integrate   | `FixNH::initial_integrate()`      |
///   | `launch()` final + end_of_step | `FixNH::final_integrate()`        |
///   | `nhc_atoms_half()`             | `FixNH::nhc_temp_integrate()`     |
///   | `nhc_baro_half()`              | `FixNH::nhc_press_integrate()`    |
///   | `update_eps_dot()`             | `FixNH::nh_omega_dot()`           |
///   | `scale_v_by_baro()`            | `FixNH::nh_v_press()`             |
///   | `kick_v_by_force()`            | `FixNH::nve_v()`                  |
///   | `propagate_x_and_box()`        | `FixNH::nve_x()` + `FixNH::remap()`|
///
/// One deliberate divergence: LAMMPS splits the central box operator
/// as `remap(dt/2) → nve_x(dt) → remap(dt/2)`; Atomiverse fuses it
/// into one `propagate_x_and_box()` using the exact analytic drift
/// factor `(e^{ε̇dt}-1)/(ε̇dt)` — equal to LAMMPS's two-remap result
/// to O((ε̇dt)²).
///
/// **Replaces (do not add alongside)** `VelocityVerlet`, `NoseHoover`,
/// and `Barostat`. Compatible with `KineticEnergy` -- this component
/// rewrites `frame.kinetic` after the second NHC update so any thermo
/// output bound at `end_of_step @ > -40` reads the post-NPT value.
class NoseHooverNPT : public cpu::ComponentBase {
public:
    /// @brief All NPT parameters, grouped so each field is named at the
    /// call site and the order of `tau_T` vs `tau_P` can't be swapped
    /// silently. Built as an aggregate (no defaulted fields) so the
    /// compiler flags any omission via `-Wmissing-field-initializers`.
    struct Config {
        double T_target;  ///< Target temperature [K].
        double tau_T;     ///< Thermostat coupling time [fs] (LAMMPS `Tdamp`).
        double P_target;  ///< Target pressure [eV/Å³], LAMMPS-standard
                          ///  sign (`+trace(σ)/(3V)`, positive ⇔
                          ///  compressive) — same convention as LAMMPS
                          ///  `fix npt iso`. 1 bar = 6.2415e-7 eV/Å³.
        double tau_P;     ///< Barostat coupling time [fs] (LAMMPS `Pdamp`).
        double dt;        ///< Timestep [fs].
        double drag = 0.0;  ///< LAMMPS-style empirical drag on the
                            ///  innermost thermostat and barostat NHC
                            ///  momenta (applied at each chain update as
                            ///  `eta_dot[0] *= 1 - drag*half_dt/tau`).
                            ///  Zero disables. `fix npt ... drag 2`.
        int mtchain = 3;    ///< Atom-NHC chain length. LAMMPS default 3.
        int mpchain = 3;    ///< Baro-NHC chain length. LAMMPS default 3.
    };

    /// @param cfg Aggregate NPT configuration. See `Config` for units.
    explicit NoseHooverNPT(Config cfg);

    std::string name() const override { return "nose_hoover_npt"; }
    /// @brief Role tag — Integrator. Although this component subsumes
    /// VelocityVerlet + NoseHoover + Barostat in one Trotter splitting,
    /// only one structural role is exposed; the integrator role is the
    /// load-bearing one for `DynamicsDriver` validation.
    ComponentRole role() const override { return ComponentRole::Integrator; }
    /// @brief Reads the kinetic-stress component of `frame.stress` to
    /// compute the internal pressure feedback. Triggers auto-injection
    /// of `VirialStress` in `MD::install_infra()` when no canonical
    /// stress provider is otherwise registered.
    Capability reads() const override { return Capability::Stress; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::initial_integrate), 0},
                {static_cast<int>(ATV::md::Phase::final_integrate),   0},
                {static_cast<int>(ATV::md::Phase::end_of_step),     -40}};
    }
    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;

private:
    // Configuration
    double T_target_;
    double tau_T_;
    double P_target_;
    double tau_P_;
    double dt_;
    double drag_;
    int    mtchain_;
    int    mpchain_;

    // Derived in build(): (1 - drag * half_dt / tau). `1.0` = no drag.
    double tdrag_factor_ = 1.0;
    double pdrag_factor_ = 1.0;

    // Derived in build()
    int    n_dof_   = 0;     ///< 3N - 3 (COM removed)
    double kT_      = 0.0;   ///< kB · T_target [eV]
    std::vector<double> Q_t_; ///< Atom-NHC masses per chain link [eV·fs²]
    std::vector<double> Q_b_; ///< Baro-NHC masses per chain link [eV·fs²]
    double W_       = 0.0;   ///< Barostat mass [eV·fs²]

    // Extended-Lagrangian state (per-step)
    std::vector<double> eta_t_dot_; ///< Atom-NHC momenta [1/fs]
    std::vector<double> eta_t_;     ///< Atom-NHC positions
    std::vector<double> eta_b_dot_; ///< Baro-NHC momenta [1/fs]
    std::vector<double> eta_b_;     ///< Baro-NHC positions
    double eps_dot_   = 0.0;        ///< Box log-volume momentum [1/fs]

    // Cached global kinetic energy carried between phases
    double KE_global_ = 0.0; ///< 0.5 · Σ m·v² (global, eV)

    // Helpers --------------------------------------------------------
    double compute_global_KE(const cpu::LaunchContext& ctx) const;
    void   nhc_atoms_half(const cpu::LaunchContext& ctx, double half_dt);
    void   nhc_baro_half(double half_dt);
    void   update_eps_dot(const cpu::LaunchContext& ctx, double half_dt);
    void   scale_v_by_baro(const cpu::LaunchContext& ctx, double half_dt);
    void   kick_v_by_force(const cpu::LaunchContext& ctx, double half_dt);
    void   propagate_x_and_box(const cpu::LaunchContext& ctx);
};

} // namespace md
} // namespace cpu
} // namespace ATV
