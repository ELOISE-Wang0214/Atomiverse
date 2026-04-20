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

/// @brief Single-chain Martyna-Tobias-Klein isotropic NPT integrator.
///
/// Mirrors LAMMPS `fix npt temp T_target T_target tau_T iso P_target
/// P_target tau_P` (the isotropic-cell variant of `fix_nh.cpp`) but with
/// **single thermostat chains** rather than the LAMMPS default of length 3.
/// For systems large enough that the chain length is statistically
/// irrelevant (which is the regime where one would compare the two
/// codes), the ensemble averages match those of LAMMPS.
///
/// Algorithm (Tuckerman 2006, *J. Phys. A* 39, 5629; Shinoda 2004,
/// *Phys. Rev. B* 69, 134103):
///
/// Trotter splitting per timestep `dt`:
///
///     [P1] NHC-atom step       (dt/2)   updates η_t, η̇_t; scales v
///     [P2] NHC-baro step       (dt/2)   updates η_b, η̇_b; scales ε̇
///     [P3] barostat-momentum   (dt/2)   ε̇ += G_ε · dt/2
///     [P4] barostat-velocity   (dt/2)   v *= exp(-(1+3/N) · ε̇ · dt/2)
///     [P5] velocity half-kick  (dt/2)   v += (dt/2) · F/m
///     [P6] position+box update (dt)     x = x · e^{ε̇dt} + v · dt · sinhc(ε̇dt/2);
///                                       box *= e^{ε̇dt}
///                                  ↓ (force evaluation; pair virial)
///     [P7] velocity half-kick  (dt/2)   v += (dt/2) · F/m
///     [P8] barostat-velocity   (dt/2)
///     [P9] barostat-momentum   (dt/2)   ε̇ += G_ε(P_new, KE_new) · dt/2
///     [P10] NHC-baro step      (dt/2)
///     [P11] NHC-atom step      (dt/2)
///
/// where the driving force is
///
///     G_ε = [3·V·(P_target - P_internal) + (6/N_dof)·KE] / W
///
/// `P_internal = -(σ_xx + σ_yy + σ_zz) / (3V)` matches the convention
/// used by `Barostat` (Berendsen), so `P_target = 0` means "atmospheric"
/// in both components.
///
/// Phase distribution:
///
///   * `initial_integrate @ 0`  -- steps P1 through P6
///   * `final_integrate   @ 0`  -- step  P7
///   * `end_of_step      @ -40` -- steps P8 through P11
///                                 (after `ReduceEnergy @ -50` so
///                                  pressure & KE are global)
///
/// **Replaces (do not add alongside)** `VelocityVerlet`, `NoseHoover`,
/// and `Barostat`. Compatible with `KineticEnergy` -- this component
/// rewrites `frame.kinetic` after the second NHC update so any thermo
/// output bound at `end_of_step @ > -40` reads the post-NPT value.
class NoseHooverNPT : public cpu::ComponentBase {
public:
    /// @param T_target Target temperature [K].
    /// @param tau_T    Thermostat coupling time [fs] (LAMMPS `Tdamp`).
    /// @param P_target Target pressure in `frame.stress` convention
    ///                 (`-trace/(3V)`); pass 0 for "atmospheric".
    /// @param tau_P    Barostat coupling time [fs] (LAMMPS `Pdamp`).
    /// @param dt       Timestep [fs].
    NoseHooverNPT(double T_target, double tau_T,
                  double P_target, double tau_P,
                  double dt);

    std::string name() const override { return "nose_hoover_npt"; }
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

    // Derived in build()
    int    n_dof_   = 0;     ///< 3N - 3 (COM removed)
    double kT_      = 0.0;   ///< kB · T_target [eV]
    double Q_t_     = 0.0;   ///< Atom-thermostat mass [eV·fs²]
    double Q_b_     = 0.0;   ///< Baro-thermostat mass [eV·fs²]
    double W_       = 0.0;   ///< Barostat mass [eV·fs²]

    // Extended-Lagrangian state (per-step)
    double eta_t_dot_ = 0.0; ///< Atom-thermostat momentum [1/fs]
    double eta_t_     = 0.0; ///< Atom-thermostat position
    double eta_b_dot_ = 0.0; ///< Baro-thermostat momentum [1/fs]
    double eta_b_     = 0.0; ///< Baro-thermostat position
    double eps_dot_   = 0.0; ///< Box log-volume momentum [1/fs]

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
