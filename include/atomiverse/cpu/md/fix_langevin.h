/// @file cpu/md/fix_langevin.h
/// @brief CPU Langevin thermostat -- LAMMPS-equivalent `fix langevin`.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "atomiverse/cpu/component.h"
#include "atomiverse/core/binding.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief Schneider-Stoll Langevin thermostat applied as an additive force.
///
/// Mirrors LAMMPS `fix langevin Tstart Tstop damp seed`
/// (`lammps/src/fix_langevin.cpp`) in the constant-temperature limit
/// (`Tstart == Tstop == T_target`). At every step, before the velocity
/// update, two forces are added per atom:
///
/// 1. Friction (Stokes drag):  F_drag = - m / (damp · ftm2v) · v
/// 2. Random kick:              F_rand =   sqrt(2 m kB T / (damp · dt · mvv2e)) / ftm2v · ξ
///
/// where `ξ ~ N(0, 1)` per Cartesian component per atom per step,
/// `ftm2v = kFtm2v` and `mvv2e = kMvv2e` (the metal-unit conversion
/// constants). The combined force satisfies the fluctuation-dissipation
/// theorem so velocity-Verlet integration of `f_total = f_phys + F_drag
/// + F_rand` produces the canonical NVT ensemble at `T_target`.
///
/// LAMMPS uses uniform random samples drawn from the Marsaglia-Tsang
/// `RanMars` generator. We use Gaussian samples drawn from a per-(seed,
/// step, atom-tag) splitmix stream. The two are statistically
/// equivalent in the limit of many particles -- ⟨T⟩ and the velocity
/// distribution match -- but individual trajectories will differ
/// (different RNGs).
///
/// Bound at `post_force @ +500` so the thermostat forces are added
/// after every physical-force evaluator has accumulated, but before
/// any user-defined post-force corrections that might depend on the
/// total. Pair this component with `VelocityVerlet` (or `NVELimit`
/// for far-from-equilibrium starts) -- the integrator picks up the
/// modified forces in the next velocity half-kick.
class FixLangevin : public cpu::ComponentBase {
public:
    /// @param T_target Target temperature [K].
    /// @param damp     Damping time constant [fs] (LAMMPS `damp`).
    ///                 Effective friction γ = 1 / damp.
    /// @param dt       Timestep [fs].
    /// @param seed     RNG seed; reproducible across rank counts.
    FixLangevin(double T_target, double damp, double dt, std::uint32_t seed);

    std::string name() const override { return "fix_langevin"; }
    ComponentRole role() const override { return ComponentRole::Thermostat; }
    /// @brief Langevin's role is Thermostat, but it functionally writes
    /// drag + random forces into fx/fy/fz — so it advertises
    /// Capability::Force. This lets a pure Langevin run (no LJ / bonded
    /// evaluator) satisfy `DynamicsDriver`'s "something writes forces"
    /// contract.
    Capability capabilities() const override { return Capability::Force; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::post_force), +500}};
    }
    void launch(const cpu::LaunchContext& ctx) override;

private:
    double T_target_;
    double damp_;
    double dt_;
    std::uint32_t seed_;
};

} // namespace md
} // namespace cpu
} // namespace ATV
