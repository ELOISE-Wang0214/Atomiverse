/// @file cpu/md/init.h
/// @brief CPU MD infrastructure components: ZeroForces, ZeroEnergy, ZeroKinetic,
///        KineticEnergy, InitVelocity.
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/component.h"
#include "atomiverse/core/binding.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief Zero all forces. Phase-agnostic — driver binds to its own phase.
class ZeroForces : public cpu::ComponentBase {
public:
    std::string name() const override { return "zero_forces"; }
    ComponentRole role() const override { return ComponentRole::None; }
    void launch(const cpu::LaunchContext& ctx) override;
};

/// @brief Zero potential energy. Phase-agnostic — driver binds to its own phase.
class ZeroEnergy : public cpu::ComponentBase {
public:
    std::string name() const override { return "zero_energy"; }
    ComponentRole role() const override { return ComponentRole::None; }
    void launch(const cpu::LaunchContext& ctx) override;
};

/// @brief Zero stress tensor. Phase-agnostic — driver binds to its own phase.
///
/// Stress accumulation by force kernels (LJ pair virial, future bonded
/// virials) requires a clean baseline at the start of each step.
class ZeroStress : public cpu::ComponentBase {
public:
    std::string name() const override { return "zero_stress"; }
    ComponentRole role() const override { return ComponentRole::None; }
    void launch(const cpu::LaunchContext& ctx) override;
};

/// @brief Zero kinetic energy. Phase-agnostic — driver binds to its own phase.
class ZeroKinetic : public cpu::ComponentBase {
public:
    std::string name() const override { return "zero_kinetic"; }
    ComponentRole role() const override { return ComponentRole::None; }
    void launch(const cpu::LaunchContext& ctx) override;
};

/// @brief Compute kinetic energy at end_of_step phase.
///
/// KE = 0.5 * kMvv2e * sum(m * v^2), with OpenMP reduction.
class KineticEnergy : public cpu::ComponentBase {
public:
    std::string name() const override { return "kinetic_energy"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::end_of_step), -100}};
    }
    void launch(const cpu::LaunchContext& ctx) override;
};

/// @brief Remove net linear COM momentum every `freq` steps.
///
/// LAMMPS `fix momentum N linear 1 1 1` equivalent. Prevents flying-ice-cube
/// drift that accumulates under Langevin + NPT coupling (the stochastic
/// Langevin kicks plus barostat scaling slowly bias the COM away from 0,
/// and fix_npt alone doesn't zero it).  Bound late in `end_of_step` so
/// the freshly-integrated velocities are zeroed before the next step.
class ZeroLinearMomentum : public cpu::ComponentBase {
public:
    /// @param freq Remove momentum every `freq` steps (LAMMPS default 1).
    explicit ZeroLinearMomentum(int freq = 1) : freq_(freq) {}

    std::string name() const override { return "zero_linear_momentum"; }
    std::vector<DefaultBinding> default_bindings() const override {
        // +25 = after KineticEnergy (@-100) / ReduceEnergy (@-50) /
        // NoseHooverNPT (@-40) so we only touch post-integration
        // velocities.  ThermoPrinter sits at +20, also behind us.
        return {{static_cast<int>(ATV::md::Phase::end_of_step), +25}};
    }
    void launch(const cpu::LaunchContext& ctx) override;

private:
    int freq_;
    long long step_ = 0;
};

/// @brief Maxwell-Boltzmann velocity initialization (build-time only).
///
/// Generates velocities from Gaussian distribution scaled by sqrt(kB*T/m),
/// then removes COM velocity. No-op during launch().
class InitVelocity : public cpu::ComponentBase {
public:
    /// @param temperature Target temperature [K].
    /// @param seed Random number generator seed.
    InitVelocity(double temperature, unsigned seed);

    std::string name() const override { return "init_velocity"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::initial_integrate), -999}};
    }
    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override {}

private:
    double temperature_;
    unsigned seed_;
};

} // namespace md
} // namespace cpu
} // namespace ATV
