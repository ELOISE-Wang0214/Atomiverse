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

/// @brief Zero all forces at pre_force phase.
class ZeroForces : public cpu::ComponentBase {
public:
    std::string name() const override { return "zero_forces"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::pre_force), -1000}};
    }
    void launch(const cpu::LaunchContext& ctx) override;
};

/// @brief Zero potential energy at pre_force phase.
class ZeroEnergy : public cpu::ComponentBase {
public:
    std::string name() const override { return "zero_energy"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::pre_force), -999}};
    }
    void launch(const cpu::LaunchContext& ctx) override;
};

/// @brief Zero stress tensor at pre_force phase.
///
/// Stress accumulation by force kernels (LJ pair virial, future bonded
/// virials) requires a clean baseline at the start of each step.
/// Bound at `pre_force @ -998`, after ZeroForces / ZeroEnergy.
class ZeroStress : public cpu::ComponentBase {
public:
    std::string name() const override { return "zero_stress"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::pre_force), -998}};
    }
    void launch(const cpu::LaunchContext& ctx) override;
};

/// @brief Zero kinetic energy at end_of_step phase.
class ZeroKinetic : public cpu::ComponentBase {
public:
    std::string name() const override { return "zero_kinetic"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::end_of_step), -1000}};
    }
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
