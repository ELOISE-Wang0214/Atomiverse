/// @file cpu/minimizer/init.h
/// @brief Minimizer infrastructure components (CPU): ZeroForces, ZeroEnergy.
///
/// These bind to minimizer::Phase::force, NOT md::Phase::pre_force.
/// Each minimizer driver auto-adds these in its constructor.
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/component.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace minimizer {

/// @brief Zero all forces at minimizer force phase.
class ZeroForces : public cpu::ComponentBase {
public:
    std::string name() const override { return "min_zero_forces"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::minimizer::Phase::force), -1000}};
    }
    void launch(const cpu::LaunchContext& ctx) override;
};

/// @brief Zero potential energy at minimizer force phase.
class ZeroEnergy : public cpu::ComponentBase {
public:
    std::string name() const override { return "min_zero_energy"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::minimizer::Phase::force), -999}};
    }
    void launch(const cpu::LaunchContext& ctx) override;
};

} // namespace minimizer
} // namespace cpu
} // namespace ATV
