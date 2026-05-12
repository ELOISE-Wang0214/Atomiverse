/// @file cpu/pw/pw_density.h
/// @brief PW_Density — build charge density from wavefunctions + Pulay mixing.
#pragma once

#include "atomiverse/cpu/component.h"
#include "atomiverse/facade/stage.h"
#include "atomiverse/cpu/pw/pw_workspace.h"

namespace ATV {
namespace cpu {
namespace pw {

class PW_Density : public cpu::ComponentBase {
    PWWorkspace* ws_ = nullptr;
    Real mix_beta_ = Real(0.4);
    int pulay_depth_ = 5;
    int iter_ = 0;
public:
    PW_Density() = default;
    PW_Density(Real beta, int pulay_depth = 5)
        : mix_beta_(beta), pulay_depth_(pulay_depth) {}
    std::string name() const override { return "PW_Density"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {DefaultBinding{static_cast<int>(ATV::scf::Phase::build_density), 0}};
    }
    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;
};

} // namespace pw
} // namespace cpu
} // namespace ATV
