/// @file cpu/pw/pw_convergence.h
/// @brief PW_Convergence — SCF convergence check.
#pragma once

#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/pw/pw_workspace.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace pw {

class PW_Convergence : public cpu::ComponentBase {
    PWWorkspace* ws_ = nullptr;
    Real tol_ = Real(1e-8);
    int warmup_ = 2;
    int iter_ = 0;

  public:
    PW_Convergence() = default;
    PW_Convergence(Real tol, int warmup = 2) : tol_(tol), warmup_(warmup) {}
    std::string name() const override { return "PW_Convergence"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {DefaultBinding{
            static_cast<int>(ATV::scf::Phase::convergence_check), 0}};
    }
    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;
};

} // namespace pw
} // namespace cpu
} // namespace ATV
