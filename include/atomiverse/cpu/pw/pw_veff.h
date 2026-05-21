/// @file cpu/pw/pw_veff.h
/// @brief PW_Veff — assemble effective potential V_eff = V_loc + V_H + V_xc.
#pragma once

#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/pw/pw_workspace.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace pw {

class PW_Veff : public cpu::ComponentBase {
    PWWorkspace* ws_ = nullptr;

  public:
    std::string name() const override { return "PW_Veff"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {DefaultBinding{
            static_cast<int>(ATV::scf::Phase::build_hamiltonian), 50}};
    }
    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;
};

} // namespace pw
} // namespace cpu
} // namespace ATV
