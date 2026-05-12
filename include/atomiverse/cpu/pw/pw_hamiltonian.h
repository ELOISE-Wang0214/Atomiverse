/// @file cpu/pw/pw_hamiltonian.h
/// @brief PW_Hamiltonian — kinetic energy + effective potential operator.
#pragma once

#include "atomiverse/cpu/component.h"
#include "atomiverse/facade/stage.h"
#include "atomiverse/cpu/pw/pw_workspace.h"

namespace ATV {
namespace cpu {
namespace pw {

class PW_Hamiltonian : public cpu::ComponentBase {
    PWWorkspace* ws_ = nullptr;
public:
    std::string name() const override { return "PW_Hamiltonian"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {DefaultBinding{static_cast<int>(ATV::scf::Phase::build_hamiltonian), 0}};
    }
    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;
};

} // namespace pw
} // namespace cpu
} // namespace ATV
