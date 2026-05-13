/// @file cpu/pw/pw_poisson.h
/// @brief PW_Poisson — Hartree potential solver in reciprocal space.
#pragma once

#include "atomiverse/cpu/component.h"
#include "atomiverse/facade/stage.h"
#include "atomiverse/cpu/pw/pw_workspace.h"

namespace ATV {
namespace cpu {
namespace pw {

class PW_Poisson : public cpu::ComponentBase {
    PWWorkspace* ws_ = nullptr;
public:
    std::string name() const override { return "PW_Poisson"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {DefaultBinding{static_cast<int>(ATV::scf::Phase::build_hamiltonian), -100}};
    }
    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;
};

} // namespace pw
} // namespace cpu
} // namespace ATV
