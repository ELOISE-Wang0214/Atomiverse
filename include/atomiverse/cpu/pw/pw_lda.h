/// @file cpu/pw/pw_lda.h
/// @brief PW_LDA — XC potential on FFT grid via evaluate_xc.
#pragma once

#include "atomiverse/cpu/component.h"
#include "atomiverse/facade/stage.h"
#include "atomiverse/cpu/pw/pw_workspace.h"

namespace ATV {
namespace cpu {
namespace pw {

class PW_LDA : public cpu::ComponentBase {
    PWWorkspace* ws_ = nullptr;
public:
    std::string name() const override { return "PW_LDA"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {DefaultBinding{static_cast<int>(ATV::scf::Phase::build_hamiltonian), -50}};
    }
    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;
};

} // namespace pw
} // namespace cpu
} // namespace ATV
