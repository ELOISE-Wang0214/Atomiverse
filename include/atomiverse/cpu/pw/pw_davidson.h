/// @file cpu/pw/pw_davidson.h
/// @brief PW_Davidson — matrix-free iterative diagonalization.
#pragma once

#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/pw/pw_workspace.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace pw {

/// @brief Test hook: one H·ψ application on the bound k-sphere of `ws`
/// (same code path PW_Davidson uses). For validation oracles only.
void pw_apply_h_for_tests(const PWWorkspace& ws,
                          const std::complex<Real>* psi_g,
                          std::complex<Real>* hpsi_g);

class PW_Davidson : public cpu::ComponentBase {
    PWWorkspace* ws_ = nullptr;
    int n_bands_ = 4;
    Real tol_ = Real(1e-8);
    int max_iter_ = 50;

  public:
    PW_Davidson() = default;
    PW_Davidson(int n_bands, Real tol = Real(1e-8), int max_iter = 50)
        : n_bands_(n_bands), tol_(tol), max_iter_(max_iter) {}
    std::string name() const override { return "PW_Davidson"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {
            DefaultBinding{static_cast<int>(ATV::scf::Phase::diagonalize), 0}};
    }
    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;
};

} // namespace pw
} // namespace cpu
} // namespace ATV
