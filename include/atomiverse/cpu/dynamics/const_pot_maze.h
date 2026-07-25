/// @file cpu/dynamics/const_pot_maze.h
/// @brief Mass-zero-like electrode charge dynamics foothold (MaZe-inspired).
#pragma once

#include <string>
#include <vector>

#include "atomiverse/core/binding.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/dynamics/const_pot_bo.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief MaZe-style path: CG + ASPC on ConstPotBO each step (adiabatic Q).
///
/// Full symplectic mass-zero constrained MD (Coretti et al.) is deferred.
/// This foothold enforces the BO manifold via CG each `pre_force` with ASPC
/// warm-start (`predictor_order=2` by default), exposing charge continuity
/// metrics for regression pins. Composes as a thin owner of `ConstPotBO`.
class ConstPotMaZe : public cpu::Component {
  public:
    explicit ConstPotMaZe(double eta);

    ConstPotMaZe& add_electrode(int type, double psi) &;
    ConstPotMaZe add_electrode(int type, double psi) &&;

    std::string name() const override { return "const_pot_maze"; }

    std::vector<DefaultBinding> default_bindings() const override {
        return {ATV::make_binding(ATV::md::Phase::pre_force, -50)};
    }

    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;
    void destroy(const DestroyContext& ctx) override;

    ConstPotBO& solver() { return bo_; }
    const ConstPotBO& solver() const { return bo_; }

    /// @brief Max |Q_i(t)−Q_i(t−1)| after last launch [e].
    double max_dQ() const { return max_dQ_; }

  private:
    ConstPotBO bo_;
    std::vector<double> Q_prev_;
    double max_dQ_ = 0.0;
    bool have_prev_ = false;
};

} // namespace md
} // namespace cpu
} // namespace ATV
