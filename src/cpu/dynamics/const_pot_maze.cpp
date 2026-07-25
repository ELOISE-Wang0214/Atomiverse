/// @file cpu/dynamics/const_pot_maze.cpp
#include "atomiverse/cpu/dynamics/const_pot_maze.h"

#include <cmath>

namespace ATV {
namespace cpu {
namespace md {

ConstPotMaZe::ConstPotMaZe(double eta) : bo_(eta) {
    bo_.use_cg = true;
    bo_.predictor_order = 2;
    bo_.electroneutrality = true;
}

ConstPotMaZe& ConstPotMaZe::add_electrode(int type, double psi) & {
    bo_.add_electrode(type, psi);
    return *this;
}
ConstPotMaZe ConstPotMaZe::add_electrode(int type, double psi) && {
    bo_.add_electrode(type, psi);
    return std::move(*this);
}

void ConstPotMaZe::build(const BuildContext& ctx) {
    bo_.build(ctx);
    Q_prev_.clear();
    have_prev_ = false;
    max_dQ_ = 0.0;
}

void ConstPotMaZe::launch(const LaunchContext& ctx) {
    bo_.launch(ctx);
    const auto& Q = bo_.charges();
    max_dQ_ = 0.0;
    if (have_prev_ && Q.size() == Q_prev_.size()) {
        for (std::size_t i = 0; i < Q.size(); ++i)
            max_dQ_ = std::max(max_dQ_, std::abs(Q[i] - Q_prev_[i]));
    }
    Q_prev_ = Q;
    have_prev_ = true;
}

void ConstPotMaZe::destroy(const DestroyContext& ctx) {
    bo_.destroy(ctx);
    Q_prev_.clear();
    have_prev_ = false;
}

} // namespace md
} // namespace cpu
} // namespace ATV
