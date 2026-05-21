/// @file cpu/pw/pw_convergence.cpp
#include "atomiverse/cpu/pw/pw_convergence.h"
#include <cmath>

namespace ATV {
namespace cpu {
namespace pw {

void PW_Convergence::build(const cpu::BuildContext& ctx) {
    ws_ = static_cast<PWWorkspace*>(ctx.pw_workspace);
}

void PW_Convergence::launch(const cpu::LaunchContext&) {
    ++iter_;
    ws_->converged = false;
    if (iter_ <= warmup_)
        return;

    Real diff = 0;
    for (int i = 0; i < ws_->n_g; ++i) {
        Real dr = std::abs(ws_->rho_G[i] - ws_->rho_G_old[i]);
        diff += dr * dr;
    }
    diff = std::sqrt(diff);
    ws_->converged = (diff < tol_);
}

} // namespace pw
} // namespace cpu
} // namespace ATV
