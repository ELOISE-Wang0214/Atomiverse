/// @file cpu/pw/pw_poisson.cpp
#include "atomiverse/cpu/pw/pw_poisson.h"

namespace ATV {
namespace cpu {
namespace pw {

void PW_Poisson::build(const cpu::BuildContext& ctx) {
    ws_ = static_cast<PWWorkspace*>(ctx.pw_workspace);
}

void PW_Poisson::launch(const cpu::LaunchContext&) {
    for (int i = 0; i < ws_->n_g; ++i) {
        Real g2 = ws_->gvec->g[i].norm_squared();
        if (g2 < Real(1e-16)) {
            ws_->V_H_G[i] = std::complex<Real>(0, 0);
        } else {
            ws_->V_H_G[i] = Real(4.0 * 3.141592653589793) * ws_->rho_G[i] / g2;
        }
    }
}

} // namespace pw
} // namespace cpu
} // namespace ATV
