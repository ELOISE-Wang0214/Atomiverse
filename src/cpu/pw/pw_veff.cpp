/// @file cpu/pw/pw_veff.cpp
#include "atomiverse/cpu/pw/pw_veff.h"

namespace ATV {
namespace cpu {
namespace pw {

void PW_Veff::build(const cpu::BuildContext& ctx) {
    ws_ = static_cast<PWWorkspace*>(ctx.pw_workspace);
}

void PW_Veff::launch(const cpu::LaunchContext&) {
    for (int i = 0; i < ws_->n_g; ++i) {
        ws_->V_eff_G[i] = ws_->V_loc_G_real[i] + ws_->V_H_G[i] + ws_->V_xc_G[i];
    }
}

} // namespace pw
} // namespace cpu
} // namespace ATV
