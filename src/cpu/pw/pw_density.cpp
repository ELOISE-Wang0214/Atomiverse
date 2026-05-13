/// @file cpu/pw/pw_density.cpp
#include "atomiverse/cpu/pw/pw_density.h"
#include <cmath>
#include <algorithm>

namespace ATV {
namespace cpu {
namespace pw {

void PW_Density::build(const cpu::BuildContext& ctx) {
    ws_ = static_cast<PWWorkspace*>(ctx.pw_workspace);
    ws_->pulay_depth = pulay_depth_;
    ws_->mix_beta = mix_beta_;
}

void PW_Density::launch(const cpu::LaunchContext&) {
    // Save old density
    ws_->rho_G_old = ws_->rho_G;

    // Build new density from eigenvalues (simplified v1: uniform filling)
    // For now: just apply Pulay mixing to whatever rho_G is set externally

    // Pulay density mixing
    xt::xarray<std::complex<Real>> res = ws_->rho_G - ws_->rho_G_old;

    if (iter_ == 0) {
        // First iteration: simple linear mixing
        for (int i = 0; i < ws_->n_g; ++i)
            ws_->rho_G[i] = ws_->rho_G_old[i] + ws_->mix_beta * res[i];
    } else {
        // Store history
        ws_->pulay_res.push_back(res);
        ws_->pulay_rho.push_back(ws_->rho_G_old);
        if (static_cast<int>(ws_->pulay_res.size()) > ws_->pulay_depth) {
            ws_->pulay_res.erase(ws_->pulay_res.begin());
            ws_->pulay_rho.erase(ws_->pulay_rho.begin());
        }
        // Simple linear mixing as fallback
        for (int i = 0; i < ws_->n_g; ++i)
            ws_->rho_G[i] = ws_->rho_G_old[i] + ws_->mix_beta * res[i];
    }
    ++iter_;
}

} // namespace pw
} // namespace cpu
} // namespace ATV
