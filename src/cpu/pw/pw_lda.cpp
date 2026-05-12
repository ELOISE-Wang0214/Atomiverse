/// @file cpu/pw/pw_lda.cpp
#include "atomiverse/cpu/pw/pw_lda.h"
#include "atomiverse/cpu/scf/xc_functional.h"
#include <cmath>

namespace ATV {
namespace cpu {
namespace pw {

void PW_LDA::build(const cpu::BuildContext& ctx) {
    ws_ = static_cast<PWWorkspace*>(ctx.pw_workspace);
}

void PW_LDA::launch(const cpu::LaunchContext&) {
    // LDA: V_xc(r) = d(e_xc*rho)/d_rho at each grid point
    // For closed-shell: rho_a = rho_b = rho(r)/2
    for (int i = 0; i < ws_->nr; ++i) {
        Real rho = ws_->rho_r[i];
        if (rho < Real(1e-12)) {
            ws_->V_xc_r[i] = 0;
        } else {
            auto out = cpu::scf::evaluate_xc(
                cpu::scf::XCType::LSDA, rho / Real(2), rho / Real(2));
            ws_->V_xc_r[i] = static_cast<Real>(out.vrho_alpha); // same for alpha and beta
        }
    }
    // FFT not done here — caller or downstream component handles FFT of V_xc_r → V_xc_G
}

} // namespace pw
} // namespace cpu
} // namespace ATV
