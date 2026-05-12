/// @file cpu/pw/pw_hamiltonian.cpp
#include "atomiverse/cpu/pw/pw_hamiltonian.h"

namespace ATV {
namespace cpu {
namespace pw {

void PW_Hamiltonian::build(const cpu::BuildContext& ctx) {
    ws_ = static_cast<PWWorkspace*>(ctx.pw_workspace);
    ws_->T_G.resize({static_cast<std::size_t>(ws_->n_g)});
    for (int i = 0; i < ws_->n_g; ++i)
        ws_->T_G[i] = ws_->gvec->ekin[i]; // 0.5*|k+G|^2 already computed
}

void PW_Hamiltonian::launch(const cpu::LaunchContext&) {
    // V_eff_G is already assembled by PW_Veff (runs at +50 after this at 0)
    // This component's launch precomputes T_G — actual H·psi is called from Davidson
}

} // namespace pw
} // namespace cpu
} // namespace ATV
