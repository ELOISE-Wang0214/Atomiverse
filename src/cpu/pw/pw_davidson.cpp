/// @file cpu/pw/pw_davidson.cpp
#include "atomiverse/cpu/pw/pw_davidson.h"
#include <cmath>
#include <algorithm>

namespace ATV {
namespace cpu {
namespace pw {

namespace {
// H·psi for a single wavefunction psi(G) → (H·psi)(G)
void apply_H(const PWWorkspace& ws,
             const std::complex<Real>* psi,
             std::complex<Real>* Hpsi) {
    for (int i = 0; i < ws.n_g; ++i) {
        // Kinetic: T(G) * psi(G)
        Hpsi[i] = ws.T_G[i] * psi[i];
        // Potential convolution: Σ_G' V_eff(G-G') * psi(G')
        // Simplified for v1: diagonal approximation (only G=G' term)
        Hpsi[i] += ws.V_eff_G[i].real() * psi[i];
    }
}
} // anonymous

void PW_Davidson::build(const cpu::BuildContext& ctx) {
    ws_ = static_cast<PWWorkspace*>(ctx.pw_workspace);
}

void PW_Davidson::launch(const cpu::LaunchContext&) {
    // Matrix-free Davidson for the lowest n_bands_ eigenvalues
    // Simplified: just compute band energies from diagonal
    int ng = ws_->n_g;
    ws_->eigenvalues.resize({static_cast<size_t>(n_bands_)});

    // For v1: sorted kinetic+potential energies as eigenvalue estimates
    std::vector<std::pair<Real, int>> diag;
    for (int i = 0; i < ng; ++i) {
        Real e = ws_->T_G[i] + ws_->V_eff_G[i].real();
        diag.push_back({e, i});
    }
    std::sort(diag.begin(), diag.end());
    int nb = std::min(n_bands_, ng);
    for (int b = 0; b < nb; ++b)
        ws_->eigenvalues[b] = diag[b].first;
}

} // namespace pw
} // namespace cpu
} // namespace ATV
