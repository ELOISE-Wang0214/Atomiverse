/// @file regressions/multiscale/qmmm_link_ch3oh/reg_qmmm_link_ch3oh.cpp
/// @brief Textbook H-link QM/MM on a methanol-like C–O cut (RHF/STO-3G).
///
/// QM = OH + link H; MM = methyl. Gate on last_qm_energy (not LJ total):
/// finite, deep negative, and sensitive to link scale g. See REFERENCE.md.
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/multiscale/qmmm_electrostatic.h"

using namespace ATV;

namespace {

cpu::Frame make_methanol_like() {
    cpu::Frame f(6);
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_x)[1] = -0.96;
    f.span(kAtom_y)[1] = 0.0;
    f.span(kAtom_z)[1] = 0.0;
    f.span(kAtom_x)[2] = 1.43;
    f.span(kAtom_y)[2] = 0.0;
    f.span(kAtom_z)[2] = 0.0;
    f.span(kAtom_x)[3] = 1.8;
    f.span(kAtom_y)[3] = 1.0;
    f.span(kAtom_z)[3] = 0.0;
    f.span(kAtom_x)[4] = 1.8;
    f.span(kAtom_y)[4] = -0.5;
    f.span(kAtom_z)[4] = 0.87;
    f.span(kAtom_x)[5] = 1.8;
    f.span(kAtom_y)[5] = -0.5;
    f.span(kAtom_z)[5] = -0.87;
    const int Z[] = {8, 1, 6, 1, 1, 1};
    for (int i = 0; i < 6; ++i) {
        f.span(kAtom_type)[i] = Z[i];
        f.span(kAtom_tag)[i] = static_cast<TagInt>(i + 1);
        f.span(kAtom_mass)[i] = 1.0;
        f.span(kAtom_fx)[i] = f.span(kAtom_fy)[i] = f.span(kAtom_fz)[i] = 0;
    }
    f.simbox = cpu::Box::cubic(30.0);
    return f;
}

struct Energies {
    double total;
    double qm;
};

Energies energy_link(double scale) {
    auto f = make_methanol_like();
    cpu::QmmmElectrostatic qmmm{{1, 2}, /*mm_q=*/0.0, 0.01, 3.4, 8.0,
                                {{/*qm=*/1, /*mm=*/3}}, scale};
    cpu::Context ctx;
    qmmm.build({ctx, f});
    cpu::LaunchContext lc{ctx, f, 0, 0};
    lc.forces_requested = false;
    f.meta.get(kEnergy) = 0;
    qmmm.launch(lc);
    if (qmmm.last_n_links() != 1)
        throw std::runtime_error("expected one link");
    return {f.meta.get(kEnergy), qmmm.last_qm_energy()};
}

} // namespace

int main() {
    try {
        const auto a = energy_link(0.709);
        const auto b = energy_link(0.55);
        std::printf("qmmm_link_ch3oh E_qm(g=0.709)=%.6f E_qm(g=0.55)=%.6f eV "
                    "dE_qm=%.6f  E_tot=%.6f\n",
                    a.qm, b.qm, a.qm - b.qm, a.total);
        if (!std::isfinite(a.qm) || !std::isfinite(b.qm) ||
            !std::isfinite(a.total)) {
            std::printf("FAIL: non-finite energy\n");
            return 1;
        }
        if (std::fabs(a.qm - b.qm) < 1e-3) {
            std::printf("FAIL: link scale did not shift QM energy\n");
            return 1;
        }
        // STO-3G closed-shell OH+H is deep negative in eV.
        if (!(a.qm < -100.0)) {
            std::printf("FAIL: E_qm=%.6f not below -100 eV gate\n", a.qm);
            return 1;
        }
        std::printf("PASS\n");
        return 0;
    } catch (const std::exception& ex) {
        std::printf("FAIL: %s\n", ex.what());
        return 1;
    }
}
