/// @file regressions/semiempirical/spin-dftb2-o2/reg_spin-dftb2-o2.cpp
/// @brief O₂ collinear spin-polarized SCC-DFTB2 — energy + Mulliken spin gate.
///
/// Method: collinear spin-DFTB2 (Elstner SCC + Köhler / DFTB+ spin channel).
/// In-tree mio-style H/C/N/O tables (not DFTB.org mio-1-1 bit-identical).
/// Absolute energy pin is the ATV Hamiltonian twin (deterministic).
/// Family oracle: DFTB+ (offline when SK fidelity matches); literature method.
///
/// Exit non-zero on failure. Public API only. <= 60 s.
#include <cmath>
#include <cstdio>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/evaluate/evaluate.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/semiempirical/dftb2.h"

using namespace ATV;

namespace {

cpu::Frame make_o2(double r_ang = 1.21) {
    cpu::Frame f(2);
    f.span(kAtom_type)[0] = f.span(kAtom_type)[1] = 8;
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_x)[1] = r_ang;
    f.span(kAtom_y)[1] = 0.0;
    f.span(kAtom_z)[1] = 0.0;
    for (int i = 0; i < 2; ++i) {
        f.span(kAtom_mass)[i] = 15.999;
        f.span(kAtom_tag)[i] = i + 1;
        f.span(kAtom_fx)[i] = f.span(kAtom_fy)[i] = f.span(kAtom_fz)[i] = 0;
    }
    f.box = cpu::Box::cubic(30.0);
    f.meta.get(kEnergy) = 0;
    return f;
}

} // namespace

int main() {
    auto frame = make_o2();

    cpu::Evaluate eval;
    cpu::semiempirical::DFTB2 dftb;
    dftb.charge = 0;
    dftb.multiplicity = 3; // ground-state triplet O₂
    dftb.with_forces = false;
    dftb.scc_tol = 1e-10;
    dftb.max_scc = 80;
    // Keep a reference to the component for Mulliken spin after run: register
    // via add_ref so ownership stays here.
    eval.on(evaluate::Phase::evaluate).add_ref(dftb);
    eval.build();
    eval.run(frame);

    const double e = frame.meta.get(kEnergy);
    std::printf("spin-DFTB2 O2 triplet: E=%.10f eV\n", e);
    if (!std::isfinite(e)) {
        std::printf("FAIL non-finite energy\n");
        return 1;
    }
    if (std::fabs(e) < 1e-3) {
        std::printf("FAIL energy near zero\n");
        return 1;
    }

    // Hard pin: ATV mio-style Hamiltonian twin (see REFERENCE.md).
    // Offline DFTB+ absolute identity needs SK fidelity (not this slug).
    constexpr double kE_pin_eV = -142.1382039484;
    const double de = std::fabs(e - kE_pin_eV);
    std::printf("  |E - pin|=%.6e eV (tol 5e-2)\n", de);
    if (de > 5e-2) {
        std::printf("FAIL energy pin\n");
        return 1;
    }

    const auto& m = dftb.last_mulliken_spin();
    if (m.size() != 2) {
        std::printf("FAIL spin population size %zu\n", m.size());
        return 1;
    }
    const double m0 = m[0];
    const double m1 = m[1];
    const double m_tot = m0 + m1;
    std::printf("  m_Mulliken = [%.6f, %.6f]  sum=%.6f\n", m0, m1, m_tot);
    if (std::fabs(m_tot - 2.0) > 0.25) {
        std::printf("FAIL total spin population (expect ~2 for triplet)\n");
        return 1;
    }
    if (std::fabs(m0) < 0.3 || std::fabs(m1) < 0.3) {
        std::printf("FAIL per-atom spin too small\n");
        return 1;
    }
    if (std::fabs(m0 - m1) > 0.2) {
        std::printf("FAIL spin asymmetry on O2\n");
        return 1;
    }

    const auto& dq = dftb.last_mulliken_dq();
    if (dq.size() == 2) {
        std::printf("  dq_Mulliken = [%.6f, %.6f]\n", dq[0], dq[1]);
        if (!std::isfinite(dq[0]) || !std::isfinite(dq[1])) {
            std::printf("FAIL non-finite charge\n");
            return 1;
        }
    }

    std::printf("PASS\n");
    return 0;
}
