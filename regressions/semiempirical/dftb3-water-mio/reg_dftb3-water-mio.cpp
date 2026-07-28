/// @file regressions/semiempirical/dftb3-water-mio/reg_dftb3-water-mio.cpp
/// @brief H2O SCC-DFTB3 with static mio-style SK tables — energy gate.
///
/// Domain: Gaus et al., JCTC 7, 931 (2011). Oracle family: DFTB+.
/// In-tree H/C/N/O tables (not full DFTB.org mio-1-1). Public API only.
/// Exit non-zero on failure. <= 60 s.
#include <cmath>
#include <cstdio>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/evaluate/evaluate.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/semiempirical/dftb2.h"
#include "atomiverse/cpu/semiempirical/dftb3.h"

using namespace ATV;

namespace {

// Provisional ATV Hamiltonian-twin pin [eV] at the geometry below
// (scc_tol=1e-10, max_scc=50). Measured 2026-07-28 Release. NOT a DFTB.org
// mio-1-1 absolute energy and NOT a self-consistent energy — see REFERENCE.md:
// the SCC Picard loop is a period-2 limit cycle here, so this value is the
// iteration-`max_scc` branch of that cycle. DFTB+ oracle owed.
// Supersedes -84.561719022489, which was recorded under a third-order SCC
// potential that was not the derivative of the reported E^(3).
constexpr double kPinEnergyEV = -107.618793920072;
constexpr double kPinTolEV = 1e-6;

cpu::Frame make_water() {
    cpu::Frame f(3);
    f.span(kAtom_type)[0] = 8;
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_type)[1] = 1;
    f.span(kAtom_x)[1] = 0.96;
    f.span(kAtom_y)[1] = 0.0;
    f.span(kAtom_z)[1] = 0.0;
    f.span(kAtom_type)[2] = 1;
    f.span(kAtom_x)[2] = -0.24;
    f.span(kAtom_y)[2] = 0.93;
    f.span(kAtom_z)[2] = 0.0;
    for (int i = 0; i < 3; ++i) {
        f.span(kAtom_mass)[i] = (i == 0) ? 15.999 : 1.008;
        f.span(kAtom_tag)[i] = i + 1;
        f.span(kAtom_fx)[i] = f.span(kAtom_fy)[i] = f.span(kAtom_fz)[i] = 0;
    }
    f.box = cpu::Box::cubic(20.0);
    f.meta.get(kEnergy) = 0;
    return f;
}

} // namespace

int main() {
    auto frame = make_water();

    cpu::Evaluate eval;
    cpu::semiempirical::DFTB3 dftb3;
    dftb3.with_forces = false;
    dftb3.scc_tol = 1e-10;
    eval.on(evaluate::Phase::evaluate).add(std::move(dftb3));
    eval.build();
    eval.run(frame);

    const double e = frame.meta.get(kEnergy);
    std::printf("DFTB3 H2O: E=%.12f eV\n", e);
    if (!std::isfinite(e)) {
        std::printf("FAIL non-finite energy\n");
        return 1;
    }
    if (std::fabs(e) < 1e-3) {
        std::printf("FAIL energy near zero\n");
        return 1;
    }

    // DFTB2 on same geometry — third-order must move the energy.
    auto frame2 = make_water();
    cpu::Evaluate eval2;
    cpu::semiempirical::DFTB2 d2;
    d2.with_forces = false;
    d2.scc_tol = 1e-10;
    eval2.on(evaluate::Phase::evaluate).add(std::move(d2));
    eval2.build();
    eval2.run(frame2);
    const double e2 = frame2.meta.get(kEnergy);
    std::printf("DFTB2 H2O: E=%.12f eV  |dE|=%.6e\n", e2, std::fabs(e - e2));
    if (std::fabs(e - e2) < 1e-6) {
        std::printf("FAIL DFTB3 energy identical to DFTB2 (third-order inactive)\n");
        return 1;
    }

    if (std::fabs(e - kPinEnergyEV) > kPinTolEV) {
        std::printf("FAIL energy pin: got %.12f expected %.12f (tol %g)\n", e,
                    kPinEnergyEV, kPinTolEV);
        return 1;
    }
    std::printf("PIN match: %.12f eV (tol %g)\n", kPinEnergyEV, kPinTolEV);

    std::printf("PASS\n");
    return 0;
}
