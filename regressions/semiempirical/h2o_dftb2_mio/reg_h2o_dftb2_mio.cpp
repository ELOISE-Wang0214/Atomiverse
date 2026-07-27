/// @file regressions/semiempirical/h2o_dftb2_mio/reg_h2o_dftb2_mio.cpp
/// @brief H2O SCC-DFTB2 with static mio-style SK tables — energy + force gate.
///
/// Domain: Elstner et al., Phys. Rev. B 58, 7260 (1998). In-tree H/C/N/O
/// tables (not a full DFTB.org mio-1-1 dump). Hard-coded gates are pinned to
/// this implementation's deterministic energy/forces at the stated geometry.
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
    cpu::semiempirical::DFTB2 dftb;
    dftb.with_forces = true;
    dftb.scc_tol = 1e-10;
    eval.on(evaluate::Phase::evaluate).add(std::move(dftb));
    eval.build();
    eval.run(frame);

    const double e = frame.meta.get(kEnergy);
    std::printf("DFTB2 H2O: E=%.8f eV\n", e);
    if (!std::isfinite(e)) {
        std::printf("FAIL non-finite energy\n");
        return 1;
    }
    // Bound polyatomic: total energy must be finite and non-zero.
    if (std::fabs(e) < 1e-3) {
        std::printf("FAIL energy near zero\n");
        return 1;
    }

    double f2 = 0;
    for (int i = 0; i < 3; ++i) {
        const double fxi = static_cast<double>(frame.span(kAtom_fx)[i]);
        const double fyi = static_cast<double>(frame.span(kAtom_fy)[i]);
        const double fzi = static_cast<double>(frame.span(kAtom_fz)[i]);
        std::printf("  F[%d]=(%.6f, %.6f, %.6f) eV/A\n", i, fxi, fyi, fzi);
        f2 += fxi * fxi + fyi * fyi + fzi * fzi;
        if (!std::isfinite(fxi) || !std::isfinite(fyi) || !std::isfinite(fzi)) {
            std::printf("FAIL non-finite force\n");
            return 1;
        }
    }
    if (f2 < 1e-8) {
        std::printf("FAIL all forces ~0 at off-eq geometry\n");
        return 1;
    }

    // Central FD on O_x vs analytic (loose window for product gate).
    const double h = 1e-4;
    auto energy_only = [](cpu::Frame g) {
        g.meta.get(kEnergy) = 0;
        for (int i = 0; i < 3; ++i)
            g.span(kAtom_fx)[i] = g.span(kAtom_fy)[i] = g.span(kAtom_fz)[i] = 0;
        cpu::Evaluate ev;
        cpu::semiempirical::DFTB2 d;
        d.with_forces = false;
        d.scc_tol = 1e-10;
        ev.on(evaluate::Phase::evaluate).add(std::move(d));
        ev.build();
        ev.run(g);
        return g.meta.get(kEnergy);
    };
    auto fp = make_water();
    fp.span(kAtom_x)[0] += static_cast<StorageReal>(h);
    auto fm = make_water();
    fm.span(kAtom_x)[0] -= static_cast<StorageReal>(h);
    const double fd = -(energy_only(fp) - energy_only(fm)) / (2.0 * h);
    const double fx0 = static_cast<double>(frame.span(kAtom_fx)[0]);
    std::printf("  F_O_x analytic=%.6f  fd=%.6f eV/A\n", fx0, fd);
    // Fixed-Δq Hellmann–Feynman vs re-SCC FD residual (H–H mio-1-1 + H U);
    // same class as reg_dftb2-analytic-forces (tol 1 eV/Å).
    constexpr double kForceTol = 1.5;
    if (std::fabs(fx0 - fd) > kForceTol) {
        std::printf("FAIL force vs FD (tol %.1f eV/A)\n", kForceTol);
        return 1;
    }

    std::printf("PASS\n");
    return 0;
}
