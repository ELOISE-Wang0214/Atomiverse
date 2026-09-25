/// @file regressions/gpw/gpw-revpbe-si/reg_gpw-revpbe-si.cpp
/// @brief END-TO-END: Γ-point GPW revPBE on Si₂ (DZVP-MOLOPT-GTH + GTH-PBE-q4)
///        dual-gate vs PBE + offline pin. See REFERENCE.md.
///
/// Public API only. Requires ATV_BUILD_GPW (+ HF + PW). Exit non-zero on
/// mismatch. ≤ 60 s, single CPU. No live CP2K.

#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/gth_potential.h"
#include "atomiverse/molrs.hpp"

namespace {

// Direct GpwDft path (public headers) — avoids dual Evaluate overhead and
// keeps the e2e case within the ≤ 60 s budget.
ATV::cpu::Frame make_si2_frame() {
    constexpr double a = 5.43; // [Å]
    // Centred in the cell: the Martyna-Tuckerman cluster boundary condition is
    // not translation invariant, and GpwDft::build rejects off-centre
    // nonperiodic systems (a Gaussian clipped at a cell face corrupts the
    // forces while leaving the energy plausible). The fragment's own midpoint
    // is at a/8, so shift by a/2 - a/8 = 3a/8.
    const double s = 3.0 * a / 8.0;
    ATV::cpu::Frame frame(2);
    for (int i = 0; i < 2; ++i) {
        const double p = s + i * a / 4.0;
        frame.span(ATV::kAtom_type)[i] = 14;
        frame.span(ATV::kAtom_x)[i] = p;
        frame.span(ATV::kAtom_y)[i] = p;
        frame.span(ATV::kAtom_z)[i] = p;
    }
    frame.box = ATV::cpu::Box::cubic(a);
    return frame;
}

ATV::cpu::gpw::GpwResult run_xc(ATV::cpu::gpw::XcKind xc,
                                ATV::cpu::Frame& frame) {
    ATV::cpu::gpw::GpwDft m{{ATV::cpu::gpw::Kind{
        molrs::Element::Si,
        ATV::cpu::gpw::GaussianBasis(molrs::Element::Si, "DZVP-MOLOPT-GTH"),
        ATV::cpu::gpw::GthPotential(molrs::Element::Si, "GTH-PBE-q4")}}};
    // 200 Ry: functional dual-gate + e2e budget; absolute pin is soft until
    // offline CP2K twin at the same settings.
    m.cutoff_ry = 200.0;
    m.n_occ = 4; // 2 Si × zion 4 → 8 e⁻
    m.max_iter = 60;
    m.tol = 1e-6;
    m.mix = 0.4;
    m.xc = xc;
    // Isolated MT Poisson: same Si₂ + DZVP-MOLOPT geometry as the CP2K twin's
    // periodic deck (twin stays periodic for offline re-pin). Periodic Γ-only
    // MOLOPT image sums exceed the ≤ 60 s e2e budget on this host; dual-gate
    // still proves revPBE activation on the GPW grid path.
    m.nonperiodic = true;
    m.build(frame);
    return m.run(frame);
}

} // namespace

int main() {
    using namespace ATV::cpu::gpw;
    // Dual-gate floor [Ha].
    constexpr double kDeFloor = 1.0e-4;
    // Soft absolute band from first-green ATV (REFERENCE.md). Replace with
    // offline CP2K twin when available; dual-gate remains the hard criterion.
    // Re-pinned 2026-07-26 on the CENTRED geometry. The previous value
    // (-13.0900394799 Ha) was taken with an atom at the cell origin under
    // Martyna-Tuckerman cluster boundaries, which clips the collocated density:
    // it was ~5.9 Ha (77%) too low. The centred magnitude is consistent with the
    // independent offline CP2K UKS-triplet Si2 reference (-7.613 Ha; see
    // regressions/gpw/gpw-fd-phonon-si/REFERENCE.md). Still a self-golden pin —
    // CP2K cannot converge this closed-shell Si2 in a 5.43 A MT box (verified
    // 2026-07-26, CP2K 2026.1) — so the gate proves functional activation and
    // grid-path wiring, not absolute Si energetics.
    constexpr double kEnergyRef = -7.2202596767; // Ha
    constexpr double kTol = 1.0e-3;              // Ha

    auto frame = make_si2_frame();
    const auto t0 = std::chrono::steady_clock::now();
    const auto r_rev = run_xc(XcKind::RevPBE, frame);
    const auto t1 = std::chrono::steady_clock::now();
    const auto r_pbe = run_xc(XcKind::PBE, frame);
    const auto t2 = std::chrono::steady_clock::now();
    const double e_rev = r_rev.e_total;
    const double e_pbe = r_pbe.e_total;
    const double wall_rev = std::chrono::duration<double>(t1 - t0).count();
    const double wall_pbe = std::chrono::duration<double>(t2 - t1).count();

    if (!r_rev.converged) {
        std::fprintf(stderr,
                     "FAIL: GPW Si revPBE did not converge (iters=%d E=%.8f)\n",
                     r_rev.iters, e_rev);
        return 1;
    }
    if (!std::isfinite(e_rev)) {
        std::fprintf(stderr, "FAIL: non-finite E_revPBE\n");
        return 1;
    }
    if (!r_pbe.converged || !std::isfinite(e_pbe)) {
        std::fprintf(stderr,
                     "FAIL: GPW Si PBE twin did not converge/finite "
                     "(iters=%d E=%.8f)\n",
                     r_pbe.iters, e_pbe);
        return 1;
    }

    const double de = std::abs(e_rev - e_pbe);
    if (!(de > kDeFloor)) {
        std::fprintf(stderr,
                     "FAIL: |E_revPBE−E_PBE|=%.6e Ha  floor=%.1e "
                     "(revPBE inactive?)\n",
                     de, kDeFloor);
        return 1;
    }
    if (std::abs(e_rev - kEnergyRef) > kTol) {
        std::fprintf(stderr,
                     "FAIL: E_revPBE got=%.10f want~%.3f +/- %.1e Ha "
                     "(update pin in REFERENCE after first green)\n",
                     e_rev, kEnergyRef, kTol);
        return 1;
    }

    std::printf("PASS: E_GPW(Si2,revPBE,DZVP-MOLOPT-GTH,MT,200Ry) = %.10f Ha  "
                "E_PBE = %.10f Ha  |dE|=%.3e Ha  iters_rev=%d  "
                "wall_rev=%.1fs wall_pbe=%.1fs  "
                "[dual-gate + offline pin band; see REFERENCE.md]\n",
                e_rev, e_pbe, de, r_rev.iters, wall_rev, wall_pbe);
    return 0;
}
