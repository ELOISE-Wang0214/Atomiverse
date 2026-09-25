/// @file regressions/gpw/gpw-optb88-vdw/reg_gpw-optb88-vdw.cpp
/// @brief Periodic water-dimer GPW optB88-vdW (optB88+PW92 + DRSLL) offline pin.
///
/// Out-of-library consumer. Semi-local XC is `XcKind::OptB88VdW`; nonlocal
/// correlation is post-SCF DRSLL via `vdw_df` (REFERENCE.md). No live oracle
/// in ctest.

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

ATV::cpu::Frame water_dimer_frame() {
    const double ang[6][3] = {
        {0.000000, 0.000000, 0.000000},  // O
        {0.758602, 0.000000, 0.504284},  // H
        {0.758602, 0.000000, -0.504284}, // H
        {3.000000, 0.000000, 0.000000},  // O
        {3.758602, 0.000000, 0.504284},  // H
        {3.758602, 0.000000, -0.504284}, // H
    };
    const int Z[6] = {8, 1, 1, 8, 1, 1};
    ATV::cpu::Frame frame(6);
    for (int i = 0; i < 6; ++i) {
        frame.span(ATV::kAtom_type)[i] = Z[i];
        frame.span(ATV::kAtom_x)[i] = ang[i][0];
        frame.span(ATV::kAtom_y)[i] = ang[i][1];
        frame.span(ATV::kAtom_z)[i] = ang[i][2];
    }
    frame.box = ATV::cpu::Box::cubic(12.0); // [Å]
    return frame;
}

ATV::cpu::gpw::GpwResult run_gpw(ATV::cpu::gpw::XcKind xc, bool vdw_df,
                                 ATV::cpu::Frame& frame) {
    using ATV::cpu::gpw::GaussianBasis;
    using ATV::cpu::gpw::GpwDft;
    using ATV::cpu::gpw::GthPotential;
    using ATV::cpu::gpw::Kind;
    GpwDft gpw{{Kind{molrs::Element::O,
                     GaussianBasis(molrs::Element::O, "DZVP-MOLOPT-SR-GTH"),
                     GthPotential(molrs::Element::O, "GTH-PBE-q6")},
                Kind{molrs::Element::H,
                     GaussianBasis(molrs::Element::H, "DZVP-MOLOPT-SR-GTH"),
                     GthPotential(molrs::Element::H, "GTH-PBE-q1")}}};
    gpw.cutoff_ry = 200.0;
    gpw.n_occ = 8;
    gpw.max_iter = 40;
    gpw.tol = 1e-7;
    gpw.mix = 0.4;
    gpw.xc = xc;
    gpw.nonperiodic = false;
    gpw.vdw_df = vdw_df;
    gpw.build(frame);
    return gpw.run(frame);
}

} // namespace

int main() {
    using namespace ATV;
    using cpu::gpw::XcKind;

    // Offline pins (REFERENCE.md) — first-green ATV freeze; replace with CP2K
    // twin digits when available offline.
    // Placeholders updated after first green run of this binary.
    constexpr double kEDispRef =
        1.5207282287995946e-01; // Ha, first-green freeze
    constexpr double kETotRef =
        -3.4481701207152035e+01;          // Ha, first-green freeze
    constexpr double kTolDisp = 5.0e-5;   // Ha
    constexpr double kTolTot = 1.0e-3;    // Ha (setup integrity until CP2K)
    constexpr double kDualFloor = 1.0e-4; // Ha |E_opt - E_pbe|

    auto frame = water_dimer_frame();

    // Dual gate: semi-local OptB88VdW (no nonlocal) must shift vs PBE.
    const auto res_pbe = run_gpw(XcKind::PBE, /*vdw_df=*/false, frame);
    const auto res_opt = run_gpw(XcKind::OptB88VdW, /*vdw_df=*/false, frame);
    if (!std::isfinite(res_pbe.e_total) || !std::isfinite(res_opt.e_total)) {
        std::fprintf(stderr,
                     "FAIL: non-finite semi-local totals "
                     "PBE=%.14e Opt=%.14e\n",
                     res_pbe.e_total, res_opt.e_total);
        return 1;
    }
    const double d_semi = std::abs(res_opt.e_total - res_pbe.e_total);
    if (d_semi <= kDualFloor) {
        std::fprintf(stderr,
                     "FAIL: OptB88VdW vs PBE dual gate |dE|=%.3e Ha "
                     "(need > %.1e); E_opt=%.14e E_pbe=%.14e\n",
                     d_semi, kDualFloor, res_opt.e_total, res_pbe.e_total);
        return 1;
    }

    // Full optB88-vdW: OptB88VdW + DRSLL.
    const auto res = run_gpw(XcKind::OptB88VdW, /*vdw_df=*/true, frame);
    if (!std::isfinite(res.e_disp) || !std::isfinite(res.e_total)) {
        std::fprintf(stderr,
                     "FAIL: non-finite optB88-vdW E_disp=%.14e "
                     "e_total=%.14e\n",
                     res.e_disp, res.e_total);
        return 1;
    }
    if (!(res.e_disp != 0.0)) {
        std::fprintf(stderr, "FAIL: E_disp is zero (vdw_df inactive?)\n");
        return 1;
    }

    if (std::abs(res.e_disp - kEDispRef) > kTolDisp) {
        std::fprintf(stderr,
                     "FAIL: E_disp optB88-vdW got=%.14e want=%.14e "
                     "tol=%.1e Ha\n",
                     res.e_disp, kEDispRef, kTolDisp);
        return 1;
    }
    if (std::abs(res.e_total - kETotRef) > kTolTot) {
        std::fprintf(stderr,
                     "FAIL: e_total optB88-vdW got=%.14e want=%.14e "
                     "tol=%.1e Ha\n",
                     res.e_total, kETotRef, kTolTot);
        return 1;
    }

    std::printf("PASS: optB88-vdW water dimer E_disp=%.14e Ha "
                "(want %.14e +/- %.1e); e_total=%.14e Ha; "
                "dual |E_opt-E_pbe|=%.6e Ha; SCF iters=%d\n",
                res.e_disp, kEDispRef, kTolDisp, res.e_total, d_semi,
                res.iters);
    return 0;
}
