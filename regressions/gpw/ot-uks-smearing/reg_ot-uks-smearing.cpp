/// @file regressions/gpw/ot-uks-smearing/reg_ot-uks-smearing.cpp
/// @brief END-TO-END: multi-k GPW UKS OT + Fermi–Dirac metallic smearing
///        (BCC Na).
///
/// See REFERENCE.md. Public headers only. Requires ATV_BUILD_GPW.
/// Gates: OT+FD vs DIIS+FD free-energy identity, electron count, frozen pin.
/// ctest never launches CP2K.

#include <chrono>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/gth_potential.h"
#include "atomiverse/molrs.hpp"

namespace {

int g_fail = 0;

void check(bool ok, const char* what) {
    if (!ok) {
        std::fprintf(stderr, "FAIL: %s\n", what);
        ++g_fail;
    } else {
        std::printf("  ok: %s\n", what);
    }
}

/// BCC conventional cubic: 2 Na at (0,0,0) and (½,½,½).
ATV::cpu::Frame na_bcc(double a_ang) {
    using namespace ATV;
    const double frac[2][3] = {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}};
    cpu::Frame frame(2);
    for (int i = 0; i < 2; ++i) {
        frame.span(kAtom_type)[i] = 11;
        frame.span(kAtom_x)[i] = frac[i][0] * a_ang;
        frame.span(kAtom_y)[i] = frac[i][1] * a_ang;
        frame.span(kAtom_z)[i] = frac[i][2] * a_ang;
    }
    frame.box = cpu::Box::cubic(a_ang);
    return frame;
}

ATV::cpu::gpw::GpwDft base_method() {
    using namespace ATV;
    cpu::gpw::GpwDft m{{cpu::gpw::Kind{
        molrs::Element::Na,
        cpu::gpw::GaussianBasis(molrs::Element::Na, "DZVP-MOLOPT-SR-GTH"),
        cpu::gpw::GthPotential(molrs::Element::Na, "GTH-PBE-q9")}}};
    m.cutoff_ry = 150.0;
    m.n_occ = 0;
    m.n_alpha = 10; // ferro seed; free under smear
    m.n_beta = 8;
    m.max_iter = 120;
    m.tol = 1e-6;
    m.eps_default = 1e-7;
    m.xc = ATV::cpu::gpw::XcKind::PADE;
    m.diis_depth = 6;
    m.mix = 0.4;
    m.k_mesh[0] = m.k_mesh[1] = m.k_mesh[2] = 2;
    m.k_time_reversal = true;
    m.smearing = ATV::cpu::gpw::Smearing::fermi_dirac;
    m.degauss = 0.02; // Ry
    m.extra_bands = 4;
    return m;
}

} // namespace

int main() {
    using namespace ATV;
    using cpu::gpw::GpwDft;
    using cpu::gpw::Smearing;

    const double a_ang = 4.225;
    const double n_elec = 18.0;
    cpu::Frame frame = na_bcc(a_ang);

    // ── A1: build accepts UKS OT + FD ─────────────────────────────────────
    {
        GpwDft m = base_method();
        m.use_ot = true;
        try {
            m.build(frame);
            check(true, "A1 build accepts UKS OT + FD smear");
        } catch (const std::exception& e) {
            std::fprintf(stderr, "A1 threw: %s\n", e.what());
            check(false, "A1 build accepts UKS OT + FD smear");
        }
    }

    // ── A2/A3: OT+FD vs DIIS+FD free-energy identity ───────────────────────
    double f_ot = 0.0;
    double f_diis = 0.0;
    double wall_ot = 0.0;
    {
        GpwDft diis = base_method();
        diis.use_ot = false;
        diis.build(frame);
        const auto& rd = diis.run(frame);
        check(rd.converged, "A2 DIIS+FD converged");
        f_diis = rd.e_total;
        check(std::fabs(rd.n_elec_occ - n_elec) < 1e-5, "A2 DIIS Ne");

        GpwDft ot = base_method();
        ot.use_ot = true;
        ot.max_iter = 200;
        ot.build(frame);
        const auto t2 = std::chrono::steady_clock::now();
        const auto& ro = ot.run(frame);
        const auto t3 = std::chrono::steady_clock::now();
        wall_ot = std::chrono::duration<double>(t3 - t2).count();
        check(ro.converged, "A2 OT+FD converged");
        f_ot = ro.e_total;
        check(std::fabs(ro.n_elec_occ - n_elec) < 1e-5, "A4 OT Ne");
        check(std::isfinite(ro.e_fermi) && std::fabs(ro.e_fermi) < 1e3,
              "A6 finite e_fermi");
        check(std::isfinite(ro.e_entropy), "A6 finite e_entropy");
        check(ro.n_kpts > 1, "A2 multi-k active");

        const double d = std::fabs(f_ot - f_diis);
        std::printf("A3: F_OT=%.10f  F_DIIS=%.10f  |Δ|=%.3e Ha  "
                    "iters ot/diis=%d/%d  wall_ot=%.3f s\n",
                    f_ot, f_diis, d, ro.iters, rd.iters, wall_ot);
        check(d <= 1e-4, "A3 |F_OT - F_DIIS| <= 1e-4 Ha");
    }

    // ── A5: frozen first-green free-energy band ────────────────────────────
    // Written after first green on this branch (REFERENCE.md).
    {
        constexpr double kFrozenFHa = -92.1263304935; // OT first green
        constexpr double kAbsTol = 1e-3;
        const double d = std::fabs(f_ot - kFrozenFHa);
        std::printf("A5: F_OT=%.10f  frozen=%.10f  |Δ|=%.3e\n", f_ot,
                    kFrozenFHa, d);
        check(d <= kAbsTol, "A5 |F_OT - frozen| <= 1e-3 Ha");
    }

    if (g_fail) {
        std::printf("FAIL: ot-uks-smearing (%d checks)\n", g_fail);
        return 1;
    }
    std::printf("PASS: ot-uks-smearing (BCC Na UKS OT+FD ≡ DIIS+FD)\n");
    return 0;
}
