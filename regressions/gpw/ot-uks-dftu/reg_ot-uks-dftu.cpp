/// @file regressions/gpw/ot-uks-dftu/reg_ot-uks-dftu.cpp
/// @brief END-TO-END: multi-k collinear UKS GPW OT + Dudarev DFT+U.
///
/// Public `GpwDft` only. Offline CP2K twin under
/// `tools/cp2k_oracle/inputs/ot-uks-dftu/` (REFERENCE.md). ctest never
/// spawns CP2K.
///
/// Gates:
///   A1 Γ-limit: explicit k=Γ OT+U ≡ multi-k DIIS+U at Γ (complex UKS path)
///   A2 multi-k OT+U ≡ multi-k DIIS+U (κ = (1/4,0,0))
///   A3 MP 2×1×1 activation |E(U)−E(0)| + e_dftu + absolute pin
///
/// Exit non-zero on failure. Target ≤ 60 s Release.

#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
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

// Compact high-spin H₂ (DZV-GTH-PADE) — same class as ot-uks-kpoints.
// Geometry authored in Bohr; the Frame carries Å.
ATV::cpu::Frame h2_toy(double L_bohr) {
    using namespace ATV;
    constexpr double b2a = units::kBohrInAngstrom;
    const double xs[2] = {0.5 * L_bohr - 0.7, 0.5 * L_bohr + 0.7};
    cpu::Frame frame(2);
    for (int i = 0; i < 2; ++i) {
        frame.span(kAtom_type)[i] = 1;
        frame.span(kAtom_x)[i] = xs[i] * b2a;
        frame.span(kAtom_y)[i] = 0.5 * L_bohr * b2a;
        frame.span(kAtom_z)[i] = 0.5 * L_bohr * b2a;
    }
    frame.box = cpu::Box::cubic(L_bohr * b2a);
    return frame;
}

std::vector<ATV::cpu::gpw::Kind> h_kinds() {
    using namespace ATV;
    return {cpu::gpw::Kind{
        molrs::Element::H,
        cpu::gpw::GaussianBasis(molrs::Element::H, "DZV-GTH-PADE"),
        cpu::gpw::GthPotential(molrs::Element::H, "GTH-PADE-q1")}};
}

ATV::cpu::gpw::GpwDft base_uks_u(double u_ev) {
    using namespace ATV;
    cpu::gpw::GpwDft cfg{h_kinds()};
    cfg.cutoff_ry = 100.0;
    cfg.n_occ = 0;
    cfg.n_alpha = 2;
    cfg.n_beta = 0;
    cfg.max_iter = 120;
    cfg.tol = 1e-7;
    cfg.mix = 0.35;
    cfg.xc = cpu::gpw::XcKind::PADE;
    cfg.plus_u.push_back(cpu::gpw::PlusU{molrs::Element::H, /*l=*/0, u_ev});
    return cfg;
}

} // namespace

int main() {
    using namespace ATV;
    using cpu::gpw::GpwDft;
    using cpu::gpw::GpwResult;
    using cpu::gpw::PlusU;
    using cpu::gpw::XcKind;
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();

    const double L = 8.0;
    cpu::Frame frame = h2_toy(L);
    constexpr double kU_ev = 4.0;

    // ── A1: explicit Γ OT+U ≡ multi-k DIIS+U at Γ ──────────────────────────
    {
        GpwDft diis = base_uks_u(kU_ev);
        diis.use_ot = false;
        diis.k_points = {{0.0, 0.0, 0.0}};
        diis.k_point_weights = {1.0};
        diis.build(frame);
        const GpwResult rd = diis.run(frame);

        GpwDft ot = base_uks_u(kU_ev);
        ot.use_ot = true;
        ot.k_points = {{0.0, 0.0, 0.0}};
        ot.k_point_weights = {1.0};
        ot.max_iter = 200;
        ot.build(frame);
        const GpwResult ro = ot.run(frame);

        check(rd.converged && ro.converged, "A1 both converge");
        const double d = std::fabs(ro.e_total - rd.e_total);
        std::printf("  A1 E_DIIS+U=%.12f  E_OT+U=%.12f  |Δ|=%.3e  "
                    "e_dftu_ot=%.6e  s2=%.6f\n",
                    rd.e_total, ro.e_total, d, ro.e_dftu, ro.s2);
        check(d <= 1e-4, "A1 |E_OT+U − E_DIIS+U| ≤ 1e-4 Ha");
        check(std::fabs(ro.e_dftu) > 1e-8, "A1 e_dftu non-trivial");
        check(std::fabs(ro.e_dftu - rd.e_dftu) <= 1e-3,
              "A1 |e_dftu OT − DIIS| ≤ 1e-3 Ha");
    }

    // ── A2: multi-k OT+U ≡ multi-k DIIS+U (κ = 1/4) ────────────────────────
    double e_ot_a2 = 0.0;
    {
        GpwDft diis = base_uks_u(kU_ev);
        diis.use_ot = false;
        diis.k_points.push_back({0.25, 0.0, 0.0});
        diis.k_point_weights.push_back(1.0);
        diis.build(frame);
        const GpwResult rd = diis.run(frame);
        check(rd.converged, "A2 multi-k DIIS+U converged");

        GpwDft ot = base_uks_u(kU_ev);
        ot.use_ot = true;
        ot.k_points.push_back({0.25, 0.0, 0.0});
        ot.k_point_weights.push_back(1.0);
        ot.max_iter = 200;
        ot.build(frame);
        const GpwResult ro = ot.run(frame);
        check(ro.converged, "A2 multi-k OT+U converged");
        check(ro.n_kpts == 1, "A2 n_kpts == 1");
        e_ot_a2 = ro.e_total;
        const double d = std::fabs(ro.e_total - rd.e_total);
        std::printf("A2: E_OT+U = %.12f  E_DIIS+U = %.12f  |Δ|=%.3e Ha  "
                    "e_dftu=%.8f  iters ot/diis=%d/%d  s2=%.6f\n",
                    ro.e_total, rd.e_total, d, ro.e_dftu, ro.iters, rd.iters,
                    ro.s2);
        check(d <= 1e-4, "A2 |E_OT+U - E_DIIS+U| <= 1e-4 Ha");
        check(std::fabs(ro.e_dftu - rd.e_dftu) <= 1e-3,
              "A2 |e_dftu OT − DIIS| <= 1e-3 Ha");
        check(std::fabs(ro.e_dftu) > 1e-8, "A2 e_dftu non-trivial");
        check(std::fabs(ro.s2 - rd.s2) <= 1e-6, "A2 S2 OT ≡ DIIS");
    }

    // ── A3: MP 2×1×1 activation + absolute pin ─────────────────────────────
    // First-green ATV absolute (REFERENCE.md, 2026-07-26).
    // re-pinned 2026-07-31: variational +U scatter (pre-fix -0.802825557373)
    constexpr double kE_pin = -0.986384296380;
    constexpr double kTolE = 1e-3; // Ha
    constexpr double kFloorAct = 1e-4;

    {
        auto make = [&](bool with_u) {
            GpwDft m{h_kinds()};
            m.cutoff_ry = 100.0;
            m.n_occ = 0;
            m.n_alpha = 2;
            m.n_beta = 0;
            m.max_iter = 200;
            m.tol = 1e-7;
            m.mix = 0.35;
            m.xc = XcKind::PADE;
            m.k_mesh[0] = 2;
            m.k_mesh[1] = 1;
            m.k_mesh[2] = 1;
            m.k_time_reversal = true;
            m.use_ot = true;
            if (with_u)
                m.plus_u.push_back(PlusU{molrs::Element::H, /*l=*/0, kU_ev});
            m.build(frame);
            return m.run(frame);
        };

        const GpwResult r0 = make(false);
        const GpwResult ru = make(true);

        check(r0.converged && ru.converged, "A3 both OT converge");
        check(ru.n_kpts >= 1, "A3 n_kpts >= 1");
        check(std::isfinite(ru.e_total) && std::isfinite(ru.e_dftu),
              "A3 energies finite");
        check(std::fabs(ru.e_dftu) > 1e-8, "A3 e_dftu non-trivial");
        const double dU = std::fabs(ru.e_total - r0.e_total);
        check(dU >= kFloorAct, "A3 |E(U)−E(0)| activation floor");
        check(std::fabs(ru.e_total - kE_pin) <= kTolE,
              "A3 |E − pin class| within tol");

        std::printf("  E(U=0, OT) = %.12f Ha  iters=%d\n", r0.e_total,
                    r0.iters);
        std::printf("  E(U=4, OT) = %.12f Ha  iters=%d  e_dftu=%.8f\n",
                    ru.e_total, ru.iters, ru.e_dftu);
        std::printf("  |E(U)−E(0)| = %.6e Ha  A2_E=%.12f\n", dU, e_ot_a2);
    }

    const auto t1 = clock::now();
    const double wall = std::chrono::duration<double>(t1 - t0).count();
    std::printf("wall = %.3f s  fails=%d\n", wall, g_fail);
    if (g_fail != 0) {
        std::fprintf(stderr, "reg_ot-uks-dftu: %d check(s) failed\n", g_fail);
        return 1;
    }
    std::printf("reg_ot-uks-dftu: PASS\n");
    return 0;
}
