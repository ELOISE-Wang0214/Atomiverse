/// @file regressions/gpw/ot-kpoints/reg_ot-kpoints.cpp
/// @brief END-TO-END: multi-k GPW RKS OT vs DIIS identities + frozen pin.
///
/// Public API only. See REFERENCE.md. Requires ATV_BUILD_GPW (+ HF + PW).
/// ctest never launches CP2K.

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
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

// Compact periodic H₂ (DZV-GTH-PADE) — same class as gpw-kpoints-si A1 toy.
// Cubic box of edge L [Bohr]; atoms ±0.7 Bohr about the centre.
ATV::cpu::Frame make_h2(double L) {
    using namespace ATV;
    cpu::Frame frame(2);
    const double c = 0.5 * L * units::kBohrInAngstrom;
    const double dx = 0.7 * units::kBohrInAngstrom;
    for (int i = 0; i < 2; ++i) {
        frame.span(kAtom_type)[i] = 1;
        frame.span(kAtom_y)[i] = c;
        frame.span(kAtom_z)[i] = c;
    }
    frame.span(kAtom_x)[0] = c - dx;
    frame.span(kAtom_x)[1] = c + dx;
    frame.box = cpu::Box::cubic(L * units::kBohrInAngstrom);
    return frame;
}

ATV::cpu::gpw::GpwDft base_method() {
    using namespace ATV;
    cpu::gpw::GpwDft m({cpu::gpw::Kind{
        molrs::Element::H,
        cpu::gpw::GaussianBasis(molrs::Element::H, "DZV-GTH-PADE"),
        cpu::gpw::GthPotential(molrs::Element::H, "GTH-PADE-q1")}});
    m.cutoff_ry = 100.0;
    m.n_occ = 1;
    m.max_iter = 120;
    m.tol = 1e-8;
    m.xc = ATV::cpu::gpw::XcKind::PADE;
    m.eps_default = 1e-7;
    return m;
}

} // namespace

int main() {
    using namespace ATV;
    using cpu::gpw::GpwDft;
    using cpu::gpw::XcKind;

    const double L = 8.0; // Bohr
    ATV::cpu::Frame frame = make_h2(L);

    // ── A1: explicit-Γ OT ≡ real Γ OT ──────────────────────────────────────
    {
        GpwDft real_ot = base_method();
        real_ot.use_ot = true;
        real_ot.max_iter = 200;
        real_ot.build(frame);
        const auto& rr = real_ot.run(frame);
        check(rr.converged, "A1 real Γ OT converged");

        GpwDft k_ot = base_method();
        k_ot.use_ot = true;
        k_ot.k_points.push_back({0.0, 0.0, 0.0});
        k_ot.k_point_weights.push_back(1.0);
        k_ot.max_iter = 200;
        k_ot.build(frame);
        const auto& rk = k_ot.run(frame);
        check(rk.converged, "A1 k-path Γ OT converged");
        const double d = std::abs(rr.e_total - rk.e_total);
        std::printf("A1: real Γ OT = %.10f  k-Γ OT = %.10f  |Δ|=%.3e Ha\n",
                    rr.e_total, rk.e_total, d);
        check(d <= 1e-6, "A1 |E_realOT - E_kGammaOT| <= 1e-6 Ha");
    }

    // ── A2: multi-k OT vs multi-k DIIS (κ = (1/4,0,0)) ─────────────────────
    double e_ot = 0.0;
    double e_diis = 0.0;
    {
        GpwDft diis = base_method();
        diis.k_points.push_back({0.25, 0.0, 0.0});
        diis.k_point_weights.push_back(1.0);
        diis.use_ot = false;
        diis.build(frame);
        const auto& rd = diis.run(frame);
        check(rd.converged, "A2 multi-k DIIS converged");
        e_diis = rd.e_total;

        GpwDft ot = base_method();
        ot.k_points.push_back({0.25, 0.0, 0.0});
        ot.k_point_weights.push_back(1.0);
        ot.use_ot = true;
        ot.max_iter = 200;
        ot.build(frame);
        const auto& ro = ot.run(frame);
        check(ro.converged, "A2 multi-k OT converged");
        e_ot = ro.e_total;
        check(ro.n_kpts == 1, "A2 n_kpts == 1");
        check(!ro.k_frac.empty() && std::abs(ro.k_frac[0][0] - 0.25) < 1e-12,
              "A2 k_frac payload");
        const double d = std::abs(e_ot - e_diis);
        std::printf("A2: E_OT = %.10f  E_DIIS = %.10f  |Δ|=%.3e Ha  "
                    "iters ot/diis=%d/%d\n",
                    e_ot, e_diis, d, ro.iters, rd.iters);
        check(d <= 1e-4, "A2 |E_OT - E_DIIS| <= 1e-4 Ha");
        check(std::isfinite(e_ot), "A2 finite E_OT");
    }

    // ── A3: Monkhorst–Pack 2×1×1 TR mesh OT vs DIIS ────────────────────────
    {
        GpwDft diis = base_method();
        diis.k_mesh[0] = 2;
        diis.k_mesh[1] = 1;
        diis.k_mesh[2] = 1;
        diis.k_time_reversal = true;
        diis.use_ot = false;
        diis.build(frame);
        const auto& rd = diis.run(frame);
        check(rd.converged, "A3 MP DIIS converged");

        GpwDft ot = base_method();
        ot.k_mesh[0] = 2;
        ot.k_mesh[1] = 1;
        ot.k_mesh[2] = 1;
        ot.k_time_reversal = true;
        ot.use_ot = true;
        ot.max_iter = 200;
        ot.build(frame);
        const auto& ro = ot.run(frame);
        check(ro.converged, "A3 MP OT converged");
        const double d = std::abs(ro.e_total - rd.e_total);
        std::printf("A3: MP OT = %.10f  MP DIIS = %.10f  |Δ|=%.3e Ha  "
                    "n_kpts=%d\n",
                    ro.e_total, rd.e_total, d, ro.n_kpts);
        check(d <= 1e-4, "A3 |E_OT - E_DIIS| MP mesh <= 1e-4 Ha");
        check(ro.n_kpts >= 1, "A3 n_kpts >= 1");
    }

    // ── A4: multi-k UKS OT is productized (reg_ot-uks-kpoints); RKS slug
    //        only checks that build accepts the combo without throw ───────
    {
        bool ok = true;
        try {
            GpwDft m = base_method();
            m.n_occ = 0;
            m.n_alpha = 2;
            m.n_beta = 0;
            m.use_ot = true;
            m.k_mesh[0] = 2;
            m.k_mesh[1] = m.k_mesh[2] = 1;
            m.build(frame);
        } catch (const std::exception& e) {
            ok = false;
            std::fprintf(stderr, "A4 unexpected: %s\n", e.what());
        }
        check(ok, "A4 multi-k UKS OT builds (product path)");
    }

    // ── A5: frozen first-green absolute band (until offline CP2K pin) ──────
    // Measured 2026-07-26; replace with CP2K OT+KPOINTS when available
    // (tools/cp2k_oracle/inputs/ot-kpoints/).
    {
        constexpr double kFrozenOtHa = -1.1080204833;
        constexpr double kAbsTol = 1e-3;
        const double d = std::abs(e_ot - kFrozenOtHa);
        std::printf("A5: multi-k OT E = %.10f Ha  frozen = %.10f  |Δ|=%.3e\n",
                    e_ot, kFrozenOtHa, d);
        check(d <= kAbsTol, "A5 |E_OT - frozen first-green| <= 1e-3 Ha");
    }

    if (g_fail) {
        std::fprintf(stderr, "FAIL ot-kpoints (%d checks)\n", g_fail);
        return 1;
    }
    std::printf("PASS: gpw ot-kpoints (multi-k OT identities)\n");
    return 0;
}
