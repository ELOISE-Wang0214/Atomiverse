/// @file regressions/gpw/gpw-uks-smearing/reg_gpw-uks-smearing.cpp
/// @brief END-TO-END: multi-k GPW collinear UKS + Fermi–Dirac metallic smearing.
///
/// See REFERENCE.md. Public headers only. Requires ATV_BUILD_GPW.
/// Gates: closed-shell UKS+smear ≡ RKS+smear (H₂ toy); BCC Na UKS FD free-energy
/// pin, electron count, finite Fermi level.

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

constexpr double kA2B = 1.0 / 0.52917720859; // Å → Bohr

int g_fail = 0;

void check(bool ok, const char* what) {
    if (!ok) {
        std::fprintf(stderr, "FAIL: %s\n", what);
        ++g_fail;
    } else {
        std::printf("  ok: %s\n", what);
    }
}

/// H₂ toy in a cubic box of edge L [Bohr]; atoms ±0.7 Bohr about the centre.
ATV::cpu::Frame h2_toy(double L_bohr) {
    using namespace ATV;
    cpu::Frame frame(2);
    const double c = 0.5 * L_bohr / kA2B;
    const double dz = 0.7 / kA2B;
    for (int i = 0; i < 2; ++i) {
        frame.span(kAtom_type)[i] = 1;
        frame.span(kAtom_y)[i] = c;
        frame.span(kAtom_z)[i] = c;
    }
    frame.span(kAtom_x)[0] = c - dz;
    frame.span(kAtom_x)[1] = c + dz;
    frame.box = cpu::Box::cubic(L_bohr / kA2B);
    return frame;
}

ATV::cpu::gpw::Kind h_kind() {
    using namespace ATV;
    return cpu::gpw::Kind{
        molrs::Element::H,
        cpu::gpw::GaussianBasis(molrs::Element::H, "DZV-GTH-PADE"),
        cpu::gpw::GthPotential(molrs::Element::H, "GTH-PADE-q1")};
}

/// BCC conventional cubic: 2 Na at (0,0,0) and (½,½,½); a in Å.
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

} // namespace

int main() {
    using namespace ATV;
    using cpu::gpw::GpwDft;
    using cpu::gpw::Smearing;
    using cpu::gpw::XcKind;

    // ── R5: closed-shell UKS+smear ≡ RKS+smear on H₂ multi-k ──────────────
    {
        const double L = 8.0;
        cpu::Frame frame = h2_toy(L);

        GpwDft rks({h_kind()});
        rks.cutoff_ry = 100.0;
        rks.n_occ = 1;
        rks.max_iter = 80;
        rks.tol = 1e-7;
        rks.xc = XcKind::PADE;
        rks.k_mesh[0] = rks.k_mesh[1] = 1;
        rks.k_mesh[2] = 2;
        rks.smearing = Smearing::fermi_dirac;
        rks.degauss = 0.02;
        rks.build(frame);
        const auto rr = rks.run(frame);

        GpwDft uks({h_kind()});
        uks.cutoff_ry = 100.0;
        uks.n_alpha = 1;
        uks.n_beta = 1;
        uks.max_iter = 80;
        uks.tol = 1e-7;
        uks.xc = XcKind::PADE;
        uks.k_mesh[0] = uks.k_mesh[1] = 1;
        uks.k_mesh[2] = 2;
        uks.smearing = Smearing::fermi_dirac;
        uks.degauss = 0.02;
        uks.build(frame);
        const auto ru = uks.run(frame);

        check(rr.converged && ru.converged, "R5 both converge");
        check(std::fabs(ru.e_total - rr.e_total) <= 5e-6,
              "R5 |F_UKS − F_RKS| ≤ 5e-6 Ha");
        check(std::fabs(ru.n_elec_occ - 2.0) < 1e-5
                  && std::fabs(rr.n_elec_occ - 2.0) < 1e-5,
              "R5 both recover Ne=2");
        std::printf("  R5 dF = %.3e Ha  F_rks=%.10f  F_uks=%.10f\n",
                    ru.e_total - rr.e_total, rr.e_total, ru.e_total);
    }

    // ── Primary pin: BCC Na multi-k UKS + FD (ferro seed) ─────────────────
    // Offline free-energy pin — frozen first green (REFERENCE.md).
    const double f_pin_ha = -92.126330608004; // Ha free energy F
    const double f_tol_ha = 1e-4;             // Ha
    const double ne_tol = 1e-5;
    const double a_ang = 4.225;
    const double n_elec = 18.0; // 2 × GTH-PBE-q9

    cpu::Frame frame = na_bcc(a_ang);

    GpwDft m({cpu::gpw::Kind{
        molrs::Element::Na,
        cpu::gpw::GaussianBasis(molrs::Element::Na, "DZVP-MOLOPT-SR-GTH"),
        cpu::gpw::GthPotential(molrs::Element::Na, "GTH-PBE-q9")}});
    m.cutoff_ry = 150.0;
    m.n_alpha = 10; // ferro seed; M free under smear
    m.n_beta = 8;
    m.max_iter = 80;
    m.tol = 1e-6;
    m.eps_default = 1e-7;
    m.xc = XcKind::PADE;
    m.diis_depth = 6;
    m.mix = 0.4;
    m.k_mesh[0] = m.k_mesh[1] = m.k_mesh[2] = 2;
    m.k_time_reversal = true;
    m.smearing = Smearing::fermi_dirac;
    m.degauss = 0.02; // Ry
    m.extra_bands = 4;

    m.build(frame);

    const auto t0 = std::chrono::steady_clock::now();
    const auto& r = m.run(frame);
    const auto t1 = std::chrono::steady_clock::now();
    const double wall_s = std::chrono::duration<double>(t1 - t0).count();

    std::printf("F=%.12f Ha  e_fermi=%.6f Ha  n_elec=%.8f  e_entropy=%.6e Ha  "
                "iters=%d  wall=%.3f s  converged=%d\n",
                r.e_total, r.e_fermi, r.n_elec_occ, r.e_entropy, r.iters,
                wall_s, r.converged ? 1 : 0);

    check(r.converged, "R1 SCF converged");
    check(std::isfinite(r.e_total) && std::isfinite(r.e_fermi),
          "finite free energy and e_fermi");
    check(std::fabs(r.n_elec_occ - n_elec) < ne_tol, "R2 electron count");
    check(std::fabs(r.e_fermi) < 1e3, "R3 Fermi level finite");
    check(std::fabs(r.e_total - f_pin_ha) < f_tol_ha, "R4 free energy pin");
    check(std::isfinite(r.e_entropy), "finite e_entropy");
    check(r.n_kpts > 1, "multi-k sampling active");

    if (g_fail) {
        std::printf("FAIL: gpw-uks-smearing (%d checks)\n", g_fail);
        return 1;
    }
    std::printf("PASS: gpw-uks-smearing (BCC Na UKS FD smear)\n");
    return 0;
}
