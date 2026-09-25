/// @file regressions/gpw/periodic-rpa-finite-q/reg_periodic-rpa-finite-q.cpp
/// @brief END-TO-END: finite-q inter-k RI-GPW dRPA on diamond Si (bct 4-atom,
///        MP 2×1×1). See REFERENCE.md. Public API only. ATV_BUILD_GPW
///        (+ HF + PW). No live CP2K. <= 60 s.

#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/gth_potential.h"
#include "atomiverse/cpu/gpw/rpa.h"
#include "atomiverse/molrs.hpp"

namespace {

// Diamond Si body-centred-tetragonal 4-atom cell (same as periodic-rpa-kpoints).
ATV::cpu::Frame si_bct(double a_ang) {
    using namespace ATV;
    const double at = a_ang / std::sqrt(2.0);
    const double frac[4][3] = {
        {0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}, {0.5, 0.0, 0.25}, {0.0, 0.5, 0.75}};
    cpu::Frame frame(4);
    for (int i = 0; i < 4; ++i) {
        frame.span(kAtom_type)[i] = 14;
        frame.span(kAtom_x)[i] = frac[i][0] * at;
        frame.span(kAtom_y)[i] = frac[i][1] * at;
        frame.span(kAtom_z)[i] = frac[i][2] * a_ang;
    }
    frame.box = cpu::Box::orthorhombic(at, at, a_ang);
    return frame;
}

ATV::cpu::gpw::Kind si_kind() {
    using namespace ATV::cpu::gpw;
    return Kind{molrs::Element::Si,
                GaussianBasis(molrs::Element::Si, "DZVP-MOLOPT-GTH"),
                GthPotential(molrs::Element::Si, "GTH-PBE-q4")};
}

} // namespace

int main() {
    using namespace ATV::cpu::gpw;
    const double a_ang = 5.43;

    // Offline pin strategy (no CP2K binary on worker host):
    //  - multi-k 2×1×1 hard gates + finite-q spectrum size N_q*N_k*nov
    //  - unit tests pin N_k=1 multi-k ≡ real plasmon
    //  - first-green absolute band (REFERENCE.md).
    // First-green pin 2026-07-26 (OMP=4); setup band covers orbital scatter.
    constexpr double kEcorr2kRef = -0.0736273682;
    constexpr double kEcorr2kTol = 5.0e-3;
    constexpr double kIdTol = 1.0e-12;

    ATV::cpu::Frame frame = si_bct(a_ang);

    GpwDft m{{si_kind()}};
    m.cutoff_ry = 80.0;
    m.n_occ = 8;
    m.max_iter = 60;
    m.tol = 1e-7;
    m.mix = 0.3;
    m.level_shift = 0.0;
    m.eps_default = 1e-7;
    m.xc = XcKind::PBE;
    m.mp2_ri = true;
    m.mp2_n_virt = 2;
    m.k_mesh[0] = 2;
    m.k_mesh[1] = 1;
    m.k_mesh[2] = 1;
    m.k_time_reversal = false;
    m.diis_depth = 6;
    m.build(frame);

    const auto t0 = std::chrono::steady_clock::now();
    const GpwResult& scf = m.run(frame);
    if (!scf.converged || scf.n_kpts != 2 || scf.C_k.size() != 2) {
        std::fprintf(stderr,
                     "FAIL: 2k SCF (conv=%d n_kpts=%d C_k=%zu iters=%d "
                     "E=%.8f)\n",
                     static_cast<int>(scf.converged), scf.n_kpts,
                     scf.C_k.size(), scf.iters, scf.e_total);
        return 1;
    }

    const RpaResult r = m.direct_rpa();
    const auto t1 = std::chrono::steady_clock::now();
    const double wall = std::chrono::duration<double>(t1 - t0).count();

    if (r.n_kpts != 2 || !r.ri) {
        std::fprintf(stderr, "FAIL: RPA flags (n_kpts=%d ri=%d)\n", r.n_kpts,
                     static_cast<int>(r.ri));
        return 1;
    }
    if (!std::isfinite(r.e_corr) || !(r.e_corr < 0.0)) {
        std::fprintf(stderr, "FAIL: E_corr=%.10f not finite/negative\n",
                     r.e_corr);
        return 1;
    }
    if (std::abs(r.e_total - (r.e_scf + r.e_corr)) > kIdTol) {
        std::fprintf(stderr, "FAIL: e_total != e_scf+e_corr\n");
        return 1;
    }
    if (r.n_occ != 8 || r.n_virt != 2) {
        std::fprintf(stderr, "FAIL: n_occ=%d n_virt=%d (want 8/2)\n", r.n_occ,
                     r.n_virt);
        return 1;
    }
    // Finite-q: n_omega = N_q * N_k * n_occ * n_virt with N_q = N_k
    const int n_omega_expect = 2 * 2 * 8 * 2;
    if (static_cast<int>(r.omega_ha.size()) != n_omega_expect) {
        std::fprintf(stderr, "FAIL: omega size %zu != %d (finite-q product)\n",
                     r.omega_ha.size(), n_omega_expect);
        return 1;
    }
    for (double om : r.omega_ha) {
        if (!(om > 0.0) || !std::isfinite(om)) {
            std::fprintf(stderr, "FAIL: non-positive Ω=%.6e\n", om);
            return 1;
        }
    }
    if (std::abs(r.e_corr - kEcorr2kRef) > kEcorr2kTol) {
        std::fprintf(stderr,
                     "FAIL: E_corr=%.10f outside loose band %.3f +/- %.3f "
                     "(update pin after first green)\n",
                     r.e_corr, kEcorr2kRef, kEcorr2kTol);
        return 1;
    }
    if (wall > 60.0) {
        std::fprintf(stderr, "FAIL: wall=%.1f s exceeds 60 s budget\n", wall);
        return 1;
    }

    std::printf(
        "PASS: finite-q inter-k RI-GPW dRPA Si4-bct (MP 2x1x1, n_virt=%d, "
        "%.0f Ry)\n"
        "  E_SCF   = %.10f Ha (iters=%d, n_kpts=%d)\n"
        "  E_corr  = %.10f Ha (Tr A=%.6f; n_occ=%d n_virt=%d)\n"
        "  E_total = %.10f Ha   wall=%.1f s\n"
        "  Ω_min   = %.6f Ha    Ω_max=%.6f Ha  (n_om=%zu)\n"
        "  [FIRST_GREEN_PIN e_corr=%.12f — paste into REFERENCE.md]\n",
        r.n_virt, m.cutoff_ry, r.e_scf, scf.iters, r.n_kpts, r.e_corr, r.tr_A,
        r.n_occ, r.n_virt, r.e_total, wall, r.omega_ha.front(),
        r.omega_ha.back(), r.omega_ha.size(), r.e_corr);
    return 0;
}
