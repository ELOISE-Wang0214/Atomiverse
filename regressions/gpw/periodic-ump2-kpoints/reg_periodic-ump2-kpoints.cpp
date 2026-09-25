/// @file regressions/gpw/periodic-ump2-kpoints/reg_periodic-ump2-kpoints.cpp
/// @brief END-TO-END: multi-k unrestricted RI-GPW MP2 (UMP2) on diamond Si
///        (bct 4-atom, MP 2×1×1, UKS nα=nβ). See REFERENCE.md.
///        Public API only. ATV_BUILD_GPW (+ HF + PW). No live CP2K. <= 60 s.

#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/gth_potential.h"
#include "atomiverse/cpu/gpw/mp2.h"
#include "atomiverse/molrs.hpp"

namespace {

// Diamond Si body-centred-tetragonal 4-atom cell (same as periodic-mp2-kpoints).
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

    // Dual-gate vs closed-shell multi-k RMP2 first-green
    // (regressions/gpw/periodic-mp2-kpoints): UKS nα=nβ → UMP2 ≡ RMP2 when
    // orbitals match. Offline CP2K re-pin when available — REFERENCE.md.
    constexpr double kEcorr2kRef = -0.0657637511; // Ha (RMP2 first-green twin)
    constexpr double kEcorr2kTol = 2.0e-3;        // Ha (UKS numerical drift)
    constexpr double kIdTol = 1.0e-12;

    ATV::cpu::Frame frame = si_bct(a_ang);

    GpwDft m{{si_kind()}};
    m.cutoff_ry = 80.0;
    m.n_alpha = 8; // 16 valence e⁻, closed-shell occupations
    m.n_beta = 8;
    m.max_iter = 60;
    m.tol = 1e-7;
    m.mix = 0.3;
    m.level_shift = 0.0; // k-path rejects nonzero
    m.eps_default = 1e-7;
    m.xc = XcKind::PBE;
    m.mp2_ri = true;
    m.mp2_n_virt = 2; // tight virt cap for ≤60 s
    m.k_mesh[0] = 2;
    m.k_mesh[1] = 1;
    m.k_mesh[2] = 1;
    m.k_time_reversal = false; // equal weights
    m.diis_depth = 6;
    m.build(frame);

    const auto t0 = std::chrono::steady_clock::now();
    const GpwResult& scf = m.run(frame);
    if (!scf.converged || scf.n_kpts != 2 || scf.C_k.size() != 2
        || scf.C_k_beta.size() != 2) {
        std::fprintf(stderr,
                     "FAIL: 2k UKS SCF (conv=%d n_kpts=%d C_k=%zu C_k_beta=%zu "
                     "iters=%d E=%.8f)\n",
                     static_cast<int>(scf.converged), scf.n_kpts,
                     scf.C_k.size(), scf.C_k_beta.size(), scf.iters,
                     scf.e_total);
        return 1;
    }
    if (scf.eps_k_beta.size() != 2) {
        std::fprintf(stderr, "FAIL: missing eps_k_beta (size=%zu)\n",
                     scf.eps_k_beta.size());
        return 1;
    }

    const Mp2Result r = m.mp2();
    const auto t1 = std::chrono::steady_clock::now();
    const double wall = std::chrono::duration<double>(t1 - t0).count();

    if (r.n_kpts != 2 || !r.ri) {
        std::fprintf(stderr, "FAIL: UMP2 flags (n_kpts=%d ri=%d)\n", r.n_kpts,
                     static_cast<int>(r.ri));
        return 1;
    }
    if (!std::isfinite(r.e_corr) || !(r.e_corr < 0.0)) {
        std::fprintf(stderr, "FAIL: E_corr=%.10f not finite/negative\n",
                     r.e_corr);
        return 1;
    }
    if (std::abs(r.e_corr_os + r.e_corr_ss - r.e_corr) > kIdTol) {
        std::fprintf(stderr,
                     "FAIL: OS+SS broken: os=%.12f ss=%.12f corr=%.12f\n",
                     r.e_corr_os, r.e_corr_ss, r.e_corr);
        return 1;
    }
    if (std::abs(r.e_total - (r.e_scf + r.e_corr)) > kIdTol) {
        std::fprintf(stderr, "FAIL: e_total != e_scf+e_corr\n");
        return 1;
    }
    // n_occ = nα+nβ = 16; n_virt = 2+2 = 4 under mp2_n_virt=2 per spin
    if (r.n_occ != 16 || r.n_virt != 4) {
        std::fprintf(stderr, "FAIL: n_occ=%d n_virt=%d (want 16/4)\n", r.n_occ,
                     r.n_virt);
        return 1;
    }
    if (std::abs(r.e_corr - kEcorr2kRef) > kEcorr2kTol) {
        std::fprintf(stderr,
                     "FAIL: E_corr=%.10f outside dual-gate band %.6f +/- %.3f "
                     "(RMP2 multi-k twin; see REFERENCE.md)\n",
                     r.e_corr, kEcorr2kRef, kEcorr2kTol);
        return 1;
    }
    if (wall > 60.0) {
        std::fprintf(stderr, "FAIL: wall=%.1f s exceeds 60 s budget\n", wall);
        return 1;
    }

    std::printf(
        "PASS: multi-k RI-GPW UMP2 Si4-bct (MP 2x1x1, UKS nα=nβ=8, n_virt=%d, "
        "%.0f Ry)\n"
        "  E_SCF   = %.10f Ha (iters=%d, n_kpts=%d, <S2>=%.6f)\n"
        "  E_corr  = %.10f Ha (OS=%.10f SS=%.10f)\n"
        "  E_total = %.10f Ha   wall=%.1f s\n"
        "  [multi-k UMP2 dual-gate vs RMP2 first-green; see REFERENCE.md]\n",
        r.n_virt, m.cutoff_ry, r.e_scf, scf.iters, r.n_kpts, scf.s2, r.e_corr,
        r.e_corr_os, r.e_corr_ss, r.e_total, wall);
    return 0;
}
