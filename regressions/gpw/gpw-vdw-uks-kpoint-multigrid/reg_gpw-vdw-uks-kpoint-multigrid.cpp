/// @file regressions/gpw/gpw-vdw-uks-kpoint-multigrid/reg_gpw-vdw-uks-kpoint-multigrid.cpp
/// @brief Multi-k collinear UKS GPW PBE + nonlocal vdW-DF on REL_CUTOFF multigrid.
///
/// Out-of-library consumer. Product pin for the triple combo already wired in
/// `build_uks_k_energy` / DIIS twin: V_nl of the k-summed total density is
/// added to both spin potentials on the fine grid before per-level restriction.
/// Gates (REFERENCE.md):
///   G1 multi-k UKS single-grid ↔ multigrid consistency (E_tot + E_disp);
///   G2 multi-k UKS multigrid vs offline CP2K 2026.1 LSD + DRSLL + MP + NGRIDS 4.
/// No live oracle in ctest.

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

// Same water dimer as gpw-vdw-uks-kpoints / gpw-vdw-multigrid [Å].
const double kAng[6][3] = {
    {0.000000, 0.000000, 0.000000},  // O
    {0.758602, 0.000000, 0.504284},  // H
    {0.758602, 0.000000, -0.504284}, // H
    {3.000000, 0.000000, 0.000000},  // O
    {3.758602, 0.000000, 0.504284},  // H
    {3.758602, 0.000000, -0.504284}, // H
};
const int kZ[6] = {8, 1, 1, 8, 1, 1};

ATV::cpu::Frame make_frame() {
    using namespace ATV;
    cpu::Frame frame(6);
    for (int i = 0; i < 6; ++i) {
        frame.span(kAtom_type)[i] = kZ[i];
        frame.span(kAtom_x)[i] = kAng[i][0];
        frame.span(kAtom_y)[i] = kAng[i][1];
        frame.span(kAtom_z)[i] = kAng[i][2];
    }
    frame.box = cpu::Box::cubic(12.0);
    return frame;
}

std::vector<ATV::cpu::gpw::Kind> water_kinds() {
    using namespace ATV::cpu::gpw;
    return {Kind{molrs::Element::O,
                 GaussianBasis(molrs::Element::O, "DZVP-MOLOPT-SR-GTH"),
                 GthPotential(molrs::Element::O, "GTH-PBE-q6")},
            Kind{molrs::Element::H,
                 GaussianBasis(molrs::Element::H, "DZVP-MOLOPT-SR-GTH"),
                 GthPotential(molrs::Element::H, "GTH-PBE-q1")}};
}

// Closed-shell UKS (nα = nβ = 8) water dimer; multi-k MP mesh fixed.
ATV::cpu::gpw::GpwDft make_uks(int ngrids, double rel_cutoff) {
    ATV::cpu::gpw::GpwDft gpw(water_kinds());
    gpw.cutoff_ry = 200.0;
    gpw.ngrids = ngrids;
    gpw.rel_cutoff = rel_cutoff;
    gpw.n_occ = 0;
    gpw.n_alpha = 8;
    gpw.n_beta = 8;
    gpw.max_iter = 60;
    gpw.tol = 1e-7;
    gpw.xc = ATV::cpu::gpw::XcKind::PBE;
    gpw.nonperiodic = false;
    gpw.vdw_df = true;
    gpw.k_mesh[0] = 2;
    gpw.k_mesh[1] = 1;
    gpw.k_mesh[2] = 1;
    // Match CP2K's full MP mesh (no TR fold) for a 1:1 twin.
    gpw.k_time_reversal = false;
    return gpw;
}

} // namespace

int main() {
    using namespace ATV::cpu::gpw;
    using Clock = std::chrono::steady_clock;

    // Offline CP2K 2026.1 (cp2k.psmp rev git:5e54ba2): water dimer LSD MULT 1
    // + NON_LOCAL DRSLL + MP 2×1×1 NGRIDS 4 / REL_CUTOFF 40 — REFERENCE.md.
    constexpr double kEDispCp2k = 0.15193559105517;   // Ha
    constexpr double kETotCp2k = -34.186333514352619; // Ha
    constexpr double kEDispTol = 1.0e-4;              // Ha
    // Same Γ-dimer twin-family systematic as gpw-vdw-kpoints / multigrid.
    constexpr double kETotTol = 1.0e-2; // Ha
    // Internal single-grid ↔ multigrid (gpw-vdw-multigrid uses 2e-5).
    constexpr double kEConsistencyTol = 2.0e-5; // Ha
    // CP2K timer TOTAL TIME on the same twin job (1 MPI × OMP=8).
    constexpr double kCp2kWall = 13.703; // s

    ATV::cpu::Frame frame = make_frame();

    // ── G1: multi-k UKS single-grid ↔ multigrid consistency ───────────────
    double e_tot_sg = 0.0, e_disp_sg = 0.0;
    {
        GpwDft sg = make_uks(/*ngrids=*/1, /*rel_cutoff=*/0.0);
        sg.build(frame);
        const GpwResult& rs = sg.run(frame);
        if (!rs.converged) {
            std::fprintf(stderr,
                         "FAIL: G1 multi-k UKS single-grid SCF did not "
                         "converge (%d)\n",
                         rs.iters);
            return 1;
        }
        e_tot_sg = rs.e_total;
        e_disp_sg = rs.e_disp;
        if (!(std::isfinite(e_disp_sg) && std::abs(e_disp_sg) > 1e-6)) {
            std::fprintf(stderr,
                         "FAIL: G1 expected non-trivial E_disp on single-grid, "
                         "got %.6e\n",
                         e_disp_sg);
            return 1;
        }
    }

    // ── G2: multigrid multi-k UKS vs offline CP2K (fair-wall twin) ────────
    double wall_g2 = 0.0;
    double e_tot_mg = 0.0, e_disp_mg = 0.0;
    {
        GpwDft mg = make_uks(/*ngrids=*/4, /*rel_cutoff=*/40.0);
        mg.build(frame);
        const auto t0 = Clock::now();
        const GpwResult& rm = mg.run(frame);
        wall_g2 = std::chrono::duration<double>(Clock::now() - t0).count();
        if (!rm.converged) {
            std::fprintf(stderr,
                         "FAIL: G2 multi-k UKS multigrid SCF did not converge "
                         "(%d)\n",
                         rm.iters);
            return 1;
        }
        e_tot_mg = rm.e_total;
        e_disp_mg = rm.e_disp;

        const double d_tot_sg = std::abs(e_tot_mg - e_tot_sg);
        const double d_disp_sg = std::abs(e_disp_mg - e_disp_sg);
        std::printf("G1: E_tot  sg=%.14f  mg=%.14f  |Δ|=%.3e Ha\n", e_tot_sg,
                    e_tot_mg, d_tot_sg);
        std::printf("G1: E_disp sg=%.14f  mg=%.14f  |Δ|=%.3e Ha\n", e_disp_sg,
                    e_disp_mg, d_disp_sg);
        if (d_tot_sg > kEConsistencyTol || d_disp_sg > kEConsistencyTol) {
            std::fprintf(stderr,
                         "FAIL: G1 mg↔sg |ΔE_tot|=%.3e |ΔE_disp|=%.3e > %.1e "
                         "Ha\n",
                         d_tot_sg, d_disp_sg, kEConsistencyTol);
            return 1;
        }

        std::printf("G2: E_tot  ATV=%.14f  CP2K=%.14f  |Δ|=%.3e Ha\n", e_tot_mg,
                    kETotCp2k, std::abs(e_tot_mg - kETotCp2k));
        std::printf("G2: E_disp ATV=%.14f  CP2K=%.14f  |Δ|=%.3e Ha\n",
                    e_disp_mg, kEDispCp2k, std::abs(e_disp_mg - kEDispCp2k));
        std::printf("G2: wall %.1f s (CP2K twin wall %.1f s OMP=8)  "
                    "T_atv/T_cp2k=%.2f\n",
                    wall_g2, kCp2kWall, wall_g2 / kCp2kWall);

        if (!(std::isfinite(e_disp_mg) && std::abs(e_disp_mg) > 1e-6)) {
            std::fprintf(stderr,
                         "FAIL: G2 expected non-trivial E_disp on multigrid, "
                         "got %.6e\n",
                         e_disp_mg);
            return 1;
        }
        if (std::abs(e_disp_mg - kEDispCp2k) > kEDispTol) {
            std::fprintf(stderr,
                         "FAIL: G2 E_disp got=%.14e want=%.14e ±%.1e Ha\n",
                         e_disp_mg, kEDispCp2k, kEDispTol);
            return 1;
        }
        if (std::abs(e_tot_mg - kETotCp2k) > kETotTol) {
            std::fprintf(stderr,
                         "FAIL: G2 E_tot got=%.14f want=%.14f ±%.1e Ha\n",
                         e_tot_mg, kETotCp2k, kETotTol);
            return 1;
        }
    }

    std::printf("PASS gpw-vdw-uks-kpoint-multigrid (G2 wall %.1f s)\n",
                wall_g2);
    return 0;
}
