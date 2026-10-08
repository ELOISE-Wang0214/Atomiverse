/// @file regressions/gpw/gpw-d4-kpoint-stress/reg_gpw-d4-kpoint-stress.cpp
/// @brief END-TO-END: multi-k GPW PBE-D4 analytic cell stress vs offline
///        CP2K Quickstep NUMERICAL-stress twins (+ dftd4 E_disp fidelity).
///
/// Capability: the k-summed fixed-\(P_k\) stress Lagrangian rebuilds the
/// geometry-only \(E_{\mathrm{disp}}^{\mathrm{D4}}\) at every strained cell
/// (`energy_fixed_P_k` in `gpw_stress.cpp`). This case adds the dedicated
/// multi-k oracle pin (deps closed: `gpw-d4-stress` Γ pin,
/// `gpw-kpoint-stress` multi-k plain stress).
///
/// Hard pins:
///   1. `E_disp` — offline dftd4 3.7.0 (ATV D4 backend; k-mesh independent).
///   2. `σ_disp` — σ(D4) − σ(plain) from two offline CP2K 2026.1 multi-k
///      `STRESS_TENSOR NUMERICAL` runs (`libdftd4` build).
/// Decisive dual gate: product analytic fixed-P_k stress (`d4 = true`,
/// `k_mesh = {1,1,2}`, `stress_use_fd = false`) vs multi-SCF strain-FD
/// (`gpw_stress_fd`) on the same multi-k method (all three diagonal axes).
/// Fair wall: product path vs the CP2K D4 multi-k NUMERICAL-stress twin.
///
/// No live CP2K. See REFERENCE.md.

#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <utility>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/evaluate/evaluate.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gpw_evaluator.h"
#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/gth_potential.h"
#include "atomiverse/cpu/scf/dftd4.h"
#include "atomiverse/molrs.hpp"

using namespace ATV;

namespace {

// ── Offline dftd4 3.7.0 unstrained E_disp (REFERENCE.md) ──────────────────
// Periodic water dimer, 12 Å cube, PBE-D4 BJ, ATM off. Bit-identical to the
// Γ `gpw-d4-water` / `gpw-d4-stress` pin (geometry-only; k-mesh independent).
constexpr double kEDispDftd4 = -0.00099458842745846774; // Ha
constexpr double kEDispTol = 1.0e-8;                    // Ha

// ── Offline CP2K 2026.1 multi-k pins (libdftd4 build; REFERENCE.md) ───────
// σ(D4) − σ(plain), diagonal, from two multi-k `STRESS| Numerical stress
// tensor [bar]` tables converted at 2.9421015697e8 bar per Ha/Bohr³.
// Twin: FACTOR_S9_TERM 0, CALCULATE_C9_TERM .FALSE., D4_CUTOFF 30 / R_CUTOFF 20.
constexpr double kSigmaDispCp2k[3] = {-7.61267484462925e-08,
                                      -1.13641449922510e-09,
                                      +2.00931506270126e-09}; // Ha/Bohr³
// Residual floor: CP2K↔dftd4 D4 parameter/cutoff parity (~2.5e-10 measured).
constexpr double kSigmaDispTol = 2.0e-9; // Ha/Bohr³

// Analytic fixed-P_k vs multi-SCF strain-FD (both d4 = true, k 1×1×2).
constexpr double kAnFdTol = 2.0e-5;     // Ha/Bohr³
constexpr double kOffDiagTol = 1.0e-12; // Ha/Bohr³ (must be exactly 0)
constexpr double kDiagFloor = 1.0e-10;  // Ha/Bohr³
constexpr double kWireTol = 2.0e-9;     // Ha/Bohr³ product σ(d4)−σ(plain)
constexpr double kWallBudget = 60.0;    // s product path
// Fair wall: CP2K internal timer on D4 multi-k NUMERICAL-stress twin
// (1 MPI × 8 OMP, 2026-08-01, private libdftd4 rebuild git:5e54ba2).
constexpr double kOracleWallS = 94.525;

int g_fail = 0;

void check(bool ok, const char* what) {
    if (!ok) {
        std::fprintf(stderr, "FAIL: %s\n", what);
        ++g_fail;
    } else {
        std::printf("  ok: %s\n", what);
    }
}

} // namespace

int main() {
    using clock = std::chrono::steady_clock;

    // Periodic water dimer, 12 Å cubic cell — same as gpw-d4-stress /
    // gpw-d3bj-kpoint-stress.
    constexpr double L = 12.0; // Å
    const int z[6] = {8, 1, 1, 8, 1, 1};
    const double xs[6] = {0.000000, 0.758602, 0.758602,
                          3.000000, 3.758602, 3.758602};
    const double ys[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    const double zs[6] = {0.000000, 0.504284, -0.504284,
                          0.000000, 0.504284, -0.504284};

    // ── (1) Product multi-k GPW evaluator with D4, analytic stress ────────
    cpu::Frame frame(6);
    for (int i = 0; i < 6; ++i) {
        frame.span(kAtom_type)[i] = z[i];
        frame.span(kAtom_x)[i] = xs[i];
        frame.span(kAtom_y)[i] = ys[i];
        frame.span(kAtom_z)[i] = zs[i];
    }
    frame.box = cpu::Box::cubic(L);

    const cpu::gpw::Kind O{
        molrs::Element::O,
        cpu::gpw::GaussianBasis(molrs::Element::O, "DZVP-MOLOPT-SR-GTH"),
        cpu::gpw::GthPotential(molrs::Element::O, "GTH-PBE-q6")};
    const cpu::gpw::Kind H{
        molrs::Element::H,
        cpu::gpw::GaussianBasis(molrs::Element::H, "DZVP-MOLOPT-SR-GTH"),
        cpu::gpw::GthPotential(molrs::Element::H, "GTH-PBE-q1")};

    cpu::gpw::GpwDft eval_method{{O, H}};
    eval_method.cutoff_ry = 120.0;
    eval_method.max_iter = 120;
    eval_method.tol = 1e-6;
    eval_method.mix = 0.4;
    eval_method.xc = cpu::gpw::XcKind::PBE;
    eval_method.nonperiodic = false;
    eval_method.k_mesh[0] = 1;
    eval_method.k_mesh[1] = 1;
    eval_method.k_mesh[2] = 2; // MP 1×1×2, TR-folded
    eval_method.d4 = true;
    cpu::gpw::GPWEvaluator eval(std::move(eval_method));
    eval.compute_forces = false;
    eval.compute_stress = true;
    eval.stress_use_fd = false;
    eval.method().stress_strain_h = 2.0e-3;

    const auto t0 = clock::now();
    cpu::Evaluate driver;
    driver.on(evaluate::Phase::evaluate).add_ref(eval);
    driver.enable_timing(false);
    driver.build();
    driver.run(frame);
    const auto t1 = clock::now();
    const double wall_product_s =
        std::chrono::duration<double>(t1 - t0).count();

    check(eval.last_result().converged, "A multi-k SCF converged");
    check(eval.last_result().n_kpts >= 1 && !eval.last_result().C_occ_k.empty(),
          "E multi-k C_occ_k payload present");
    const double e_disp = eval.last_result().e_disp;
    std::printf("A E_disp(multi-k) = %.14e Ha (dftd4 pin %.14e, |d|=%.2e)\n",
                e_disp, kEDispDftd4, std::abs(e_disp - kEDispDftd4));
    check(std::isfinite(e_disp) && std::abs(e_disp - kEDispDftd4) <= kEDispTol,
          "A E_disp equals the dftd4 pin (k-mesh independent)");

    const auto& san = eval.last_stress_ha_bohr3();
    for (int k = 3; k < 6; ++k)
        check(std::abs(san[static_cast<std::size_t>(k)]) <= kOffDiagTol,
              "E orthorhombic off-diagonal exactly 0");

    // ── (2) Multi-SCF strain-FD dual gate on the same multi-k method ──────

    auto make_cfg = [&](bool d4) {
        cpu::gpw::GpwDft m{{O, H}};
        m.cutoff_ry = 120.0;
        m.n_occ = 8;
        m.max_iter = 120;
        m.tol = 1e-6;
        m.mix = 0.4;
        m.xc = cpu::gpw::XcKind::PBE;
        m.nonperiodic = false;
        m.local_pp_on_grid = true;
        m.k_mesh[0] = 1;
        m.k_mesh[1] = 1;
        m.k_mesh[2] = 2;
        m.d4 = d4;
        return m;
    };

    const auto t2 = clock::now();
    // gpw_stress_fd returns eV/Å³; this twin compares in Ha/Bohr³.
    const auto sfd_ev_a3 =
        cpu::gpw::gpw_stress_fd(frame, make_cfg(true), 2.0e-3);
    std::array<double, 6> sfd{};
    for (std::size_t k = 0; k < 6; ++k)
        sfd[k] = sfd_ev_a3[k] / units::kHartreePerBohr3InEvPerAngstrom3;
    const auto t3 = clock::now();
    const double wall_fd_s = std::chrono::duration<double>(t3 - t2).count();

    double max_diag_err = 0.0;
    double max_diag = 0.0;
    const char* names[3] = {"xx", "yy", "zz"};
    for (int a = 0; a < 3; ++a) {
        const double an = san[static_cast<std::size_t>(a)];
        const double fd = sfd[static_cast<std::size_t>(a)];
        check(std::isfinite(an) && std::isfinite(fd),
              "C finite stress components");
        max_diag_err = std::max(max_diag_err, std::abs(an - fd));
        max_diag = std::max(max_diag, std::abs(an));
        std::printf("  sigma_%s  analytic=% .8e  fd=% .8e  |d|=%.3e Ha/Bohr3\n",
                    names[a], an, fd, std::abs(an - fd));
    }
    check(max_diag > kDiagFloor, "C stress diagonal nonzero");
    std::printf("C max|sigma_an - sigma_fd| = %.3e Ha/Bohr3\n", max_diag_err);
    check(max_diag_err <= kAnFdTol, "C analytic vs multi-SCF FD diagonal");

    // ── (3) Wiring: σ(d4) − σ(plain) vs the offline CP2K multi-k Δσ pin ───
    cpu::gpw::GpwDft plain = make_cfg(false);
    plain.stress_strain_h = 2.0e-3;
    cpu::Frame plain_frame = frame; // stress() writes kStress / kEnergy
    plain.build(plain_frame);
    std::array<double, 6> s_plain = plain.stress(plain_frame);
    for (auto& v : s_plain)
        v /= units::kHartreePerBohr3InEvPerAngstrom3; // eV/Å³ → Ha/Bohr³
    double max_wire = 0.0;
    double max_cp2k = 0.0;
    for (int a = 0; a < 3; ++a) {
        const double shift = san[static_cast<std::size_t>(a)]
                             - s_plain[static_cast<std::size_t>(a)];
        const double d_wire = std::abs(shift - kSigmaDispCp2k[a]);
        max_wire = std::max(max_wire, d_wire);
        max_cp2k = std::max(max_cp2k, d_wire);
        std::printf("  sigma_%s  shift=% .6e  cp2k_pin=% .6e  |d|=%.3e\n",
                    names[a], shift, kSigmaDispCp2k[a], d_wire);
    }
    check(max_cp2k <= kSigmaDispTol,
          "D sigma(d4)-sigma(plain) equals the CP2K multi-k Delta-sigma pin");
    check(max_wire <= kWireTol, "D wiring residual within product tol");

    // ── (4) Frame kStress unit conversion (eV/Å³) ─────────────────────────
    // kStress is the extensive Voigt virial [eV]: the method's Ha/Bohr^3
    // return stays intensive, so the published value is sigma * V.
    const double kStressVolA3 = frame.box.volume();
    const auto& st = frame.meta.get(kStress);
    const double ha_bohr3_to_ev_a3 =
        units::kHartreeInEv
        / (units::kBohrInAngstrom * units::kBohrInAngstrom
           * units::kBohrInAngstrom);
    bool conv_ok = true;
    for (int k = 0; k < 3; ++k) {
        const double expect =
            san[static_cast<std::size_t>(k)] * ha_bohr3_to_ev_a3 * kStressVolA3;
        if (std::abs(st[static_cast<std::size_t>(k)] - expect)
            > 1e-9 * kStressVolA3)
            conv_ok = false;
    }
    check(conv_ok, "E kStress unit conversion");

    check(wall_product_s <= kWallBudget, "product wall <= 60 s");
    const double ratio = wall_product_s / kOracleWallS;
    std::printf("fair_wall ratio T_atv/T_cp2k = %.3f  (cp2k=%.3f s  "
                "product=%.3f s  fd_twin=%.3f s)\n",
                ratio, kOracleWallS, wall_product_s, wall_fd_s);
    check(ratio <= 10.0, "fair wall ratio T_atv/T_cp2k <= 10");

    std::printf("gpw-d4-kpoint-stress  |dE_disp|=%.2e Ha  "
                "max|an-fd|=%.3e  max|wire-cp2k|=%.3e Ha/Bohr3  iters=%d  "
                "nk=%d  wall_product=%.2f s  ratio=%.3f\n",
                std::abs(e_disp - kEDispDftd4), max_diag_err, max_cp2k,
                eval.last_result().iters, eval.last_result().n_kpts,
                wall_product_s, ratio);

    if (g_fail) {
        std::fprintf(stderr, "FAILED: %d gate(s)\n", g_fail);
        return 1;
    }
    std::printf("PASS: gpw-d4-kpoint-stress\n");
    return 0;
}
