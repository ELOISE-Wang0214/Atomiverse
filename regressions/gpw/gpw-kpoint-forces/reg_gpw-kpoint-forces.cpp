/// @file regressions/gpw/gpw-kpoint-forces/reg_gpw-kpoint-forces.cpp
/// @brief END-TO-END: GPW analytic forces on a Monkhorst–Pack k mesh.
///
/// Periodic H₂ (DZV-GTH-PADE, PADE XC), k_mesh = 1×1×2, compressed bond.
/// See REFERENCE.md. Public: `GpwDft::forces(frame)`.

#include <array>
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

// Frozen first-green ATV (Release, OMP=8). Offline FD on related geometry:
// max|F_an − F_fd| = 5.98e-6 Ha/Bohr (h = 1e-3). See REFERENCE.md.
constexpr double kEPin = -1.057153829426;  // Ha
constexpr double kF0xPin = -0.23541914472; // Ha/Bohr
constexpr double kF1xPin = 0.23541914472;  // Ha/Bohr
constexpr double kETol = 1e-5;             // Ha
constexpr double kForceTol = 1e-5;         // Ha/Bohr
constexpr double kSumFTol = 1e-6;          // Ha/Bohr (symmetry residual)
constexpr double kNewtonTol = 1e-6;

} // namespace

int main() {
    using namespace ATV;
    using namespace ATV::cpu::gpw;

    const double rA[3] = {5.5, 6.0, 6.0}; // Bohr
    const double rB[3] = {6.5, 6.0, 6.0}; // 1.0 Bohr compressed bond
    cpu::Frame frame(2);
    for (int i = 0; i < 2; ++i) {
        const double* r = (i == 0) ? rA : rB;
        frame.span(kAtom_type)[i] = 1;
        frame.span(kAtom_x)[i] = r[0] * units::kBohrInAngstrom;
        frame.span(kAtom_y)[i] = r[1] * units::kBohrInAngstrom;
        frame.span(kAtom_z)[i] = r[2] * units::kBohrInAngstrom;
    }
    frame.box = cpu::Box::cubic(12.0 * units::kBohrInAngstrom);

    auto make_cfg = []() {
        GpwDft m{{Kind{molrs::Element::H,
                       GaussianBasis(molrs::Element::H, "DZV-GTH-PADE"),
                       GthPotential(molrs::Element::H, "GTH-PADE-q1")}}};
        m.cutoff_ry = 100.0;
        m.n_occ = 1;
        m.tol = 1e-8;
        m.max_iter = 60;
        m.nonperiodic = false;
        m.local_pp_on_grid = true;
        m.xc = XcKind::PADE;
        m.eps_default = 1e-8;
        m.k_mesh[0] = 1;
        m.k_mesh[1] = 1;
        m.k_mesh[2] = 2;
        return m;
    };

    GpwDft energy = make_cfg();
    energy.build(frame);
    const auto& er = energy.run(frame);
    check(er.converged, "A1 k-mesh SCF converged");
    check(std::isfinite(er.e_total), "A1 finite energy");
    std::printf("E(k 1x1x2) = %.12f Ha  iters=%d  nk=%d\n", er.e_total,
                er.iters, er.n_kpts);
    check(std::abs(er.e_total - kEPin) <= kETol, "A3 |E - pin| <= 1e-5 Ha");

    // Fresh method: SCF + analytic forces; Frame forces [eV/Å] → Ha/Bohr.
    GpwDft fm = make_cfg();
    fm.build(frame);
    fm.forces(frame);
    std::vector<std::array<double, 3>> Fa(2);
    for (int i = 0; i < 2; ++i)
        Fa[static_cast<std::size_t>(i)] = {
            frame.span(kAtom_fx)[i] / units::kHartreePerBohrInEvPerAngstrom,
            frame.span(kAtom_fy)[i] / units::kHartreePerBohrInEvPerAngstrom,
            frame.span(kAtom_fz)[i] / units::kHartreePerBohrInEvPerAngstrom};
    check(Fa.size() == 2, "A1 two force vectors");
    for (int a = 0; a < 2; ++a)
        for (int k = 0; k < 3; ++k)
            check(std::isfinite(Fa[static_cast<std::size_t>(a)][k]),
                  "A1 finite force component");

    std::printf("  F0 = (% .10e, % .10e, % .10e)\n", Fa[0][0], Fa[0][1],
                Fa[0][2]);
    std::printf("  F1 = (% .10e, % .10e, % .10e)\n", Fa[1][0], Fa[1][1],
                Fa[1][2]);

    check(Fa[0][0] < 0.0 && Fa[1][0] > 0.0, "A1 restoring bond force signs");
    check(std::abs(Fa[0][0] - kF0xPin) <= kForceTol, "A2 |F0_x - pin|");
    check(std::abs(Fa[1][0] - kF1xPin) <= kForceTol, "A2 |F1_x - pin|");
    check(std::abs(Fa[0][0] + Fa[1][0]) <= kNewtonTol, "A4 Newton-3 F_x");
    check(std::abs(Fa[0][1]) <= 1e-6 && std::abs(Fa[0][2]) <= 1e-6
              && std::abs(Fa[1][1]) <= 1e-6 && std::abs(Fa[1][2]) <= 1e-6,
          "A2 F_y/F_z ~ 0 by symmetry");

    double sum[3] = {Fa[0][0] + Fa[1][0], Fa[0][1] + Fa[1][1],
                     Fa[0][2] + Fa[1][2]};
    const double sum_mag =
        std::sqrt(sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]);
    std::printf("A4: |sum F| = %.3e Ha/Bohr\n", sum_mag);
    check(sum_mag <= kSumFTol, "A4 |sum F| <= 1e-6 Ha/Bohr");

    if (g_fail) {
        std::fprintf(stderr, "FAILED: %d gate(s)\n", g_fail);
        return 1;
    }
    std::printf("PASS: gpw-kpoint-forces  E=%.10f  F0x=%.8f  |sumF|=%.3e\n",
                er.e_total, Fa[0][0], sum_mag);
    return 0;
}
