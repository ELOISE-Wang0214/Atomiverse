/// @file tests/cpu/scf/test_dispersion_gradient.cpp
/// @brief The molecular dispersion evaluators are conservative: `forces()` is
///        the analytic gradient of the SAME functional `energy()` reports.
///
/// This is the unit-level half of the "BOMD+d3bj/d4 forces are
/// non-conservative" fix. `MeanFieldMethod` folds E_disp into the total energy
/// it reports (`e_total_ += e_disp_`); `internal::ScfGradient` now folds
/// −∂E_disp/∂R into the force it reports. Whether the two describe one surface
/// reduces entirely to whether `D3BJ`/`D4` `energy()` and `forces()` do — which
/// is what these tests measure, with no SCF anywhere (unit tier: one method per
/// test, milliseconds).
///
/// The gate is analytic-vs-central-difference of the same functional, not a
/// literature number: `regressions/scf/d3bj-molecular-pbe` and
/// `regressions/scf/dft-d4-h2o-pbe` pin the *energies* against offline
/// simple-dftd3 / dftd4 oracles, and the end-to-end conservation witness is
/// `regressions/scf/bomd-d3bj-ar2-nve`.
///
/// Note on flavours: D3(BJ) has a closed-form analytic gradient (two-body
/// C6/C8 BJ terms plus the dC6/dCN·dCN/dR coordination-number chain) shared
/// with GPW. D4 two-body (`s9 = 0`) now has the matching closed-form path:
/// BJ pair + CN_d4 weight chain + EEQ charge-response adjoint. ATM (`s9 ≠ 0`)
/// still uses central differences of the full energy. The D4 rows below are
/// analytic-vs-FD of the same functional (dual-gated against dftd4 3.7.0 `-g`
/// in the GPW unit suite).

#ifdef ATV_HAS_HF

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>
#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/dftd3.h"
#include "atomiverse/cpu/scf/dftd4.h"
#include "atomiverse/cpu/scf/gradient_result.h"
#include "atomiverse/cpu/scf/internal/dispersion.h"

using ATV::cpu::Frame;
using ATV::cpu::scf::D3BJ;
using ATV::cpu::scf::D4;
using ATV::cpu::scf::GradientResult;

namespace {

using ATV::kAtom_type;
using ATV::kAtom_x;
using ATV::kAtom_y;
using ATV::kAtom_z;

/// @brief Water dimer [Å] — the geometry of the GPW D3 CP2K reference case.
/// @return Six-atom Frame carrying Z in `kAtom_type` and positions in Å.
Frame water_dimer() {
    const int z[6] = {8, 1, 1, 8, 1, 1};
    const double xyz[6][3] = {
        {0.000000, 0.000000, 0.000000},  // O
        {0.758602, 0.000000, 0.504284},  // H
        {0.758602, 0.000000, -0.504284}, // H
        {3.000000, 0.000000, 0.000000},  // O
        {3.758602, 0.000000, 0.504284},  // H
        {3.758602, 0.000000, -0.504284}, // H
    };
    Frame f(6);
    for (int i = 0; i < 6; ++i) {
        f.span(kAtom_type)[i] = z[i];
        f.span(kAtom_x)[i] = xyz[i][0];
        f.span(kAtom_y)[i] = xyz[i][1];
        f.span(kAtom_z)[i] = xyz[i][2];
    }
    return f;
}

/// @brief Displace one Cartesian component of one atom [Å].
/// @param f Frame to mutate.
/// @param atom Atom index.
/// @param d Cartesian direction 0/1/2.
/// @param delta Displacement [Å].
void displace(Frame& f, int atom, int d, double delta) {
    if (d == 0)
        f.span(kAtom_x)[atom] += delta;
    else if (d == 1)
        f.span(kAtom_y)[atom] += delta;
    else
        f.span(kAtom_z)[atom] += delta;
}

/// @brief Worst |F_analytic − F_central_difference| over every atom/direction.
///
/// Both sides are read from the SAME evaluator object, so this measures only
/// whether its `forces()` differentiates its `energy()`.
/// @param method  Dispersion evaluator (`D3BJ` or `D4`).
/// @param h_bohr  Central-difference half-step [Bohr].
/// @return Worst absolute residual [Hartree/Bohr].
template <typename Method>
double worst_fd_residual(const Method& method, double h_bohr) {
    Frame f = water_dimer();
    const std::vector<std::array<double, 3>> fa = method.forces(f);
    const double h_ang = h_bohr * ATV::units::kBohrInAngstrom;
    double worst = 0.0;
    for (int a = 0; a < 6; ++a) {
        for (int d = 0; d < 3; ++d) {
            displace(f, a, d, h_ang);
            const double ep = method.energy(f);
            displace(f, a, d, -2.0 * h_ang);
            const double em = method.energy(f);
            displace(f, a, d, h_ang);
            // F = −dE/dR, in Hartree/Bohr
            const double fd = -(ep - em) / (2.0 * h_bohr);
            worst = std::max(
                worst,
                std::abs(fd - fa[static_cast<std::size_t>(a)]
                                 [static_cast<std::size_t>(d)]));
        }
    }
    return worst;
}

/// @brief Largest force component of a dispersion evaluation [Hartree/Bohr].
/// @param f Forces as returned by an evaluator.
/// @return max |F| over atoms and directions.
double max_abs(const std::vector<std::array<double, 3>>& f) {
    double m = 0.0;
    for (const auto& t : f)
        for (int d = 0; d < 3; ++d)
            m = std::max(m, std::abs(t[static_cast<std::size_t>(d)]));
    return m;
}

/// @brief eV/Å per Hartree/Bohr — the public-surface scale of a residual.
constexpr double kHaBohrInEvAng = 1.0 / ATV::units::kEvPerAngstromInHartreePerBohr;

} // namespace

// ── D3(BJ): closed-form analytic gradient ─────────────────────────────────

/// The whole conservativeness claim for `d3bj`: the analytic force (two-body
/// BJ terms + the dC6/dCN chain) reproduces a central difference of the energy
/// this same object reports. A dropped CN chain shows up here at ~1e-4 Ha/Bohr,
/// four orders above this gate.
TEST(ScfDispersionGradient, D3BjForcesMatchFiniteDifference) {
    const D3BJ d3;
    const double worst = worst_fd_residual(d3, 1e-4);
    std::printf("[d3bj] worst |F_analytic - F_fd| = %.3e Ha/Bohr = %.3e eV/A\n",
                worst, worst * kHaBohrInEvAng);
    // 1e-9 Ha/Bohr = 5.1e-8 eV/Å — comfortably inside the 1e-6 eV/Å ask.
    EXPECT_LT(worst, 1e-9);
}

/// A gradient that is zero everywhere would pass the FD gate trivially. The
/// dispersion force on this dimer is a real, O(1e-4) Ha/Bohr force.
TEST(ScfDispersionGradient, D3BjForcesAreNotTriviallyZero) {
    const Frame f = water_dimer();
    EXPECT_GT(max_abs(D3BJ{}.forces(f)), 1e-5);
}

/// E_disp depends on interatomic distances only, so a rigid translation must
/// leave it alone: Σ_A F_A = 0 to machine precision.
TEST(ScfDispersionGradient, D3BjForcesSumToZero) {
    const Frame f = water_dimer();
    const auto forces = D3BJ{}.forces(f);
    for (int d = 0; d < 3; ++d) {
        double net = 0.0;
        for (const auto& t : forces)
            net += t[static_cast<std::size_t>(d)];
        EXPECT_NEAR(net, 0.0, 1e-14);
    }
}

/// `forces()` must read the same members `energy()` does — otherwise a tuned
/// s8 would move the surface but not the force on it.
TEST(ScfDispersionGradient, D3BjForcesFollowTheMembers) {
    D3BJ tuned;
    tuned.s8 = 1.5; // != the 0.7875 PBE default
    tuned.a1 = 0.30;
    const double worst = worst_fd_residual(tuned, 1e-4);
    EXPECT_LT(worst, 1e-9);
    const Frame f = water_dimer();
    EXPECT_GT(std::abs(max_abs(tuned.forces(f)) - max_abs(D3BJ{}.forces(f))),
              1e-7);
}

/// No fallback: an element outside the in-tree D3 table throws on the force
/// path exactly as it does on the energy path (a silently-zero dispersion
/// force is the defect this file exists to prevent).
TEST(ScfDispersionGradient, D3BjForcesThrowOnUnsupportedElement) {
    Frame f(1);
    f.span(kAtom_type)[0] = 100; // Fm — outside the table
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    EXPECT_THROW((void)D3BJ{}.forces(f), std::runtime_error);
}

// ── D4: analytic gradient of the D4 functional (two-body + EEQ adjoint) ───

/// D4's `forces()` is −∇ of D4's `energy()` (CN_d4 + EEQ charges + BJ pair
/// all move with the nuclei). Gate is analytic-vs-central-difference of the
/// same functional; dropping the EEQ adjoint or CN chain fails here.
TEST(ScfDispersionGradient, D4ForcesMatchFiniteDifference) {
    const D4 d4;
    const double worst = worst_fd_residual(d4, 1e-4);
    std::printf("[d4] worst |F_analytic - F_fd| = %.3e Ha/Bohr = %.3e eV/A\n",
                worst, worst * kHaBohrInEvAng);
    // Match D3(BJ) dual-gate level; offline dftd4 pin is 1e-7 Ha/Bohr.
    EXPECT_LT(worst, 1e-7);
}

/// Same non-triviality guard as D3(BJ).
TEST(ScfDispersionGradient, D4ForcesAreNotTriviallyZero) {
    const Frame f = water_dimer();
    EXPECT_GT(max_abs(D4{}.forces(f)), 1e-5);
}

/// Translation invariance of the D4 surface.
TEST(ScfDispersionGradient, D4ForcesSumToZero) {
    const Frame f = water_dimer();
    const auto forces = D4{}.forces(f);
    for (int d = 0; d < 3; ++d) {
        double net = 0.0;
        for (const auto& t : forces)
            net += t[static_cast<std::size_t>(d)];
        EXPECT_NEAR(net, 0.0, 1e-11);
    }
}

/// The ATM three-body term is part of the functional when `s9 != 0`, so it
/// must be part of the gradient too (FD of the full energy when ATM is on).
/// The first assertion keeps the row from going vacuous: `s9` has to actually
/// reach the energy for its gradient to mean anything.
TEST(ScfDispersionGradient, D4AtmForcesMatchFiniteDifference) {
    D4 d4;
    d4.s9 = 1.0; // PBE-D4-ATM
    const Frame f = water_dimer();
    EXPECT_GT(std::abs(d4.energy(f) - D4{}.energy(f)), 1e-8);
    const double worst = worst_fd_residual(d4, 5e-4);
    std::printf("[d4-atm] worst |F - F_fd| = %.3e Ha/Bohr\n", worst);
    // ATM path is still central-difference; gate step-independence vs 5e-4.
    EXPECT_LT(worst, 1e-9);
}

/// No fallback on the D4 force path either.
TEST(ScfDispersionGradient, D4ForcesThrowOnUnsupportedElement) {
    Frame f(1);
    f.span(kAtom_type)[0] = 95; // Am — beyond the in-tree D4 table (Z <= 94)
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    EXPECT_THROW((void)D4{}.forces(f), std::runtime_error);
}

// ── The shared flag → evaluator switch ────────────────────────────────────
//
// `internal/dispersion.h` is the single place that turns (d3bj, d4) into an
// evaluator. Four call sites depend on it agreeing with itself: the energy
// fold in `drive_scf_core_`, `ScfGradient`, `mp2_forces` and `ccsd_forces`.
// The defect these rows guard against is a force and an energy that answer the
// flags differently.

using ATV::cpu::scf::internal::add_dispersion_forces;
using ATV::cpu::scf::internal::dispersion_energy;
using ATV::cpu::scf::internal::dispersion_forces;

/// Flags off ⇒ exactly zero energy and an empty force list, so a vacuum run
/// is bit-identical to one that never asked about dispersion.
TEST(ScfDispersionSwitch, NoFlagsIsExactlyNothing) {
    const Frame f = water_dimer();
    EXPECT_EQ(dispersion_energy(f, false, false), 0.0);
    EXPECT_TRUE(dispersion_forces(f, false, false).empty());
    GradientResult gr;
    gr.forces = xt::zeros<double>({6, 3});
    gr.forces(0, 0) = 1.25;
    add_dispersion_forces(f, false, false, gr);
    EXPECT_EQ(gr.forces(0, 0), 1.25);
}

/// Energy and force must select the same flavour from the same flags — that
/// identity is the whole reason the switch is a shared function.
TEST(ScfDispersionSwitch, EnergyAndForceSelectTheSameFlavour) {
    const Frame f = water_dimer();
    EXPECT_EQ(dispersion_energy(f, true, false), D3BJ{}.energy(f));
    EXPECT_EQ(dispersion_forces(f, true, false)[0][2], D3BJ{}.forces(f)[0][2]);
    EXPECT_EQ(dispersion_energy(f, false, true), D4{}.energy(f));
    EXPECT_EQ(dispersion_forces(f, false, true)[0][2], D4{}.forces(f)[0][2]);
}

/// Both flags set: D4 wins on both paths (the method documents them as
/// mutually exclusive; the precedence exists so they cannot disagree).
TEST(ScfDispersionSwitch, D4WinsOverD3BjOnBothPaths) {
    const Frame f = water_dimer();
    EXPECT_EQ(dispersion_energy(f, true, true), D4{}.energy(f));
    EXPECT_EQ(dispersion_forces(f, true, true)[0][2], D4{}.forces(f)[0][2]);
}

/// `add_dispersion_forces` is the post-HF entry (mp2_forces / ccsd_forces):
/// it must add F = −∂E/∂R in the public eV/Å currency, not the a.u. it reads.
TEST(ScfDispersionSwitch, AddsForcesInEvPerAngstromWithTheRightSign) {
    const Frame f = water_dimer();
    GradientResult gr;
    gr.forces = xt::zeros<double>({6, 3});
    add_dispersion_forces(f, true, false, gr);
    const auto ha_bohr = D3BJ{}.forces(f);
    for (int a = 0; a < 6; ++a)
        for (int d = 0; d < 3; ++d)
            EXPECT_NEAR(gr.forces(a, d),
                        ha_bohr[static_cast<std::size_t>(a)]
                               [static_cast<std::size_t>(d)]
                            * kHaBohrInEvAng,
                        1e-14);
}

/// A gradient sized for a different system is a wiring bug, not something to
/// paper over with a partial add.
TEST(ScfDispersionSwitch, ThrowsOnAtomCountMismatch) {
    const Frame f = water_dimer();
    GradientResult gr;
    gr.forces = xt::zeros<double>({3, 3});
    EXPECT_THROW(add_dispersion_forces(f, true, false, gr),
                 std::runtime_error);
}

#endif // ATV_HAS_HF
