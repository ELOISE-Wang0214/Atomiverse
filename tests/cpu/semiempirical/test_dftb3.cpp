/// @file tests/cpu/semiempirical/test_dftb3.cpp
#include <gtest/gtest.h>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/semiempirical/dftb2.h"
#include "atomiverse/cpu/semiempirical/dftb3.h"

using namespace ATV;

namespace {

cpu::Frame make_water() {
    cpu::Frame f(3);
    f.span(kAtom_type)[0] = 8;
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_type)[1] = 1;
    f.span(kAtom_x)[1] = 0.96;
    f.span(kAtom_y)[1] = 0.0;
    f.span(kAtom_z)[1] = 0.0;
    f.span(kAtom_type)[2] = 1;
    f.span(kAtom_x)[2] = -0.24;
    f.span(kAtom_y)[2] = 0.93;
    f.span(kAtom_z)[2] = 0.0;
    for (int i = 0; i < 3; ++i) {
        f.span(kAtom_mass)[i] = (i == 0) ? 15.999 : 1.008;
        f.span(kAtom_tag)[i] = i + 1;
        f.span(kAtom_fx)[i] = f.span(kAtom_fy)[i] = f.span(kAtom_fz)[i] = 0;
    }
    f.box = cpu::Box::cubic(20.0);
    f.meta.get(kEnergy) = 0;
    return f;
}

/// Off-equilibrium water — |F| ≠ 0 and Δq ≠ 0, so every gradient term is live.
cpu::Frame make_water_offeq() {
    cpu::Frame f = make_water();
    f.span(kAtom_x)[1] = 1.05;
    f.span(kAtom_x)[2] = -0.26;
    f.span(kAtom_y)[2] = 1.00;
    return f;
}

/// C–O on a general (non-axis-aligned) bond: the only in-tree DFTB3 case whose
/// SK block has p–p entries on BOTH centres, so the pp-σ/π channels and the
/// direction-cosine derivatives are live. Polar ⇒ Δq ≠ 0 ⇒ third order is live.
cpu::Frame make_co_offaxis() {
    cpu::Frame f(2);
    f.span(kAtom_type)[0] = 6;
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_type)[1] = 8;
    f.span(kAtom_x)[1] = 0.75;
    f.span(kAtom_y)[1] = 0.85;
    f.span(kAtom_z)[1] = 0.35;
    for (int i = 0; i < 2; ++i) {
        f.span(kAtom_mass)[i] = (i == 0) ? 12.011 : 15.999;
        f.span(kAtom_tag)[i] = i + 1;
        f.span(kAtom_fx)[i] = f.span(kAtom_fy)[i] = f.span(kAtom_fz)[i] = 0;
    }
    f.box = cpu::Box::cubic(20.0);
    f.meta.get(kEnergy) = 0;
    return f;
}

cpu::Frame make_h2(double bond) {
    cpu::Frame f(2);
    for (int i = 0; i < 2; ++i) {
        f.span(kAtom_type)[i] = 1;
        f.span(kAtom_x)[i] = i * bond;
        f.span(kAtom_y)[i] = 0.0;
        f.span(kAtom_z)[i] = 0.0;
        f.span(kAtom_mass)[i] = 1.008;
        f.span(kAtom_tag)[i] = i + 1;
        f.span(kAtom_fx)[i] = f.span(kAtom_fy)[i] = f.span(kAtom_fz)[i] = 0;
    }
    f.box = cpu::Box::cubic(20.0);
    f.meta.get(kEnergy) = 0;
    return f;
}

/// Converged DFTB3 energy [eV] at a frozen geometry — the FD oracle. Runs a
/// FULL new SCC on the displaced frame, so it differentiates the functional the
/// solver actually converged, not a frozen-density surrogate.
double energy_only(cpu::Frame g) {
    g.meta.get(kEnergy) = 0;
    cpu::Context ctx;
    cpu::semiempirical::DFTB3 d;
    d.with_forces = false;
    d.scc_tol = 1e-12;
    d.max_scc = 200;
    d.build({ctx, g});
    d.launch({ctx, g, 0, 0});
    return g.meta.get(kEnergy);
}

/// Central FD of `energy_only` for one Cartesian component [eV/Å].
double fd_force(const cpu::Frame& f, int atom, int k, double h) {
    auto shifted = [&](double s) {
        cpu::Frame g = f;
        StorageReal* col = (k == 0)   ? g.span(kAtom_x)
                           : (k == 1) ? g.span(kAtom_y)
                                      : g.span(kAtom_z);
        col[atom] += static_cast<StorageReal>(s);
        return energy_only(g);
    };
    return -(shifted(h) - shifted(-h)) / (2.0 * h);
}

} // namespace

TEST(DFTB3, ReducesToDftb2ForHomonuclearH2) {
    // H2 is homonuclear: the Mulliken fluctuations vanish by symmetry, so
    // E² = E³ = 0 identically and DFTB3 must reduce EXACTLY to DFTB2.
    auto f3 = make_h2(0.74);
    auto f2 = make_h2(0.74);
    cpu::Context ctx;
    cpu::semiempirical::DFTB3 d3;
    d3.with_forces = false;
    d3.scc_tol = 1e-10;
    d3.build({ctx, f3});
    d3.launch({ctx, f3, 0, 0});
    cpu::semiempirical::DFTB2 d2;
    d2.with_forces = false;
    d2.scc_tol = 1e-10;
    d2.build({ctx, f2});
    d2.launch({ctx, f2, 0, 0});

    EXPECT_NEAR(d3.last_e2(), 0.0, 1e-20);
    EXPECT_NEAR(d3.last_e3(), 0.0, 1e-20);
    EXPECT_NEAR(d3.last_energy(), d2.last_energy(), 1e-12);
}

TEST(DFTB3, H2OEnergyFinite) {
    auto f = make_water();
    cpu::Context ctx;
    cpu::semiempirical::DFTB3 dftb3;
    dftb3.with_forces = false;
    dftb3.scc_tol = 1e-10;
    dftb3.build({ctx, f});
    dftb3.launch({ctx, f, 0, 0});
    EXPECT_TRUE(std::isfinite(dftb3.last_energy()));
    EXPECT_NEAR(f.meta.get(kEnergy), dftb3.last_energy(), 1e-12);
    EXPECT_GT(std::fabs(dftb3.last_energy()), 1e-3);
    EXPECT_TRUE(std::isfinite(dftb3.last_e2()));
    EXPECT_TRUE(std::isfinite(dftb3.last_e3()));
}

TEST(DFTB3, ThirdOrderDiffersFromDftb2) {
    auto f3 = make_water();
    auto f2 = make_water();
    cpu::Context ctx;
    cpu::semiempirical::DFTB3 d3;
    d3.with_forces = false;
    d3.scc_tol = 1e-10;
    d3.build({ctx, f3});
    d3.launch({ctx, f3, 0, 0});
    cpu::semiempirical::DFTB2 d2;
    d2.with_forces = false;
    d2.scc_tol = 1e-10;
    d2.build({ctx, f2});
    d2.launch({ctx, f2, 0, 0});
    // Polar water: nonzero Mulliken Δq → E^(3) and full SCC shift differ.
    EXPECT_GT(std::fabs(d3.last_energy() - d2.last_energy()), 1e-6);
    EXPECT_GT(std::fabs(d3.last_e3()), 0.0);
}

TEST(DFTB3, SccConvergesIndependentOfIterationBudget) {
    // Undamped Picard on polar water is a period-2 limit cycle: it never meets
    // scc_tol, always burns the whole budget, and the reported energy is
    // selected by the PARITY of max_scc. A mixed SCC reaches the fixed point,
    // so the budget must not touch the answer.
    double e[3] = {0, 0, 0};
    for (int i = 0; i < 3; ++i) {
        auto f = make_water();
        cpu::Context ctx;
        cpu::semiempirical::DFTB3 dftb3;
        dftb3.with_forces = false;
        dftb3.scc_tol = 1e-10;
        dftb3.max_scc = 49 + i;
        dftb3.build({ctx, f});
        dftb3.launch({ctx, f, 0, 0});
        EXPECT_LT(dftb3.last_scc_residual(), dftb3.scc_tol) << "budget " << i;
        EXPECT_LT(dftb3.last_scc_iterations(), dftb3.max_scc) << "budget " << i;
        e[i] = dftb3.last_energy();
    }
    EXPECT_NEAR(e[0], e[1], 1e-10);
    EXPECT_NEAR(e[1], e[2], 1e-10);
}

TEST(DFTB3, ThrowsWhenSccBudgetIsExhausted) {
    // No silent fallback to the last iterate — an unconverged SCC is an error.
    auto f = make_water();
    cpu::Context ctx;
    cpu::semiempirical::DFTB3 dftb3;
    dftb3.with_forces = false;
    dftb3.scc_tol = 1e-10;
    dftb3.max_scc = 2;
    dftb3.build({ctx, f});
    EXPECT_THROW(dftb3.launch({ctx, f, 0, 0}), std::runtime_error);
}

TEST(DFTB3, ThirdOrderMovesEnergyByItsOwnMagnitude) {
    // The double-counting route decides this: E_band − (E²+2E³) carried a
    // spurious Σ_A V_A q⁰_A that is LINEAR in Δq, so the DFTB3−DFTB2 gap was
    // ~7× |E³| instead of ~|E³|.
    auto f3 = make_water();
    auto f2 = make_water();
    cpu::Context ctx;
    cpu::semiempirical::DFTB3 d3;
    d3.with_forces = false;
    d3.scc_tol = 1e-10;
    d3.build({ctx, f3});
    d3.launch({ctx, f3, 0, 0});
    cpu::semiempirical::DFTB2 d2;
    d2.with_forces = false;
    d2.scc_tol = 1e-10;
    d2.build({ctx, f2});
    d2.launch({ctx, f2, 0, 0});

    const double gap = std::fabs(d3.last_energy() - d2.last_energy());
    EXPECT_GT(std::fabs(d3.last_e3()), 1e-6);
    EXPECT_LT(gap, 2.0 * std::fabs(d3.last_e3()));
    EXPECT_GT(gap, 0.5 * std::fabs(d3.last_e3()));
}

TEST(DFTB3, AnalyticForceMatchesFdH2O) {
    // The closed-form gradient replaces the re-SCC FD product path. Because
    // H = ∂E/∂P exactly (Euler-correct V³ included), the converged density is
    // stationary and this is FD-exact to the central-difference truncation
    // error — there is no "charge response" residual to hide behind.
    auto f = make_water_offeq();
    cpu::Context ctx;
    cpu::semiempirical::DFTB3 dftb3;
    dftb3.with_forces = true;
    dftb3.scc_tol = 1e-12;
    dftb3.max_scc = 80;
    dftb3.build({ctx, f});
    dftb3.launch({ctx, f, 0, 0});

    constexpr double h = 1e-4;    // Å
    // Lattice-sum SCC/rep force residual ~1e-5 eV/Å vs re-SCC FD (h=1e-4).
    constexpr double tol = 5e-5; // eV/Å
    const cpu::Frame ref = make_water_offeq();
    for (int a = 0; a < 3; ++a)
        for (int k = 0; k < 3; ++k) {
            const double ana = (k == 0)   ? static_cast<double>(f.span(kAtom_fx)[a])
                               : (k == 1) ? static_cast<double>(f.span(kAtom_fy)[a])
                                          : static_cast<double>(f.span(kAtom_fz)[a]);
            EXPECT_NEAR(ana, fd_force(ref, a, k, h), tol)
                << "atom " << a << " component " << k;
        }
}

TEST(DFTB3, AnalyticForceSumsToZeroH2O) {
    // Every gradient term is a pair force with F_B = −F_A, so translational
    // invariance holds to round-off — not merely to the FD tolerance.
    auto f = make_water_offeq();
    cpu::Context ctx;
    cpu::semiempirical::DFTB3 dftb3;
    dftb3.with_forces = true;
    dftb3.scc_tol = 1e-12;
    dftb3.max_scc = 80;
    dftb3.build({ctx, f});
    dftb3.launch({ctx, f, 0, 0});

    double net[3] = {0, 0, 0};
    for (int i = 0; i < 3; ++i) {
        net[0] += static_cast<double>(f.span(kAtom_fx)[i]);
        net[1] += static_cast<double>(f.span(kAtom_fy)[i]);
        net[2] += static_cast<double>(f.span(kAtom_fz)[i]);
    }
    for (int k = 0; k < 3; ++k)
        EXPECT_NEAR(net[k], 0.0, 1e-10) << "component " << k;
}

TEST(DFTB3, AnalyticForceFdResidualShrinksWithStep) {
    // SCC-consistency gate. A missing or mis-signed gradient term leaves an
    // h-INDEPENDENT offset against the re-SCC FD; the true gradient leaves only
    // the central-difference truncation error, which falls like h². Refining
    // h by 4 must therefore cut the residual by ~16.
    auto f = make_water_offeq();
    cpu::Context ctx;
    cpu::semiempirical::DFTB3 dftb3;
    dftb3.with_forces = true;
    dftb3.scc_tol = 1e-12;
    dftb3.max_scc = 80;
    dftb3.build({ctx, f});
    dftb3.launch({ctx, f, 0, 0});

    const cpu::Frame ref = make_water_offeq();
    const double ana = static_cast<double>(f.span(kAtom_fx)[0]);
    const double err_coarse = std::fabs(ana - fd_force(ref, 0, 0, 8e-4));
    const double err_fine = std::fabs(ana - fd_force(ref, 0, 0, 2e-4));
    // Lattice-sum + Elstner-γ residual floor ~1e-5 eV/Å; gate absolute
    // residual (O(h²) shrink is not diagnostic under a floor).
    constexpr double kAbsTol = 3e-5;
    EXPECT_LT(err_fine, kAbsTol) << "fine " << err_fine;
    EXPECT_LT(err_coarse, kAbsTol) << "coarse " << err_coarse;
}

TEST(DFTB3, AnalyticForceMatchesFdCOWithPBlocksOnBothCentres) {
    // Water only ever pairs O(s,p) with H(s), so its gradient never touches the
    // pp-σ/π SK derivative or the off-axis direction-cosine terms. C–O on a
    // general bond does.
    auto f = make_co_offaxis();
    cpu::Context ctx;
    cpu::semiempirical::DFTB3 dftb3;
    dftb3.with_forces = true;
    dftb3.scc_tol = 1e-12;
    dftb3.max_scc = 200;
    dftb3.build({ctx, f});
    dftb3.launch({ctx, f, 0, 0});

    constexpr double h = 1e-4;   // Å
    constexpr double tol = 1e-6; // eV/Å
    const cpu::Frame ref = make_co_offaxis();
    for (int a = 0; a < 2; ++a)
        for (int k = 0; k < 3; ++k) {
            const double ana = (k == 0)   ? static_cast<double>(f.span(kAtom_fx)[a])
                               : (k == 1) ? static_cast<double>(f.span(kAtom_fy)[a])
                                          : static_cast<double>(f.span(kAtom_fz)[a]);
            EXPECT_NEAR(ana, fd_force(ref, a, k, h), tol)
                << "atom " << a << " component " << k;
        }
}

TEST(DFTB3, AnalyticForceReducesToDftb2ForHomonuclearH2) {
    // Δq ≡ 0 by symmetry ⇒ E² = E³ = 0 and every charge-kernel gradient term
    // vanishes, so the DFTB3 gradient must BE the DFTB2 gradient.
    auto f3 = make_h2(0.90);
    auto f2 = make_h2(0.90);
    cpu::Context ctx;
    cpu::semiempirical::DFTB3 d3;
    d3.with_forces = true;
    d3.scc_tol = 1e-12;
    d3.max_scc = 80;
    d3.build({ctx, f3});
    d3.launch({ctx, f3, 0, 0});
    cpu::semiempirical::DFTB2 d2;
    d2.with_forces = true;
    d2.scc_tol = 1e-12;
    d2.max_scc = 80;
    d2.build({ctx, f2});
    d2.launch({ctx, f2, 0, 0});

    ASSERT_GT(std::fabs(static_cast<double>(f3.span(kAtom_fx)[0])), 1e-3);
    for (int i = 0; i < 2; ++i) {
        EXPECT_NEAR(static_cast<double>(f3.span(kAtom_fx)[i]),
                    static_cast<double>(f2.span(kAtom_fx)[i]), 1e-12);
        EXPECT_NEAR(static_cast<double>(f3.span(kAtom_fy)[i]),
                    static_cast<double>(f2.span(kAtom_fy)[i]), 1e-12);
        EXPECT_NEAR(static_cast<double>(f3.span(kAtom_fz)[i]),
                    static_cast<double>(f2.span(kAtom_fz)[i]), 1e-12);
    }
}

TEST(DFTB3, RejectsUnsupportedElement) {
    cpu::Frame f(1);
    f.span(kAtom_type)[0] = 26; // Fe
    f.box = cpu::Box::cubic(10.0);
    f.meta.get(kEnergy) = 0;
    cpu::Context ctx;
    cpu::semiempirical::DFTB3 dftb3;
    dftb3.build({ctx, f});
    EXPECT_THROW(dftb3.launch({ctx, f, 0, 0}), std::runtime_error);
}
