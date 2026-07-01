/// @file tests/cpu/scf/test_gradient_vs_libcint.cpp
/// @brief Regression net: native 2e gradient must match libcint to tight tol.
///
/// libcint is treated as ground truth (same engine PySCF uses via cint2e_ip1).
/// When `ATV_HAS_LIBCINT` is defined these tests run the SCF twice-for-once:
///   pass 1: libcint is active — capture the reference force field
///   pass 2: `SCFConfig::force_native_two_electron_deriv = true` forces the
///           native Rys-derivative path — capture the same molecule's force
///           and compare element-wise.
///
/// Tolerances (per component, eV/Å):
///   H2O  STO-3G (7 AOs)  — 1e-5
///   CH4  STO-3G (9 AOs)  — 1e-5
///   aspirin STO-3G (73 AOs) — 1e-4 (larger basis accumulates FP error)
/// These match the documented Rys-vs-libcint ERI baseline (~5e-7 a.u.)
/// propagated through the gradient contraction — absolute limit of the
/// native Rys engine's agreement with libcint; structural bugs would
/// show up at magnitudes far above these.
///
/// Status: Tests are DISABLED initially because the current legacy native
/// path is not only slow but also numerically incorrect on p-type shells
/// (CH4 diverges from libcint by ~10% — never caught because existing
/// gradient tests cover only H2/HeH which have no p-shell derivatives).
/// Phase 2 Stage 3 replaces the legacy derivative engine with a Rys-based
/// one that matches libcint to FP precision; at that point the DISABLED_
/// prefix is removed and these tests become the permanent regression net.
///
/// Without `ATV_HAS_LIBCINT`, the file compiles to a no-op — the reference
/// backend is unavailable so there is nothing to compare against.
#if defined(ATV_HAS_RHF) && defined(ATV_HAS_LIBCINT)

#include <cmath>
#include <gtest/gtest.h>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/gradient.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_solve.h"

namespace {

using namespace ATV;

// ── Molecule builders ────────────────────────────────────────────────

cpu::Frame make_h2o() {
    cpu::Frame f(3);
    // Near-equilibrium geometry (Angstrom)
    f.atoms.x[0] = 0.0;
    f.atoms.y[0] = 0.0;
    f.atoms.z[0] = 0.117;
    f.atoms.x[1] = 0.0;
    f.atoms.y[1] = 0.757;
    f.atoms.z[1] = -0.468;
    f.atoms.x[2] = 0.0;
    f.atoms.y[2] = -0.757;
    f.atoms.z[2] = -0.468;
    f.atoms.type[0] = 8;
    f.atoms.mass[0] = 15.999;
    f.atoms.type[1] = 1;
    f.atoms.mass[1] = 1.008;
    f.atoms.type[2] = 1;
    f.atoms.mass[2] = 1.008;
    f.box = Mat3::cubic(100.0);
    return f;
}

cpu::Frame make_ch4() {
    cpu::Frame f(5);
    // Tetrahedral CH4 (Angstrom)
    f.atoms.x[0] = 0.0;
    f.atoms.y[0] = 0.0;
    f.atoms.z[0] = 0.0;
    f.atoms.x[1] = 0.629;
    f.atoms.y[1] = 0.629;
    f.atoms.z[1] = 0.629;
    f.atoms.x[2] = -0.629;
    f.atoms.y[2] = -0.629;
    f.atoms.z[2] = 0.629;
    f.atoms.x[3] = -0.629;
    f.atoms.y[3] = 0.629;
    f.atoms.z[3] = -0.629;
    f.atoms.x[4] = 0.629;
    f.atoms.y[4] = -0.629;
    f.atoms.z[4] = -0.629;
    f.atoms.type[0] = 6;
    f.atoms.mass[0] = 12.011;
    for (int i = 1; i < 5; ++i) {
        f.atoms.type[i] = 1;
        f.atoms.mass[i] = 1.008;
    }
    f.box = Mat3::cubic(100.0);
    return f;
}

/// @brief Aspirin C9H8O4 (21 atoms, 73 AOs in STO-3G).
/// Same geometry as examples/tutorial/05_geom_opt.cpp — mildly distorted
/// PubChem conformer so the gradient is non-zero.
cpu::Frame make_aspirin() {
    cpu::Frame f(21);
    const int types[21] = {
        6, 6, 6, 6, 6, 6, 6, 8, 8, 8, 6, 8, 6, 1, 1, 1, 1, 1, 1, 1, 1,
    };
    const double coords[21][3] = {
        {1.2879, 0.1068, 0.0000},    {0.1768, 1.0827, 0.0000},
        {0.5171, 2.4290, 0.0000},    {1.8519, 2.7993, 0.0000},
        {2.9293, 1.8524, 0.0000},    {2.6191, 0.4890, 0.0000},
        {-1.2076, 0.6698, 0.0000},   {-1.5541, -0.4994, 0.0000},
        {-2.0632, 1.6805, 0.0000},   {0.9973, -1.2285, 0.0000},
        {-0.2830, -1.7447, 0.0000},  {-0.3696, -2.9537, 0.0000},
        {-1.3593, -0.7411, 0.3000},  {-0.2980, 3.1365, 0.0000},
        {2.0791, 3.8586, 0.0000},    {3.9560, 2.1939, 0.0000},
        {3.4546, -0.2019, 0.0000},   {-2.9828, 1.3621, 0.0000},
        {-1.0697, -0.1755, 1.1920},  {-2.3482, -1.1583, 0.5001},
        {-1.4267, -0.2091, -0.6502},
    };
    for (int i = 0; i < 21; ++i) {
        f.atoms.type[i] = types[i];
        f.atoms.mass[i] = (types[i] == 6)   ? 12.011
                          : (types[i] == 8) ? 15.999
                                            : 1.008;
        f.atoms.x[i] = coords[i][0];
        f.atoms.y[i] = coords[i][1];
        f.atoms.z[i] = coords[i][2];
    }
    f.box = Mat3::cubic(100.0);
    return f;
}

// ── Shared machinery ─────────────────────────────────────────────────

struct GradCompare {
    double max_abs_diff = 0.0;
    double max_rel_diff = 0.0;
    int worst_atom = -1;
    int worst_axis = -1;
    double ref_value_at_worst = 0.0;
    double test_value_at_worst = 0.0;
};

GradCompare diff_forces(const xt::xtensor<double, 2>& ref,
                        const xt::xtensor<double, 2>& test) {
    GradCompare r;
    const int natm = static_cast<int>(ref.shape()[0]);
    for (int a = 0; a < natm; ++a) {
        for (int d = 0; d < 3; ++d) {
            const double rv = ref(a, d);
            const double tv = test(a, d);
            const double abs_diff = std::abs(rv - tv);
            const double denom = std::max(std::abs(rv), 1.0);
            const double rel_diff = abs_diff / denom;
            if (abs_diff > r.max_abs_diff) {
                r.max_abs_diff = abs_diff;
                r.max_rel_diff = rel_diff;
                r.worst_atom = a;
                r.worst_axis = d;
                r.ref_value_at_worst = rv;
                r.test_value_at_worst = tv;
            }
        }
    }
    return r;
}

/// Run `scf_solve` + `scf_gradient` on `frame`, returning the force field.
/// When `force_native` is true, routes through the native 2e-derivative path
/// even if libcint is compiled in.
xt::xtensor<double, 2> compute_gradient(const cpu::Frame& frame,
                                        bool force_native) {
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G,
                                               /*charge=*/0);
    auto config = cpu::scf::rhf_config(integrals);
    config.tol = 1e-10; // tight SCF so the gradient has no SCF noise
    config.max_iter = 200;
    config.force_native_two_electron_deriv = force_native;

    auto result = cpu::scf::scf_solve(integrals, config);
    EXPECT_TRUE(result.converged) << "SCF did not converge";

    auto grad = cpu::scf::scf_gradient(frame, integrals, result, config);
    return grad.forces;
}

// ── Tests ────────────────────────────────────────────────────────────

TEST(GradientNativeVsLibcint, H2O_ToTenMicroHartree) {
    auto frame = make_h2o();
    auto f_ref = compute_gradient(frame, /*force_native=*/false);
    auto f_test = compute_gradient(frame, /*force_native=*/true);
    auto cmp = diff_forces(f_ref, f_test);
    EXPECT_LT(cmp.max_abs_diff, 1e-5)
        << "native vs libcint diverge at atom " << cmp.worst_atom << " axis "
        << cmp.worst_axis << " ref=" << cmp.ref_value_at_worst
        << " test=" << cmp.test_value_at_worst;
}

TEST(GradientNativeVsLibcint, CH4_ToTenMicroHartree) {
    auto frame = make_ch4();
    auto f_ref = compute_gradient(frame, /*force_native=*/false);
    auto f_test = compute_gradient(frame, /*force_native=*/true);
    auto cmp = diff_forces(f_ref, f_test);
    EXPECT_LT(cmp.max_abs_diff, 1e-5)
        << "native vs libcint diverge at atom " << cmp.worst_atom << " axis "
        << cmp.worst_axis << " ref=" << cmp.ref_value_at_worst
        << " test=" << cmp.test_value_at_worst;
}

// Enabled after Stage 3: the Rys-based native gradient now runs
// aspirin in seconds (vs. 15 min for the legacy path), and its
// numerical agreement with libcint is within the Rys ERI baseline
// propagated through the gradient contraction.
TEST(GradientNativeVsLibcint, Aspirin_ToHundredMicroHartree) {
    auto frame = make_aspirin();
    auto f_ref = compute_gradient(frame, /*force_native=*/false);
    auto f_test = compute_gradient(frame, /*force_native=*/true);
    auto cmp = diff_forces(f_ref, f_test);
    EXPECT_LT(cmp.max_abs_diff, 1e-4)
        << "native vs libcint diverge at atom " << cmp.worst_atom << " axis "
        << cmp.worst_axis << " ref=" << cmp.ref_value_at_worst
        << " test=" << cmp.test_value_at_worst;
}

} // namespace

#endif // ATV_HAS_RHF && ATV_HAS_LIBCINT
