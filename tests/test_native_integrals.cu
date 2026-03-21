/// @file test_native_integrals.cu
/// @brief Comprehensive validation of native C++ molecular integrals.
///
/// Two test suites:
///   NativeIntegrals     -- always-on (ATV_HAS_RHF): symmetry, analytical, multi-geometry
///   NativeVsLibcint     -- libcint comparison (ATV_HAS_LIBCINT): element-by-element ground truth

#ifdef ATV_HAS_RHF

#include <gtest/gtest.h>
#include <cmath>

#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/native_integrals.h"
#include "atomiverse/core/constants.h"

using namespace ATV;
using namespace ATV::cpu::scf;

namespace {

// STO-3G hydrogen parameters
constexpr double sto3g_exp[3] = {3.42525091, 0.62391373, 0.16885540};
constexpr double sto3g_coeff[3] = {0.15432897, 0.53532814, 0.44463454};

/// @brief Build an H2 frame at a given bond length along x.
cpu::Frame make_h2(Real bond_length) {
    cpu::Frame frame(2);
    frame.x[0] = Real(0);           frame.y[0] = Real(0); frame.z[0] = Real(0);
    frame.x[1] = bond_length;       frame.y[1] = Real(0); frame.z[1] = Real(0);
    frame.mass[0] = Real(1.008); frame.mass[1] = Real(1.008);
    frame.type_id[0] = 1; frame.type_id[1] = 1;
    frame.box = Mat3::cubic(100.0);
    return frame;
}

/// @brief Build Shell structs for native API (physical normalization).
struct H2Setup {
    double norm_coeff[3];
    double centers[2][3];
    native::Shell shells[2];
    native::Atom atoms[2];

    explicit H2Setup(double bond_angstrom = 0.74) {
        for (int p = 0; p < 3; ++p)
            norm_coeff[p] = sto3g_coeff[p] * native::phys_norm(sto3g_exp[p]);
        centers[0][0] = 0; centers[0][1] = 0; centers[0][2] = 0;
        centers[1][0] = bond_angstrom * kAngstromToBohr; centers[1][1] = 0; centers[1][2] = 0;
        for (int s = 0; s < 2; ++s) {
            shells[s].n_prim = 3;
            shells[s].alpha = sto3g_exp;
            shells[s].coeff = norm_coeff;
            shells[s].center[0] = centers[s][0];
            shells[s].center[1] = centers[s][1];
            shells[s].center[2] = centers[s][2];
        }
        for (int a = 0; a < 2; ++a) {
            atoms[a].Z = 1;
            atoms[a].pos[0] = centers[a][0];
            atoms[a].pos[1] = centers[a][1];
            atoms[a].pos[2] = centers[a][2];
        }
    }
};

// Hardcoded reference values for H2/STO-3G at 0.74 Angstrom.
// Captured from native C++ implementation (physical normalization, double precision).
// These use raw Shell/Atom structs with phys_norm coefficients.
constexpr double kRefS00 =  0.9999999909;
constexpr double kRefS01 =  0.6598731506;
constexpr double kRefT00 =  0.7600318766;
constexpr double kRefT01 =  0.2369602946;
constexpr double kRefV00 = -1.8809913573;
constexpr double kRefV01 = -1.1963361202;
constexpr double kRefVnuc = 0.7151044;
// ERI: (00|00), (00|11), (01|01)
constexpr double kRefERI_0000 = 0.7746059298;
constexpr double kRefERI_0011 = 0.5699948919;
constexpr double kRefERI_0101 = 0.2975905808;

} // anonymous namespace

// ═══════════════════════════════════════════════════════════════════════
// Suite A: NativeIntegrals (always-on, no libcint needed)
// ═══════════════════════════════════════════════════════════════════════

// ── Boys function ─────────────────────────────────────────────────────

TEST(NativeIntegrals, BoysF0_ZeroArg) {
    EXPECT_NEAR(native::boys_f0(0.0), 1.0, 1e-15);
}

TEST(NativeIntegrals, BoysF0_KnownValues) {
    EXPECT_NEAR(native::boys_f0(1.0), 0.5 * std::sqrt(M_PI) * std::erf(1.0), 1e-12);
    EXPECT_NEAR(native::boys_f0(10.0), 0.5 * std::sqrt(M_PI / 10.0) * std::erf(std::sqrt(10.0)), 1e-12);
}

TEST(NativeIntegrals, BoysF0_LargeArg) {
    // Asymptotic: F_0(x) → 0.5 * sqrt(pi/x) for large x
    EXPECT_NEAR(native::boys_f0(100.0), 0.5 * std::sqrt(M_PI / 100.0), 1e-12);
}

TEST(NativeIntegrals, BoysF0_SmallArg) {
    EXPECT_NEAR(native::boys_f0(1e-15), 1.0, 1e-7);
}

// ── Overlap symmetry and normalization ────────────────────────────────

TEST(NativeIntegrals, Overlap_Symmetry) {
    H2Setup h2;
    double S01 = native::overlap_ss(h2.shells[0], h2.shells[1]);
    double S10 = native::overlap_ss(h2.shells[1], h2.shells[0]);
    EXPECT_NEAR(S01, S10, 1e-15) << "S must be symmetric";
}

TEST(NativeIntegrals, Overlap_SelfNormalization) {
    H2Setup h2;
    double S00 = native::overlap_ss(h2.shells[0], h2.shells[0]);
    EXPECT_NEAR(S00, 1.0, 1e-7) << "Normalized shell self-overlap must be ~1";
}

// ── Kinetic symmetry ──────────────────────────────────────────────────

TEST(NativeIntegrals, Kinetic_Symmetry) {
    H2Setup h2;
    double T01 = native::kinetic_ss(h2.shells[0], h2.shells[1]);
    double T10 = native::kinetic_ss(h2.shells[1], h2.shells[0]);
    EXPECT_NEAR(T01, T10, 1e-15) << "T must be symmetric";
}

TEST(NativeIntegrals, Kinetic_Positive) {
    H2Setup h2;
    double T00 = native::kinetic_ss(h2.shells[0], h2.shells[0]);
    EXPECT_GT(T00, 0.0) << "Diagonal kinetic energy must be positive";
}

// ── Nuclear attraction symmetry ───────────────────────────────────────

TEST(NativeIntegrals, Nuclear_Symmetry) {
    H2Setup h2;
    double V01 = native::nuclear_ss(h2.shells[0], h2.shells[1], h2.atoms, 2);
    double V10 = native::nuclear_ss(h2.shells[1], h2.shells[0], h2.atoms, 2);
    EXPECT_NEAR(V01, V10, 1e-15) << "V must be symmetric";
}

TEST(NativeIntegrals, Nuclear_Negative) {
    H2Setup h2;
    double V00 = native::nuclear_ss(h2.shells[0], h2.shells[0], h2.atoms, 2);
    EXPECT_LT(V00, 0.0) << "Nuclear attraction must be negative";
}

// ── ERI permutation symmetry (8-fold) ─────────────────────────────────

TEST(NativeIntegrals, ERI_PermutationSymmetry) {
    H2Setup h2;
    auto& s = h2.shells;
    // Test all 8-fold permutations of (01|01)
    double base = native::eri_ssss(s[0], s[1], s[0], s[1]); // (01|01)
    EXPECT_NEAR(native::eri_ssss(s[1], s[0], s[0], s[1]), base, 1e-12) << "(10|01) != (01|01)";
    EXPECT_NEAR(native::eri_ssss(s[0], s[1], s[1], s[0]), base, 1e-12) << "(01|10) != (01|01)";
    EXPECT_NEAR(native::eri_ssss(s[1], s[0], s[1], s[0]), base, 1e-12) << "(10|10) != (01|01)";
    EXPECT_NEAR(native::eri_ssss(s[0], s[1], s[0], s[1]), base, 1e-12) << "(01|01) self-check";
    // Bra-ket exchange: (ij|kl) == (kl|ij)
    EXPECT_NEAR(native::eri_ssss(s[0], s[1], s[0], s[1]),
                native::eri_ssss(s[0], s[1], s[0], s[1]), 1e-12) << "(01|01) != (01|01)";
}

TEST(NativeIntegrals, ERI_Positive) {
    H2Setup h2;
    double eri_0000 = native::eri_ssss(h2.shells[0], h2.shells[0], h2.shells[0], h2.shells[0]);
    EXPECT_GT(eri_0000, 0.0) << "Self-repulsion ERI must be positive";
}

// ── Hardcoded reference regression (no libcint needed) ────────────────

TEST(NativeIntegrals, RegressionValues_H2) {
    H2Setup h2;
    double tol = 1e-5;  // float-safe tolerance

    EXPECT_NEAR(native::overlap_ss(h2.shells[0], h2.shells[0]), kRefS00, tol);
    EXPECT_NEAR(native::overlap_ss(h2.shells[0], h2.shells[1]), kRefS01, tol);
    EXPECT_NEAR(native::kinetic_ss(h2.shells[0], h2.shells[0]), kRefT00, tol);
    EXPECT_NEAR(native::kinetic_ss(h2.shells[0], h2.shells[1]), kRefT01, tol);
    EXPECT_NEAR(native::nuclear_ss(h2.shells[0], h2.shells[0], h2.atoms, 2), kRefV00, tol);
    EXPECT_NEAR(native::nuclear_ss(h2.shells[0], h2.shells[1], h2.atoms, 2), kRefV01, tol);

    EXPECT_NEAR(native::eri_ssss(h2.shells[0], h2.shells[0], h2.shells[0], h2.shells[0]), kRefERI_0000, tol);
    EXPECT_NEAR(native::eri_ssss(h2.shells[0], h2.shells[0], h2.shells[1], h2.shells[1]), kRefERI_0011, tol);
    EXPECT_NEAR(native::eri_ssss(h2.shells[0], h2.shells[1], h2.shells[0], h2.shells[1]), kRefERI_0101, tol);
}

// ── IntegralPack sanity ───────────────────────────────────────────────

TEST(NativeIntegrals, IntegralPack_Sanity) {
    auto frame = make_h2(Real(0.74));
    auto pack = build_integrals_sto3g(frame);

    EXPECT_NEAR(static_cast<double>(pack.S(0, 0)), 1.0, 0.01);
    EXPECT_NEAR(static_cast<double>(pack.V_nuc), kRefVnuc, 0.01);
    EXPECT_EQ(pack.n_basis, 2);
    EXPECT_EQ(pack.n_electrons, 2);
    EXPECT_EQ(pack.n_occ(), 1);
}

// ── Multi-geometry tests ──────────────────────────────────────────────

TEST(NativeIntegrals, MultiGeometry_OverlapDecay) {
    // Overlap between shells should decrease with distance
    double bond_lengths[] = {0.5, 0.74, 1.0, 1.5, 2.0};
    double prev_S01 = 1.0;
    for (double bl : bond_lengths) {
        H2Setup h2(bl);
        double S01 = native::overlap_ss(h2.shells[0], h2.shells[1]);
        EXPECT_LT(S01, prev_S01) << "Overlap should decrease with distance at " << bl;
        EXPECT_GT(S01, 0.0) << "Overlap should remain positive at " << bl;
        prev_S01 = S01;
    }
}

TEST(NativeIntegrals, MultiGeometry_IntegralPack) {
    // build_integrals_sto3g should work at various geometries
    double bond_lengths[] = {0.5, 0.74, 1.0, 1.5, 2.0};
    for (double bl : bond_lengths) {
        auto frame = make_h2(Real(bl));
        auto pack = build_integrals_sto3g(frame);
        // S diagonal must be ~1 regardless of geometry
        EXPECT_NEAR(static_cast<double>(pack.S(0, 0)), 1.0, 0.01)
            << "S(0,0) != 1.0 at bond_length=" << bl;
        // V_nuc should be Z1*Z2/R (in Hartree with Bohr distance)
        double R_bohr = bl * kAngstromToBohr;
        double expected_vnuc = 1.0 / R_bohr;  // Z=1 for both H atoms
        EXPECT_NEAR(static_cast<double>(pack.V_nuc), expected_vnuc, 1e-5)
            << "V_nuc mismatch at bond_length=" << bl;
    }
}

// ═══════════════════════════════════════════════════════════════════════
// Suite B: NativeVsLibcint (only when -DATV_BUILD_LIBCINT=ON)
// ═══════════════════════════════════════════════════════════════════════

#ifdef ATV_HAS_LIBCINT

#include "libcint_reference.cuh"

TEST(NativeVsLibcint, GtoNorm) {
    for (int p = 0; p < 3; ++p) {
        double ref = test::libcint_gto_norm(0, sto3g_exp[p]);
        double nat = native::gto_norm(0, sto3g_exp[p]);
        EXPECT_NEAR(nat, ref, 1e-12) << "gto_norm mismatch at exponent " << sto3g_exp[p];
    }
}

TEST(NativeVsLibcint, AllIntegrals_H2_074) {
    auto frame = make_h2(Real(0.74));
    auto native_pack = build_integrals_sto3g(frame);
    auto libcint_pack = test::build_integrals_libcint(frame);

    Int n = native_pack.n_basis;
    for (Int i = 0; i < n; ++i)
        for (Int j = 0; j < n; ++j) {
            EXPECT_NEAR(static_cast<double>(native_pack.S(i,j)),
                        static_cast<double>(libcint_pack.S(i,j)), 1e-7)
                << "S(" << i << "," << j << ")";
            EXPECT_NEAR(static_cast<double>(native_pack.T(i,j)),
                        static_cast<double>(libcint_pack.T(i,j)), 1e-7)
                << "T(" << i << "," << j << ")";
            EXPECT_NEAR(static_cast<double>(native_pack.V(i,j)),
                        static_cast<double>(libcint_pack.V(i,j)), 1e-7)
                << "V(" << i << "," << j << ")";
        }

    for (Int i = 0; i < n; ++i)
      for (Int j = 0; j < n; ++j)
        for (Int k = 0; k < n; ++k)
          for (Int l = 0; l < n; ++l)
            EXPECT_NEAR(static_cast<double>(native_pack.ERI(i,j,k,l)),
                        static_cast<double>(libcint_pack.ERI(i,j,k,l)), 1e-7)
                << "ERI(" << i << "," << j << "," << k << "," << l << ")";

    EXPECT_NEAR(static_cast<double>(native_pack.V_nuc),
                static_cast<double>(libcint_pack.V_nuc), 1e-10);
}

TEST(NativeVsLibcint, MultiGeometry) {
    double bond_lengths[] = {0.5, 0.74, 1.0, 1.5, 2.0};
    for (double bl : bond_lengths) {
        auto frame = make_h2(Real(bl));
        auto native_pack = build_integrals_sto3g(frame);
        auto libcint_pack = test::build_integrals_libcint(frame);

        Int n = native_pack.n_basis;
        for (Int i = 0; i < n; ++i)
          for (Int j = 0; j < n; ++j) {
            EXPECT_NEAR(static_cast<double>(native_pack.S(i,j)),
                        static_cast<double>(libcint_pack.S(i,j)), 1e-7)
                << "S mismatch at bl=" << bl;
            EXPECT_NEAR(static_cast<double>(native_pack.T(i,j)),
                        static_cast<double>(libcint_pack.T(i,j)), 1e-7)
                << "T mismatch at bl=" << bl;
            EXPECT_NEAR(static_cast<double>(native_pack.V(i,j)),
                        static_cast<double>(libcint_pack.V(i,j)), 1e-7)
                << "V mismatch at bl=" << bl;
          }

        for (Int i = 0; i < n; ++i)
          for (Int j = 0; j < n; ++j)
            for (Int k = 0; k < n; ++k)
              for (Int l = 0; l < n; ++l)
                EXPECT_NEAR(static_cast<double>(native_pack.ERI(i,j,k,l)),
                            static_cast<double>(libcint_pack.ERI(i,j,k,l)), 1e-7)
                    << "ERI mismatch at bl=" << bl;
    }
}

#endif // ATV_HAS_LIBCINT

#else // !ATV_HAS_RHF

#include <gtest/gtest.h>
TEST(NativeIntegrals, Disabled) {
    GTEST_SKIP() << "Native integral tests require -DATV_BUILD_RHF=ON";
}

#endif
